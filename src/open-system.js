// Dinámica de un sistema cuántico abierto de pocos niveles.
//  · Ecuación maestra de Lindblad: dρ/dt = −i[H,ρ] + Σₖ γₖ·(LₖρLₖ† − ½{Lₖ†Lₖ, ρ}), con ħ = 1.
//    Se resuelve de forma exacta para H y Lₖ constantes: ρ(t+Δt) = exp(𝓛·Δt)·ρ(t), con 𝓛 el superoperador.
//  · Canales de Kraus: ρ → Σₖ KₖρKₖ†.
//  · Trayectorias de saltos cuánticos: realizaciones puras cuyo promedio reproduce la ecuación maestra.

import { zeros, identity, multiply, dagger, conjugate, add, scale, kron, expm, applyToVector, distance } from './cmatrix.js';

// Superoperador de Lindblad sobre ρ escrita por filas: vec(AρB) = (A ⊗ Bᵀ)·vec(ρ).
export function liouvillian(hamiltonian, collapse = []) {
  const n = hamiltonian.n;
  const one = identity(n);
  const transpose = (a) => conjugate(dagger(a));
  // −i·(H ⊗ I − I ⊗ Hᵀ)
  let result = scale(add(kron(hamiltonian, one), scale(kron(one, transpose(hamiltonian)), -1)), 0, -1);
  for (const { operator, rate } of collapse) {
    if (!(rate > 0)) continue;
    const product = multiply(dagger(operator), operator);
    const jump = kron(operator, conjugate(operator));
    const loss = scale(add(kron(product, one), kron(one, transpose(product))), -0.5);
    result = add(result, scale(add(jump, loss), rate));
  }
  return result;
}

// Propagador exacto de un paso Δt para un superoperador constante.
export function superPropagator(superoperator, dt) {
  return expm(scale(superoperator, dt));
}

// Aplica un superoperador (o su propagador) a una matriz densidad.
export function applySuper(propagator, rho) {
  const n = rho.n;
  const moved = applyToVector(propagator, { re: rho.re, im: rho.im });
  return { n, re: moved.re, im: moved.im };
}

export function applyKraus(rho, operators) {
  let result = zeros(rho.n);
  for (const operator of operators) result = add(result, multiply(multiply(operator, rho), dagger(operator)));
  return result;
}

// Desviación de Σ Kₖ†Kₖ respecto de la identidad: cero si el canal conserva la traza.
export function krausDefect(operators) {
  let total = zeros(operators[0].n);
  for (const operator of operators) total = add(total, multiply(dagger(operator), operator));
  return distance(total, identity(operators[0].n));
}

// Canales de un qubit. «strength» es la probabilidad del suceso en cada uso del canal.
export const channel = Object.freeze({
  // Pérdida de fase: las coherencias se multiplican por √(1 − λ) y las poblaciones no cambian.
  dephasing(strength) {
    const keep = Math.sqrt(1 - strength);
    return [{ n: 2, re: Float64Array.of(1, 0, 0, keep), im: new Float64Array(4) }, { n: 2, re: Float64Array.of(0, 0, 0, Math.sqrt(strength)), im: new Float64Array(4) }];
  },
  // Relajación: |1⟩ cae a |0⟩ con probabilidad γ.
  amplitudeDamping(strength) {
    return [{ n: 2, re: Float64Array.of(1, 0, 0, Math.sqrt(1 - strength)), im: new Float64Array(4) }, { n: 2, re: Float64Array.of(0, Math.sqrt(strength), 0, 0), im: new Float64Array(4) }];
  },
  // Despolarización: el vector de Bloch se encoge por igual en todas las direcciones.
  depolarizing(strength) {
    const side = Math.sqrt(strength / 4);
    return [
      { n: 2, re: Float64Array.of(Math.sqrt(1 - (3 * strength) / 4), 0, 0, Math.sqrt(1 - (3 * strength) / 4)), im: new Float64Array(4) },
      { n: 2, re: Float64Array.of(0, side, side, 0), im: new Float64Array(4) },
      { n: 2, re: new Float64Array(4), im: Float64Array.of(0, -side, side, 0) },
      { n: 2, re: Float64Array.of(side, 0, 0, -side), im: new Float64Array(4) }
    ];
  }
});

// Trayectorias de saltos cuánticos por el método del tiempo de espera: entre saltos el ket evoluciona
// con H_ef = H − (i/2)·Σ γₖLₖ†Lₖ, cuya norma decrece; cuando cae por debajo de un umbral sorteado, salta.
export function createJumpUnraveling(hamiltonian, collapse, dt) {
  const active = collapse.filter((item) => item.rate > 0);
  let effective = scale(hamiltonian, 0, -1);
  for (const { operator, rate } of active) effective = add(effective, scale(multiply(dagger(operator), operator), -0.5 * rate));
  const propagator = expm(scale(effective, dt));
  const norm = (psi) => {
    let total = 0;
    for (let index = 0; index < psi.re.length; index += 1) total += psi.re[index] ** 2 + psi.im[index] ** 2;
    return total;
  };
  const rescale = (psi, factor) => {
    for (let index = 0; index < psi.re.length; index += 1) {
      psi.re[index] *= factor;
      psi.im[index] *= factor;
    }
    return psi;
  };
  return {
    // Estado de una trayectoria: el ket sin normalizar desde el último salto y el umbral sorteado.
    start(psi, rng) {
      return { psi: { re: Float64Array.from(psi.re), im: Float64Array.from(psi.im) }, threshold: rng.uniform(), jumps: 0, last: -1 };
    },
    // Avanza un paso. Devuelve el índice del salto ocurrido o −1.
    step(trajectory, rng) {
      trajectory.psi = applyToVector(propagator, trajectory.psi);
      trajectory.last = -1;
      if (norm(trajectory.psi) > trajectory.threshold || active.length === 0) return -1;
      const candidates = active.map(({ operator }) => applyToVector(operator, trajectory.psi));
      const weights = candidates.map((candidate, index) => active[index].rate * norm(candidate));
      const total = weights.reduce((sum, value) => sum + value, 0);
      if (!(total > 1e-300)) {
        // Ningún canal puede actuar sobre este ket: se renormaliza y se sigue sin salto.
        trajectory.psi = rescale(trajectory.psi, 1 / Math.sqrt(norm(trajectory.psi)));
        trajectory.threshold = rng.uniform();
        return -1;
      }
      let pick = rng.uniform() * total;
      let chosen = 0;
      for (let index = 0; index < weights.length; index += 1) {
        chosen = index;
        pick -= weights[index];
        if (pick <= 0) break;
      }
      trajectory.psi = rescale(candidates[chosen], 1 / Math.sqrt(norm(candidates[chosen])));
      trajectory.threshold = rng.uniform();
      trajectory.jumps += 1;
      trajectory.last = chosen;
      return chosen;
    },
    // Ket normalizado de la trayectoria en este instante.
    state(trajectory) {
      const copy = { re: Float64Array.from(trajectory.psi.re), im: Float64Array.from(trajectory.psi.im) };
      return rescale(copy, 1 / Math.sqrt(norm(copy)));
    }
  };
}
