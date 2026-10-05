// Estados de sistemas de pocos niveles: kets, matrices densidad, vector de Bloch y magnitudes derivadas.
// Convención del qubit: |0⟩ = espín arriba en z (primer elemento) y ρ = ½·(I + r·σ).

import { zeros, identity, matrix, multiply, dagger, trace, add, scale, hermitianEigenvalues } from './cmatrix.js';

export const PAULI = Object.freeze({
  I: identity(2),
  X: matrix([[0, 1], [1, 0]]),
  Y: matrix([[0, [0, -1]], [[0, 1], 0]]),
  Z: matrix([[1, 0], [0, -1]])
});

// Ket normalizado a partir de amplitudes; cada una es un número real o un par [real, imaginaria].
export function ket(amplitudes) {
  const re = Float64Array.from(amplitudes, (value) => (Array.isArray(value) ? value[0] : value));
  const im = Float64Array.from(amplitudes, (value) => (Array.isArray(value) ? value[1] : 0));
  let total = 0;
  for (let index = 0; index < re.length; index += 1) total += re[index] ** 2 + im[index] ** 2;
  const norm = Math.sqrt(total);
  if (!(norm > 0)) throw new RangeError('El ket no puede ser nulo.');
  for (let index = 0; index < re.length; index += 1) {
    re[index] /= norm;
    im[index] /= norm;
  }
  return { re, im };
}

export function ketNorm(psi) {
  let total = 0;
  for (let index = 0; index < psi.re.length; index += 1) total += psi.re[index] ** 2 + psi.im[index] ** 2;
  return total;
}

// ρ = |ψ⟩⟨ψ|.
export function densityFromKet(psi) {
  const n = psi.re.length;
  const rho = zeros(n);
  for (let row = 0; row < n; row += 1) {
    for (let column = 0; column < n; column += 1) {
      rho.re[row * n + column] = psi.re[row] * psi.re[column] + psi.im[row] * psi.im[column];
      rho.im[row * n + column] = psi.im[row] * psi.re[column] - psi.re[row] * psi.im[column];
    }
  }
  return rho;
}

// Estado puro de un qubit con ángulos polares (θ, φ): cos(θ/2)|0⟩ + e^{iφ}·sen(θ/2)|1⟩.
export function qubitKet(theta, phi = 0) {
  return ket([Math.cos(theta / 2), [Math.sin(theta / 2) * Math.cos(phi), Math.sin(theta / 2) * Math.sin(phi)]]);
}

export function blochToDensity([x, y, z]) {
  return matrix([[0.5 * (1 + z), [0.5 * x, -0.5 * y]], [[0.5 * x, 0.5 * y], 0.5 * (1 - z)]]);
}

export function blochVector(rho) {
  return [2 * rho.re[1], -2 * rho.im[1], rho.re[0] - rho.re[3]];
}

export function unitVector([x, y, z]) {
  const norm = Math.hypot(x, y, z);
  if (!(norm > 0)) throw new RangeError('El eje no puede ser nulo.');
  return [x / norm, y / norm, z / norm];
}

// Eje a partir de sus ángulos polares en grados.
export function axisFromAngles(thetaDegrees, phiDegrees = 0) {
  const theta = (thetaDegrees * Math.PI) / 180;
  const phi = (phiDegrees * Math.PI) / 180;
  return [Math.sin(theta) * Math.cos(phi), Math.sin(theta) * Math.sin(phi), Math.cos(theta)];
}

// n·σ para un eje unitario.
export function spinAlong(axis) {
  const [x, y, z] = unitVector(axis);
  return matrix([[z, [x, -y]], [[x, y], -z]]);
}

// Rotación de un qubit: U = exp(−i·ángulo·n·σ/2). Gira el vector de Bloch «ángulo» radianes alrededor de n.
export function rotation(axis, angle) {
  const [x, y, z] = unitVector(axis);
  const c = Math.cos(angle / 2);
  const s = Math.sin(angle / 2);
  return matrix([[[c, -s * z], [-s * y, -s * x]], [[s * y, -s * x], [c, s * z]]]);
}

// UρU†.
export function evolve(rho, unitary) {
  return multiply(multiply(unitary, rho), dagger(unitary));
}

// Valor esperado Tr(ρA) de un observable hermítico.
export function expectation(rho, operator) {
  return trace(multiply(rho, operator))[0];
}

export function purity(rho) {
  return trace(multiply(rho, rho))[0];
}

// Entropía de von Neumann en bits.
export function entropy(rho) {
  let total = 0;
  for (const value of hermitianEigenvalues(rho)) if (value > 1e-14) total -= value * Math.log2(value);
  return total;
}

// Fidelidad ⟨ψ|ρ|ψ⟩ con un estado puro.
export function fidelity(rho, psi) {
  const n = psi.re.length;
  let total = 0;
  for (let row = 0; row < n; row += 1) {
    for (let column = 0; column < n; column += 1) {
      // conj(ψ_row)·ρ_{row,column}·ψ_column, parte real.
      const ar = psi.re[row] * rho.re[row * n + column] + psi.im[row] * rho.im[row * n + column];
      const ai = psi.re[row] * rho.im[row * n + column] - psi.im[row] * rho.re[row * n + column];
      total += ar * psi.re[column] - ai * psi.im[column];
    }
  }
  return total;
}

// Distancia de traza ½·‖ρ − σ‖₁.
export function traceDistance(a, b) {
  let total = 0;
  for (const value of hermitianEigenvalues(add(a, scale(b, -1)))) total += Math.abs(value);
  return 0.5 * total;
}

// Mezcla estadística Σ pₖ·ρₖ.
export function mixture(terms) {
  let result = zeros(terms[0][1].n);
  for (const [weight, rho] of terms) result = add(result, scale(rho, weight));
  return result;
}

// Traza parcial de un sistema compuesto de dimensiones «dims», quedándose con el subsistema «keep».
export function partialTrace(rho, dims, keep) {
  const total = dims.reduce((product, value) => product * value, 1);
  if (total !== rho.n) throw new RangeError('Las dimensiones no coinciden con la matriz.');
  const size = dims[keep];
  const result = zeros(size);
  const strides = dims.map((_, index) => dims.slice(index + 1).reduce((product, value) => product * value, 1));
  const others = total / size;
  for (let rest = 0; rest < others; rest += 1) {
    // Reparte «rest» entre los índices de los subsistemas que se trazan.
    let remainder = rest;
    let base = 0;
    for (let index = dims.length - 1; index >= 0; index -= 1) {
      if (index === keep) continue;
      base += (remainder % dims[index]) * strides[index];
      remainder = Math.floor(remainder / dims[index]);
    }
    for (let row = 0; row < size; row += 1) {
      for (let column = 0; column < size; column += 1) {
        const source = (base + row * strides[keep]) * total + (base + column * strides[keep]);
        result.re[row * size + column] += rho.re[source];
        result.im[row * size + column] += rho.im[source];
      }
    }
  }
  return result;
}
