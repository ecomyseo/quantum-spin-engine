// Propagador de la ecuación de Schrödinger 1D por operador partido (Strang) con Fourier.
//   ψ(t+dt) = e^{−iV·dt/2ħ} · F⁻¹ e^{−iħk²·dt/2m} F · e^{−iV·dt/2ħ} ψ(t)
// Cada factor es unitario, de modo que la norma se conserva hasta el redondeo.
// Un paso con dt negativo es el inverso exacto del paso con dt positivo.

export function createPropagator(grid, { mass = 1, hbar = 1 } = {}) {
  const { size, k, x, fft } = grid;
  const kineticCache = new Map();
  const potentialCache = new WeakMap();

  function kineticFactors(dt) {
    let entry = kineticCache.get(dt);
    if (!entry) {
      const cos = new Float64Array(size);
      const sin = new Float64Array(size);
      for (let index = 0; index < size; index += 1) {
        const angle = (-hbar * k[index] * k[index] * dt) / (2 * mass);
        cos[index] = Math.cos(angle);
        sin[index] = Math.sin(angle);
      }
      entry = { cos, sin };
      if (kineticCache.size > 12) kineticCache.clear();
      kineticCache.set(dt, entry);
    }
    return entry;
  }

  function potentialFactors(potential, dt, cache) {
    if (cache) {
      const cached = potentialCache.get(potential);
      if (cached && cached.dt === dt) return cached;
    }
    const cos = new Float64Array(size);
    const sin = new Float64Array(size);
    for (let index = 0; index < size; index += 1) {
      const angle = (-potential[index] * dt) / (2 * hbar);
      cos[index] = Math.cos(angle);
      sin[index] = Math.sin(angle);
    }
    const entry = { dt, cos, sin };
    if (cache) potentialCache.set(potential, entry);
    return entry;
  }

  function rotate(state, cos, sin) {
    const { re, im } = state;
    for (let index = 0; index < size; index += 1) {
      const real = re[index];
      re[index] = real * cos[index] - im[index] * sin[index];
      im[index] = real * sin[index] + im[index] * cos[index];
    }
  }

  // Paso en tiempo real. «cache: false» es obligatorio si el potencial cambia en sitio entre pasos.
  function step(state, potential, dt, { cache = true } = {}) {
    const halfPotential = potential ? potentialFactors(potential, dt, cache) : null;
    if (halfPotential) rotate(state, halfPotential.cos, halfPotential.sin);
    const kinetic = kineticFactors(dt);
    fft.forward(state.re, state.im);
    rotate(state, kinetic.cos, kinetic.sin);
    fft.inverse(state.re, state.im);
    if (halfPotential) rotate(state, halfPotential.cos, halfPotential.sin);
    return state;
  }

  // Paso en tiempo imaginario: aplica e^{−H·dτ/ħ}, que no es unitario, y renormaliza.
  // Devuelve la norma previa a la renormalización, útil para estimar la energía.
  function stepImaginary(state, potential, dtau) {
    const { re, im } = state;
    const damp = (factor) => {
      for (let index = 0; index < size; index += 1) {
        const scale = Math.exp(factor(index));
        re[index] *= scale;
        im[index] *= scale;
      }
    };
    if (potential) damp((index) => (-potential[index] * dtau) / (2 * hbar));
    fft.forward(re, im);
    damp((index) => (-hbar * k[index] * k[index] * dtau) / (2 * mass));
    fft.inverse(re, im);
    if (potential) damp((index) => (-potential[index] * dtau) / (2 * hbar));
    let total = 0;
    for (let index = 0; index < size; index += 1) total += re[index] * re[index] + im[index] * im[index];
    total *= grid.dx;
    const scale = total > 0 ? 1 / Math.sqrt(total) : 1;
    for (let index = 0; index < size; index += 1) {
      re[index] *= scale;
      im[index] *= scale;
    }
    return total;
  }

  // Impulso uniforme e^{iκx}: desplaza todo el espectro de momentos en ħκ.
  function kick(state, kappa) {
    const { re, im } = state;
    for (let index = 0; index < size; index += 1) {
      const angle = kappa * x[index];
      const cos = Math.cos(angle);
      const sin = Math.sin(angle);
      const real = re[index];
      re[index] = real * cos - im[index] * sin;
      im[index] = real * sin + im[index] * cos;
    }
    return state;
  }

  // Fase local e^{i·escala·fases(x)}; se usa para potenciales estocásticos de un solo paso.
  function applyPhase(state, phases, scale = 1) {
    const { re, im } = state;
    for (let index = 0; index < size; index += 1) {
      const angle = scale * phases[index];
      const cos = Math.cos(angle);
      const sin = Math.sin(angle);
      const real = re[index];
      re[index] = real * cos - im[index] * sin;
      im[index] = real * sin + im[index] * cos;
    }
    return state;
  }

  return Object.freeze({ mass, hbar, step, stepImaginary, kick, applyPhase });
}
