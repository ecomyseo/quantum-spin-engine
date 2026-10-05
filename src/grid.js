import { createFft } from './fft.js';

// Rejilla periódica uniforme en x ∈ [−L/2, L/2) y su rejilla conjugada de números de onda.
export function createGrid({ size = 512, length = 40 } = {}) {
  const fft = createFft(size);
  const dx = length / size;
  const dk = (2 * Math.PI) / length;
  const x = new Float64Array(size);
  const k = new Float64Array(size);
  for (let index = 0; index < size; index += 1) {
    x[index] = -length / 2 + index * dx;
    k[index] = dk * (index < size / 2 ? index : index - size);
  }
  return Object.freeze({ size, length, dx, dk, x, k, fft });
}

export function createState(grid) {
  return { re: new Float64Array(grid.size), im: new Float64Array(grid.size) };
}

export function cloneState(state) {
  return { re: Float64Array.from(state.re), im: Float64Array.from(state.im) };
}

export function copyInto(target, source) {
  target.re.set(source.re);
  target.im.set(source.im);
  return target;
}

export function stateNorm(grid, state) {
  let total = 0;
  for (let index = 0; index < grid.size; index += 1) {
    total += state.re[index] * state.re[index] + state.im[index] * state.im[index];
  }
  return total * grid.dx;
}

export function scaleState(state, factor) {
  for (let index = 0; index < state.re.length; index += 1) {
    state.re[index] *= factor;
    state.im[index] *= factor;
  }
  return state;
}

// Normaliza en sitio y devuelve la norma que tenía el estado antes de normalizar.
export function normalizeState(grid, state) {
  const previous = stateNorm(grid, state);
  if (previous > 0) scaleState(state, 1 / Math.sqrt(previous));
  return previous;
}

// Paquete gaussiano ψ ∝ exp(−(x−x0)²/(4σ²))·exp(i·k0·(x−x0) + i·φ); σ es la desviación típica de |ψ|².
export function gaussianPacket(grid, { center = 0, sigma = 1, waveNumber = 0, phase = 0 } = {}) {
  const state = createState(grid);
  for (let index = 0; index < grid.size; index += 1) {
    const offset = grid.x[index] - center;
    const envelope = Math.exp(-(offset * offset) / (4 * sigma * sigma));
    const angle = waveNumber * offset + phase;
    state.re[index] = envelope * Math.cos(angle);
    state.im[index] = envelope * Math.sin(angle);
  }
  normalizeState(grid, state);
  return state;
}

// Combinación lineal compleja Σ c_j·ψ_j, con c_j = re + i·im.
export function combineStates(grid, terms, { normalize = true } = {}) {
  const state = createState(grid);
  for (const term of terms) {
    const cr = term.re ?? 1;
    const ci = term.im ?? 0;
    for (let index = 0; index < grid.size; index += 1) {
      state.re[index] += cr * term.state.re[index] - ci * term.state.im[index];
      state.im[index] += cr * term.state.im[index] + ci * term.state.re[index];
    }
  }
  if (normalize) normalizeState(grid, state);
  return state;
}

export function stateFromReal(values) {
  return { re: Float64Array.from(values), im: new Float64Array(values.length) };
}
