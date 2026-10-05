// Espinor de espín ½ sobre una rejilla 1D: dos funciones de onda, una por componente del espín a lo largo
// de z, cada una con su propio potencial. Es lo que hace falta para un imán de Stern–Gerlach, donde la
// fuerza tiene sentido contrario para «arriba» y «abajo» y el espín acaba entrelazado con la posición.

import { createState } from './grid.js';
import { createPropagator } from './split-step.js';
import { zeros } from './cmatrix.js';

// Espinor producto: la misma forma espacial multiplicada por las amplitudes de espín (a, b),
// cada una un par [real, imaginaria]. Se supone |a|² + |b|² = 1 y el paquete normalizado.
export function createSpinor(grid, packet, [a, b]) {
  const up = createState(grid);
  const down = createState(grid);
  for (let index = 0; index < grid.size; index += 1) {
    up.re[index] = a[0] * packet.re[index] - a[1] * packet.im[index];
    up.im[index] = a[0] * packet.im[index] + a[1] * packet.re[index];
    down.re[index] = b[0] * packet.re[index] - b[1] * packet.im[index];
    down.im[index] = b[0] * packet.im[index] + b[1] * packet.re[index];
  }
  return { up, down };
}

// Propagador: cada componente avanza con su potencial y la misma energía cinética.
export function createSpinorPropagator(grid, { mass = 1 } = {}) {
  const propagator = createPropagator(grid, { mass });
  return {
    step(spinor, potentialUp, potentialDown, dt) {
      propagator.step(spinor.up, potentialUp, dt);
      propagator.step(spinor.down, potentialDown, dt);
      return spinor;
    }
  };
}

// Densidades de probabilidad: total y de cada componente.
export function spinorDensity(grid, spinor) {
  const up = new Float64Array(grid.size);
  const down = new Float64Array(grid.size);
  const total = new Float64Array(grid.size);
  for (let index = 0; index < grid.size; index += 1) {
    up[index] = spinor.up.re[index] ** 2 + spinor.up.im[index] ** 2;
    down[index] = spinor.down.re[index] ** 2 + spinor.down.im[index] ** 2;
    total[index] = up[index] + down[index];
  }
  return { up, down, total };
}

// Peso, posición media, anchura y momento medio de una componente (ħ = 1).
export function componentMoments(grid, state) {
  const { size, x, dx } = grid;
  let weight = 0;
  let first = 0;
  let second = 0;
  for (let index = 0; index < size; index += 1) {
    const density = state.re[index] ** 2 + state.im[index] ** 2;
    weight += density;
    first += density * x[index];
    second += density * x[index] * x[index];
  }
  weight *= dx;
  if (!(weight > 1e-300)) return { weight: 0, mean: 0, sigma: 0, momentum: 0 };
  const mean = (first * dx) / weight;
  // Momento medio en el espacio de Fourier, exacto sobre la rejilla: Σ k·|φ(k)|² / Σ |φ(k)|².
  const re = Float64Array.from(state.re);
  const im = Float64Array.from(state.im);
  grid.fft.forward(re, im);
  let spectrum = 0;
  let momentum = 0;
  for (let index = 0; index < size; index += 1) {
    const power = re[index] ** 2 + im[index] ** 2;
    spectrum += power;
    momentum += power * grid.k[index];
  }
  return { weight, mean, sigma: Math.sqrt(Math.max(0, (second * dx) / weight - mean * mean)), momentum: momentum / spectrum };
}

// Estado del espín solo, trazando la posición: ρ₀₀ = ∫|ψ↑|², ρ₁₁ = ∫|ψ↓|², ρ₀₁ = ∫ψ↑·ψ↓*.
// Cuando los dos haces dejan de solaparse, ρ₀₁ se anula: el espín, por sí solo, queda en una mezcla.
export function reducedSpin(grid, spinor) {
  const rho = zeros(2);
  let coherenceRe = 0;
  let coherenceIm = 0;
  let up = 0;
  let down = 0;
  for (let index = 0; index < grid.size; index += 1) {
    up += spinor.up.re[index] ** 2 + spinor.up.im[index] ** 2;
    down += spinor.down.re[index] ** 2 + spinor.down.im[index] ** 2;
    coherenceRe += spinor.up.re[index] * spinor.down.re[index] + spinor.up.im[index] * spinor.down.im[index];
    coherenceIm += spinor.up.im[index] * spinor.down.re[index] - spinor.up.re[index] * spinor.down.im[index];
  }
  rho.re[0] = up * grid.dx;
  rho.re[3] = down * grid.dx;
  rho.re[1] = coherenceRe * grid.dx;
  rho.im[1] = coherenceIm * grid.dx;
  rho.re[2] = rho.re[1];
  rho.im[2] = -rho.im[1];
  return rho;
}
