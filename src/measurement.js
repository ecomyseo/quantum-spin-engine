// Medición de un qubit: proyectiva a lo largo de cualquier eje, estadística de repeticiones y tomografía.

import { PAULI, spinAlong, unitVector, expectation } from './density.js';
import { add, scale, multiply } from './cmatrix.js';

// Proyectores sobre «arriba» y «abajo» a lo largo de un eje: P± = ½·(I ± n·σ).
export function axisProjectors(axis) {
  const along = spinAlong(axis);
  return { up: scale(add(PAULI.I, along), 0.5), down: scale(add(PAULI.I, scale(along, -1)), 0.5) };
}

// Medición proyectiva a lo largo de un eje: probabilidades y estado posterior a cada resultado.
export function measureAlong(rho, axis) {
  const [x, y, z] = unitVector(axis);
  const { up, down } = axisProjectors(axis);
  const along = expectation(rho, spinAlong(axis));
  const pUp = Math.min(1, Math.max(0, 0.5 * (1 + along)));
  return {
    pUp,
    pDown: 1 - pUp,
    // Estados tras el resultado: el vector de Bloch queda sobre el eje, en un sentido u otro.
    upState: [x, y, z],
    downState: [-x, -y, -z],
    // Estado si la medición se hace y no se lee: mezcla de los dos resultados.
    unread: add(multiply(multiply(up, rho), up), multiply(multiply(down, rho), down))
  };
}

// Número de resultados «arriba» en «shots» repeticiones independientes.
export function sampleUpCount(probability, shots, rng) {
  let count = 0;
  for (let index = 0; index < shots; index += 1) if (rng.uniform() < probability) count += 1;
  return count;
}

// Tomografía de un qubit: se mide σx, σy y σz sobre «shots» copias cada uno y se estima el vector de Bloch.
export function tomography(rho, shots, rng) {
  const axes = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };
  const estimate = {};
  const counts = {};
  for (const [name, axis] of Object.entries(axes)) {
    const { pUp } = measureAlong(rho, axis);
    counts[name] = sampleUpCount(pUp, shots, rng);
    estimate[name] = (2 * counts[name]) / shots - 1;
  }
  const length = Math.hypot(estimate.x, estimate.y, estimate.z);
  return {
    shots,
    counts,
    vector: [estimate.x, estimate.y, estimate.z],
    length,
    // Con pocas repeticiones la estimación puede salirse de la esfera: no sería un estado físico.
    physical: length <= 1
  };
}
