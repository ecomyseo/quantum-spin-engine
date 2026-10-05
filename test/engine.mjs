// Motor de sistemas cuánticos de pocos niveles: álgebra de matrices complejas, estados,
// ecuación maestra de Lindblad, canales de Kraus, trayectorias de saltos, medición y espinor de Stern–Gerlach.
// Cada resultado numérico se contrasta con una ley exacta o con un segundo método independiente.
import test from 'node:test';
import assert from 'node:assert/strict';
import { zeros, identity, matrix, add, scale, multiply, dagger, trace, kron, expm, distance, hermitianEigenvalues, commutator } from '../src/cmatrix.js';
import { PAULI, ket, qubitKet, densityFromKet, blochVector, blochToDensity, rotation, evolve, expectation, purity, entropy, fidelity, traceDistance, partialTrace, spinAlong, axisFromAngles, mixture } from '../src/density.js';
import { liouvillian, superPropagator, applySuper, applyKraus, krausDefect, channel, createJumpUnraveling } from '../src/open-system.js';
import { measureAlong, axisProjectors, sampleUpCount, tomography } from '../src/measurement.js';
import { createSpinor, createSpinorPropagator, spinorDensity, componentMoments, reducedSpin } from '../src/spinor.js';
import { createGrid, gaussianPacket } from '../src/grid.js';
import { createRng } from '../src/rng.js';

const near = (actual, expected, tolerance, message) => assert.ok(Math.abs(actual - expected) <= tolerance, `${message || ''} ${actual} frente a ${expected}`);

test('la exponencial de matrices reproduce las rotaciones de espín y es unitaria', () => {
  for (const [axis, angle] of [[[1, 0, 0], 0.7], [[0, 1, 0], 2.9], [[0, 0, 1], -5.1], [[1, 2, -3], 11.3]]) {
    const generator = scale(spinAlong(axis), 0, -angle / 2);
    const exponential = expm(generator);
    assert.ok(distance(exponential, rotation(axis, angle)) < 1e-13);
    assert.ok(distance(multiply(dagger(exponential), exponential), identity(2)) < 1e-13);
  }
  // Matriz no normal: exp([[0,1],[0,0]]) = [[1,1],[0,1]].
  assert.ok(distance(expm(matrix([[0, 1], [0, 0]])), matrix([[1, 1], [0, 1]])) < 1e-15);
  // Norma grande: exige reducir y elevar al cuadrado.
  const big = expm(scale(PAULI.Y, 0, -40));
  assert.ok(distance(big, rotation([0, 1, 0], 80)) < 1e-10);
});

test('los autovalores hermíticos coinciden con los casos exactos y con la traza', () => {
  assert.deepEqual(hermitianEigenvalues(PAULI.Y).map((value) => Math.round(value * 1e12) / 1e12), [-1, 1]);
  const known = hermitianEigenvalues(matrix([[2, [1, -1]], [[1, 1], 3]]));
  near(known[0], 1, 1e-12);
  near(known[1], 4, 1e-12);
  const rng = createRng('autovalores');
  for (const size of [3, 4, 8]) {
    const random = zeros(size);
    for (let index = 0; index < size * size; index += 1) {
      random.re[index] = rng.gaussian();
      random.im[index] = rng.gaussian();
    }
    const hermitian = add(random, dagger(random));
    const values = hermitianEigenvalues(hermitian);
    near(values.reduce((sum, value) => sum + value, 0), trace(hermitian)[0], 1e-10, 'traza');
    near(values.reduce((sum, value) => sum + value * value, 0), trace(multiply(hermitian, hermitian))[0], 1e-9, 'traza del cuadrado');
    assert.ok(values.every((value, index) => index === 0 || value >= values[index - 1]));
  }
});

test('el vector de Bloch gira como dicta la fórmula de Rodrigues', () => {
  const rodrigues = ([vx, vy, vz], [kx, ky, kz], angle) => {
    const norm = Math.hypot(kx, ky, kz);
    const k = [kx / norm, ky / norm, kz / norm];
    const dot = k[0] * vx + k[1] * vy + k[2] * vz;
    const cross = [k[1] * vz - k[2] * vy, k[2] * vx - k[0] * vz, k[0] * vy - k[1] * vx];
    return [0, 1, 2].map((index) => [vx, vy, vz][index] * Math.cos(angle) + cross[index] * Math.sin(angle) + k[index] * dot * (1 - Math.cos(angle)));
  };
  const start = axisFromAngles(37, 112);
  let rho = blochToDensity(start);
  for (const [axis, angle] of [[[0, 0, 1], 1.1], [[1, 1, 0], 2.3], [[0.2, -0.7, 0.4], -0.9]]) {
    const before = blochVector(rho);
    rho = evolve(rho, rotation(axis, angle));
    const expected = rodrigues(before, axis, angle);
    blochVector(rho).forEach((value, index) => near(value, expected[index], 1e-13));
  }
  near(purity(rho), 1, 1e-13);
  const psi = qubitKet(0.8, 2.1);
  const fromKet = blochVector(densityFromKet(psi));
  near(fromKet[0], Math.sin(0.8) * Math.cos(2.1), 1e-14);
  near(fromKet[1], Math.sin(0.8) * Math.sin(2.1), 1e-14);
  near(fromKet[2], Math.cos(0.8), 1e-14);
  near(fidelity(densityFromKet(psi), psi), 1, 1e-14);
  assert.ok(distance(commutator(PAULI.X, PAULI.Y), scale(PAULI.Z, 0, 2)) < 1e-15);
});

test('la ecuación maestra reproduce T1, T2 y el estado estacionario', () => {
  const lower = matrix([[0, 1], [0, 0]]);
  const start = densityFromKet(qubitKet(1.9, 0.6));
  const gamma = 0.37;
  const phase = 0.21;
  const superoperator = liouvillian(zeros(2), [{ operator: lower, rate: gamma }, { operator: PAULI.Z, rate: phase }]);
  for (const time of [0.3, 1.7, 6]) {
    const rho = applySuper(superPropagator(superoperator, time), start);
    near(rho.re[3], start.re[3] * Math.exp(-gamma * time), 1e-12, 'población excitada');
    near(Math.hypot(rho.re[1], rho.im[1]), Math.hypot(start.re[1], start.im[1]) * Math.exp(-(gamma / 2 + 2 * phase) * time), 1e-12, 'coherencia');
    near(trace(rho)[0], 1, 1e-12, 'traza');
    assert.ok(hermitianEigenvalues(rho)[0] > -1e-12);
  }
  // Dos pasos de Δt equivalen a uno de 2Δt, y un campo sobre z hace girar la fase de la coherencia.
  const field = liouvillian(scale(PAULI.Z, 0.5 * 1.3), [{ operator: PAULI.Z, rate: phase }]);
  const twice = applySuper(superPropagator(field, 0.4), applySuper(superPropagator(field, 0.4), start));
  const once = applySuper(superPropagator(field, 0.8), start);
  assert.ok(distance(twice, once) < 1e-13);
  near(Math.atan2(once.im[1], once.re[1]), Math.atan2(start.im[1], start.re[1]) - 1.3 * 0.8, 1e-12, 'fase de la coherencia');
  // Qubit conducido con relajación: población estacionaria Ω²/(γ² + 2Ω²).
  const omega = 0.9;
  const driven = liouvillian(scale(PAULI.X, 0.5 * omega), [{ operator: lower, rate: gamma }]);
  const steady = applySuper(superPropagator(driven, 400), start);
  near(steady.re[3], (omega * omega) / (gamma * gamma + 2 * omega * omega), 1e-10, 'estado estacionario');
});

test('los canales de Kraus conservan la traza y coinciden con la ecuación maestra', () => {
  const lower = matrix([[0, 1], [0, 0]]);
  const start = densityFromKet(qubitKet(2.2, -1));
  for (const strength of [0, 0.3, 1]) {
    for (const kraus of [channel.dephasing(strength), channel.amplitudeDamping(strength), channel.depolarizing(strength)]) assert.ok(krausDefect(kraus) < 1e-15);
  }
  const time = 0.9;
  const gamma = 0.6;
  const damped = applySuper(superPropagator(liouvillian(zeros(2), [{ operator: lower, rate: gamma }]), time), start);
  assert.ok(distance(damped, applyKraus(start, channel.amplitudeDamping(1 - Math.exp(-gamma * time)))) < 1e-13);
  const dephased = applySuper(superPropagator(liouvillian(zeros(2), [{ operator: PAULI.Z, rate: gamma }]), time), start);
  assert.ok(distance(dephased, applyKraus(start, channel.dephasing(1 - Math.exp(-4 * gamma * time)))) < 1e-13);
  const depolarized = blochVector(applyKraus(start, channel.depolarizing(0.4)));
  blochVector(start).forEach((value, index) => near(depolarized[index], 0.6 * value, 1e-14));
});

test('el promedio de trayectorias de saltos reproduce la ecuación maestra', () => {
  const lower = matrix([[0, 1], [0, 0]]);
  const hamiltonian = scale(PAULI.X, 0.5 * 1.4);
  const collapse = [{ operator: lower, rate: 0.5 }, { operator: PAULI.Z, rate: 0.15 }];
  const dt = 0.02;
  const steps = 120;
  const unraveling = createJumpUnraveling(hamiltonian, collapse, dt);
  const rng = createRng('trayectorias');
  const psi = qubitKet(0, 0);
  const total = 3000;
  let average = zeros(2);
  let jumps = 0;
  for (let run = 0; run < total; run += 1) {
    const trajectory = unraveling.start(psi, rng);
    for (let step = 0; step < steps; step += 1) unraveling.step(trajectory, rng);
    jumps += trajectory.jumps;
    const state = unraveling.state(trajectory);
    near(state.re[0] ** 2 + state.im[0] ** 2 + state.re[1] ** 2 + state.im[1] ** 2, 1, 1e-12);
    average = add(average, scale(densityFromKet(state), 1 / total));
  }
  const exact = applySuper(superPropagator(liouvillian(hamiltonian, collapse), dt * steps), densityFromKet(psi));
  // Error estadístico de 3000 trayectorias: del orden de 1/√3000 ≈ 0,018.
  assert.ok(distance(average, exact) < 0.05, `distancia ${distance(average, exact)}`);
  assert.ok(jumps > 0);
});

test('la medición a lo largo de un eje sigue cos²(γ/2) y deja el estado sobre el eje', () => {
  const state = axisFromAngles(50, 20);
  const rho = blochToDensity(state);
  for (const axis of [[0, 0, 1], [1, 0, 0], axisFromAngles(130, 250)]) {
    const result = measureAlong(rho, axis);
    const cosine = state[0] * axis[0] + state[1] * axis[1] + state[2] * axis[2];
    near(result.pUp, 0.5 * (1 + cosine), 1e-14);
    near(result.pUp + result.pDown, 1, 1e-15);
    near(purity(result.unread), 0.5 * (1 + cosine * cosine), 1e-14);
    const { up, down } = axisProjectors(axis);
    assert.ok(distance(multiply(up, up), up) < 1e-15);
    assert.ok(distance(multiply(up, down), zeros(2)) < 1e-15);
    near(expectation(blochToDensity(result.upState), spinAlong(axis)), 1, 1e-14);
  }
  const rng = createRng('disparos');
  const count = sampleUpCount(0.3, 20000, rng);
  assert.ok(Math.abs(count / 20000 - 0.3) < 4 * Math.sqrt((0.3 * 0.7) / 20000));
});

test('la tomografía recupera el vector de Bloch dentro de su error estadístico', () => {
  const vector = [0.5, -0.3, 0.6];
  const rho = blochToDensity(vector);
  const shots = 4000;
  const result = tomography(rho, shots, createRng('tomografia'));
  result.vector.forEach((value, index) => assert.ok(Math.abs(value - vector[index]) < 5 * Math.sqrt((1 - vector[index] ** 2) / shots)));
  assert.equal(JSON.stringify(result), JSON.stringify(tomography(rho, shots, createRng('tomografia'))));
  assert.notEqual(JSON.stringify(result), JSON.stringify(tomography(rho, shots, createRng('otra'))));
});

test('entropía, traza parcial y distancia de traza en casos exactos', () => {
  const bell = ket([1, 0, 0, 1]);
  const pair = densityFromKet(bell);
  for (const keep of [0, 1]) {
    const reduced = partialTrace(pair, [2, 2], keep);
    assert.ok(distance(reduced, scale(identity(2), 0.5)) < 1e-15);
    near(entropy(reduced), 1, 1e-12);
  }
  near(entropy(pair), 0, 1e-9);
  const product = kron(blochToDensity([0, 0, 1]), blochToDensity([1, 0, 0]));
  assert.ok(distance(partialTrace(product, [2, 2], 0), blochToDensity([0, 0, 1])) < 1e-15);
  assert.ok(distance(partialTrace(product, [2, 2], 1), blochToDensity([1, 0, 0])) < 1e-15);
  near(traceDistance(blochToDensity([0, 0, 1]), blochToDensity([0, 0, -1])), 1, 1e-12);
  near(traceDistance(blochToDensity([0.2, 0, 0]), blochToDensity([-0.4, 0, 0])), 0.3, 1e-12);
  near(purity(mixture([[0.5, blochToDensity([0, 0, 1])], [0.5, blochToDensity([0, 0, -1])]])), 0.5, 1e-15);
});

test('un imán de Stern–Gerlach separa los haces según f·t²/m y entrelaza espín y posición', () => {
  const grid = createGrid({ size: 1024, length: 200 });
  const sigma = 1.5;
  const force = 0.2;
  const mass = 1;
  const theta = 1.1;
  const spinor = createSpinor(grid, gaussianPacket(grid, { center: 0, sigma }), [[Math.cos(theta / 2), 0], [Math.sin(theta / 2), 0]]);
  const up = Float64Array.from(grid.x, (position) => -force * position);
  const down = Float64Array.from(grid.x, (position) => force * position);
  const propagator = createSpinorPropagator(grid, { mass });
  const before = reducedSpin(grid, spinor);
  near(purity(before), 1, 1e-10, 'espín puro al principio');
  const dt = 0.02;
  const steps = 500;
  for (let step = 0; step < steps; step += 1) propagator.step(spinor, up, down, dt);
  const time = dt * steps;
  const upper = componentMoments(grid, spinor.up);
  const lower = componentMoments(grid, spinor.down);
  near(upper.weight, Math.cos(theta / 2) ** 2, 1e-9, 'peso arriba');
  near(lower.weight, Math.sin(theta / 2) ** 2, 1e-9, 'peso abajo');
  near(upper.mean, (force * time * time) / (2 * mass), 2e-3, 'centro del haz de arriba');
  near(lower.mean, -(force * time * time) / (2 * mass), 2e-3, 'centro del haz de abajo');
  near(upper.momentum, force * time, 2e-3, 'momento del haz de arriba');
  const density = spinorDensity(grid, spinor);
  near(density.total.reduce((sum, value) => sum + value, 0) * grid.dx, 1, 1e-9, 'norma');
  // Solape de los dos haces: exp(−[(2ft)²σ²/2 + (ft²/m)²/(8σ²)]) para un paquete gaussiano mínimo.
  const reduced = reducedSpin(grid, spinor);
  const coherence = Math.hypot(reduced.re[1], reduced.im[1]) / Math.sqrt(upper.weight * lower.weight);
  const expected = Math.exp(-(((2 * force * time) ** 2 * sigma * sigma) / 2 + ((force * time * time) / mass) ** 2 / (8 * sigma * sigma)));
  near(coherence, expected, 1e-6, 'coherencia del espín reducido');
  assert.ok(purity(reduced) < 0.9);
});
