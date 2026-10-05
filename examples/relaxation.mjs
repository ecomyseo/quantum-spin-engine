// Energy relaxation (T₁) and dephasing (T₂) of a qubit with the Lindblad master equation,
// and the same decay recovered as the average of quantum-jump trajectories.
// Run: node examples/relaxation.mjs
import { PAULI, ket, densityFromKet, blochVector, matrix, zeros, add, scale, trace, liouvillian, superPropagator, applySuper, createJumpUnraveling, createRng } from '../src/index.js';

const T1 = 6;
const Tphi = 8;
const T2 = 1 / (1 / (2 * T1) + 1 / Tphi);
const dt = 0.05;
const steps = 120;
const lower = matrix([[0, 1], [0, 0]]); // |0⟩⟨1|
const hamiltonian = zeros(2);
const collapse = [{ operator: lower, rate: 1 / T1 }, { operator: PAULI.Z, rate: 1 / (2 * Tphi) }];

// Exact propagation of the ensemble: ρ(t + dt) = exp(𝓛·dt)·ρ(t).
const propagator = superPropagator(liouvillian(hamiltonian, collapse), dt);
const start = ket([1, 1]); // |+⟩
let rho = densityFromKet(start);
for (let index = 0; index < steps; index += 1) rho = applySuper(propagator, rho);
const time = steps * dt;
const [x, y, z] = blochVector(rho);
console.log(`T2 = ${T2.toFixed(3)}  (always ≤ 2·T1 = ${2 * T1})`);
console.log(`t = ${time}: coherence = ${Math.hypot(x, y).toFixed(6)}   law e^(-t/T2) = ${Math.exp(-time / T2).toFixed(6)}`);
console.log(`t = ${time}: P(|1>)    = ${(0.5 * (1 - z)).toFixed(6)}   law ½·e^(-t/T1) = ${(0.5 * Math.exp(-time / T1)).toFixed(6)}`);
console.log(`trace = ${trace(rho)[0].toFixed(12)}`);

// Unravelling: each trajectory is a pure state with sudden jumps; their average is ρ(t).
const unraveling = createJumpUnraveling(hamiltonian, collapse, dt);
const rng = createRng('relaxation-example');
const runs = 4000;
let average = zeros(2);
for (let run = 0; run < runs; run += 1) {
  const trajectory = unraveling.start(start, rng);
  for (let index = 0; index < steps; index += 1) unraveling.step(trajectory, rng);
  average = add(average, scale(densityFromKet(unraveling.state(trajectory)), 1 / runs));
}
const sampled = 0.5 * (1 - blochVector(average)[2]);
console.log(`average of ${runs} jump trajectories: P(|1>) = ${sampled.toFixed(4)}  (statistical error ≈ ${(1 / Math.sqrt(runs)).toFixed(3)})`);
if (Math.abs(Math.hypot(x, y) - Math.exp(-time / T2)) > 1e-9 || Math.abs(sampled - 0.5 * (1 - z)) > 0.04) process.exit(1);
