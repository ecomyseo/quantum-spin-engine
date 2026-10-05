// Rabi oscillations of a driven qubit with detuning, in the rotating frame:
//   H = (Δ/2)·σz + (Ω/2)·σx   →   P₁(t) = (Ω/Ω′)²·sin²(Ω′t/2),  Ω′ = √(Ω² + Δ²)
// Run: node examples/rabi.mjs
import { blochToDensity, blochVector, rotation, evolve } from '../src/index.js';

const omega = 1.5;
const delta = 0.8;
const dt = 0.02;
const generalized = Math.hypot(omega, delta);
const step = rotation([omega, 0, delta], generalized * dt);

let rho = blochToDensity([0, 0, 1]); // |0⟩
let worst = 0;
for (let index = 1; index <= 600; index += 1) {
  rho = evolve(rho, step);
  const time = index * dt;
  const population = 0.5 * (1 - blochVector(rho)[2]);
  const law = (omega / generalized) ** 2 * Math.sin((generalized * time) / 2) ** 2;
  worst = Math.max(worst, Math.abs(population - law));
  if (index % 100 === 0) console.log(`t = ${time.toFixed(1)}   P(|1>) = ${population.toFixed(6)}   law = ${law.toFixed(6)}`);
}
console.log(`maximum transfer Ω²/Ω′² = ${((omega / generalized) ** 2).toFixed(6)}`);
console.log(`largest deviation from the Rabi law over 600 steps: ${worst.toExponential(2)}`);
if (worst > 1e-9) process.exit(1);
