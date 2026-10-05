// A Stern–Gerlach magnet in one dimension: the gradient pushes spin-up and spin-down in opposite
// directions, the packet splits in two, and the spin on its own decoheres because it is now
// entangled with the position.
//   H = p²/2m − f·z·σz        beam separation = f·t²/m
// Run: node examples/stern-gerlach.mjs
import { createGrid, gaussianPacket, createSpinor, createSpinorPropagator, componentMoments, reducedSpin, blochVector, purity } from '../src/index.js';

const grid = createGrid({ size: 1024, length: 200 });
const sigma = 1.5;
const force = 0.2;
const mass = 1;
const tilt = (60 * Math.PI) / 180; // spin tilted 60° from z: weights cos²30° and sin²30°
const spinor = createSpinor(grid, gaussianPacket(grid, { center: 0, sigma }), [[Math.cos(tilt / 2), 0], [Math.sin(tilt / 2), 0]]);
const up = Float64Array.from(grid.x, (z) => -force * z);
const down = Float64Array.from(grid.x, (z) => force * z);
const propagator = createSpinorPropagator(grid, { mass });

const dt = 0.025;
let worst = 0;
for (let index = 1; index <= 480; index += 1) {
  propagator.step(spinor, up, down, dt);
  if (index % 120 !== 0) continue;
  const time = index * dt;
  const upper = componentMoments(grid, spinor.up);
  const lower = componentMoments(grid, spinor.down);
  const spin = reducedSpin(grid, spinor);
  const separation = upper.mean - lower.mean;
  worst = Math.max(worst, Math.abs(separation - (force * time * time) / mass));
  console.log(`t = ${time.toFixed(0).padStart(2)}  weights ${upper.weight.toFixed(4)} / ${lower.weight.toFixed(4)}  separation ${separation.toFixed(3)} (law ${((force * time * time) / mass).toFixed(3)})  |<σ>| = ${Math.hypot(...blochVector(spin)).toFixed(4)}  purity ${purity(spin).toFixed(4)}`);
}
console.log(`largest deviation from f·t²/m: ${worst.toExponential(2)}`);
if (worst > 5e-3) process.exit(1);
