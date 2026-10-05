// Live demo of quantum-spin-engine. Two experiments, both computed frame by frame with the engine:
//   · a driven qubit with detuning, relaxation (T1) and pure dephasing (Tφ);
//   · a Stern–Gerlach magnet acting on a spin-½ wave packet.
// The panel on the right compares each result with its exact law while the simulation runs.

import { PAULI, blochToDensity, blochVector, purity } from '../src/density.js';
import { matrix, add, scale, trace } from '../src/cmatrix.js';
import { liouvillian, superPropagator, applySuper } from '../src/open-system.js';
import { createGrid, gaussianPacket } from '../src/grid.js';
import { createSpinor, createSpinorPropagator, spinorDensity, componentMoments, reducedSpin } from '../src/spinor.js';
import { createRng } from '../src/rng.js';
import { createBlochRenderer } from '../src/bloch-renderer.js';

const LABELS = { caption: ['Bloch sphere: a map of states, not a spinning object'], captionNarrow: ['A map of states,', 'not a spinning object'], time: 'time →', height: 'z ↑', absorbed: '% absorbed', decimal: '.' };
const FOREVER = 40; // slider value that means "no decay"
const LOWER = matrix([[0, 1], [0, 0]]);
const SCREEN = Float64Array.from({ length: 257 }, (_, index) => -3.2 + (6.4 * index) / 256);

const $ = (id) => document.getElementById(id);
const renderer = createBlochRenderer($('stage'), { labels: LABELS });
renderer.overlay.setAttribute('aria-label', 'Bloch sphere, beam map and analyser screen computed by the engine');

const MODES = {
  qubit: {
    explain: 'Rotating frame: H = (Δ/2)σz + (Ω/2)σx, with decay |1⟩→|0⟩ at rate 1/T1 and pure dephasing at rate 1/Tφ. The ensemble is propagated exactly with exp(𝓛·dt). The screen shows a read-out in the computational basis.',
    duration: 24,
    controls: [
      { id: 'omega', label: 'Rabi frequency Ω', min: 0.2, max: 4, step: 0.05, value: 1.5 },
      { id: 'delta', label: 'Detuning Δ', min: -4, max: 4, step: 0.05, value: 0 },
      { id: 't1', label: 'Relaxation time T1', min: 1, max: FOREVER, step: 0.5, value: FOREVER, forever: true },
      { id: 'tphi', label: 'Pure dephasing time Tφ', min: 1, max: FOREVER, step: 0.5, value: FOREVER, forever: true }
    ],
    setup(p) {
      const dt = 0.03;
      const hamiltonian = add(scale(PAULI.Z, p.delta / 2), scale(PAULI.X, p.omega / 2));
      const collapse = [];
      if (p.t1 < FOREVER) collapse.push({ operator: LOWER, rate: 1 / p.t1 });
      if (p.tphi < FOREVER) collapse.push({ operator: PAULI.Z, rate: 1 / (2 * p.tphi) });
      const rho = blochToDensity([0, 0, 1]);
      return { dt, rho, propagator: superPropagator(liouvillian(hamiltonian, collapse), dt), unitary: collapse.length === 0, trail: [blochVector(rho)] };
    },
    step(state) {
      state.rho = applySuper(state.propagator, state.rho);
      state.trail.push(blochVector(state.rho));
      if (state.trail.length > 900) state.trail.shift();
    },
    frame(state, p) {
      const vector = blochVector(state.rho);
      const up = 0.5 * (1 + vector[2]);
      const density = Float64Array.from(SCREEN, (z) => (up * Math.exp(-((z - 1.5) ** 2) / 0.2312) + (1 - up) * Math.exp(-((z + 1.5) ** 2) / 0.2312)) / (0.34 * Math.sqrt(2 * Math.PI)));
      return { x: SCREEN, density, view: [-3.2, 3.2], screenLabels: ['|0⟩', '|1⟩'], bloch: { vector, trail: state.trail, axis: [p.omega, 0, p.delta], measure: [0, 0, 1] } };
    },
    readouts(state, p, time) {
      const [x, y, z] = blochVector(state.rho);
      const population = 0.5 * (1 - z);
      const generalized = Math.hypot(p.omega, p.delta);
      const rows = [['time', time.toFixed(2)], ['P(|1⟩)', population.toFixed(5)], ['|⟨σ⟩|', Math.hypot(x, y, z).toFixed(5)], ['purity Tr ρ²', purity(state.rho).toFixed(5)]];
      rows.push(check('trace = 1', Math.abs(trace(state.rho)[0] - 1), 1e-10));
      if (state.unitary) {
        rows.push(check('Rabi law (Ω/Ω′)²·sin²(Ω′t/2)', Math.abs(population - (p.omega / generalized) ** 2 * Math.sin((generalized * time) / 2) ** 2), 1e-8));
        rows.push(check('pure state: |⟨σ⟩| = 1', Math.abs(Math.hypot(x, y, z) - 1), 1e-9));
      } else {
        rows.push(['T2 = 1/(1/2T1 + 1/Tφ)', (1 / ((p.t1 < FOREVER ? 1 / (2 * p.t1) : 0) + (p.tphi < FOREVER ? 1 / p.tphi : 0))).toFixed(3)]);
        rows.push(check('mixed state: |⟨σ⟩| ≤ 1', Math.max(0, Math.hypot(x, y, z) - 1), 1e-9));
      }
      return rows;
    }
  },
  sg: {
    explain: 'One dimension along the field gradient: H = p²/2m − f·z·σz. Each spin component is propagated with a Fourier split-step. As the two beams stop overlapping, the spin alone becomes a mixture: the arrow shrinks although nothing has been measured.',
    duration: 12,
    controls: [
      { id: 'force', label: 'Gradient force f', min: 0.05, max: 0.5, step: 0.01, value: 0.2 },
      { id: 'tilt', label: 'Initial spin tilt θ (degrees)', min: 0, max: 180, step: 1, value: 90 }
    ],
    setup(p) {
      const grid = createGrid({ size: 1024, length: 200 });
      const angle = (p.tilt * Math.PI) / 180;
      return {
        dt: 0.025,
        grid,
        spinor: createSpinor(grid, gaussianPacket(grid, { center: 0, sigma: 1.5 }), [[Math.cos(angle / 2), 0], [Math.sin(angle / 2), 0]]),
        up: Float64Array.from(grid.x, (z) => -p.force * z),
        down: Float64Array.from(grid.x, (z) => p.force * z),
        propagator: createSpinorPropagator(grid, { mass: 1 }),
        rows: [],
        times: [],
        count: 0
      };
    },
    step(state, time) {
      if (state.count % 4 === 0) record(state, time);
      state.propagator.step(state.spinor, state.up, state.down, state.dt);
      state.count += 1;
    },
    frame(state) {
      const density = spinorDensity(state.grid, state.spinor);
      return {
        x: state.grid.x,
        density: density.total,
        view: [-50, 50],
        ghosts: [{ values: density.up }, { values: density.down }],
        bloch: { vector: blochVector(reducedSpin(state.grid, state.spinor)), measure: [0, 0, 1] },
        carpet: { columns: 200, rows: state.rows, times: state.times }
      };
    },
    readouts(state, p, time) {
      const upper = componentMoments(state.grid, state.spinor.up);
      const lower = componentMoments(state.grid, state.spinor.down);
      const spin = reducedSpin(state.grid, state.spinor);
      const angle = (p.tilt * Math.PI) / 180;
      const rows = [['time', time.toFixed(2)], ['weight of the ↑ beam', upper.weight.toFixed(5)], ['|⟨σ⟩| of the spin alone', Math.hypot(...blochVector(spin)).toFixed(5)]];
      rows.push(check('norm = 1', Math.abs(upper.weight + lower.weight - 1), 1e-8));
      rows.push(check('weight ↑ = cos²(θ/2)', Math.abs(upper.weight - Math.cos(angle / 2) ** 2), 1e-8));
      if (upper.weight > 1e-6 && lower.weight > 1e-6) {
        rows.push(['beam separation', (upper.mean - lower.mean).toFixed(3)]);
        rows.push(check('separation = f·t²/m', Math.abs(upper.mean - lower.mean - p.force * time * time), 5e-3));
        const overlap = Math.hypot(spin.re[1], spin.im[1]) / Math.sqrt(upper.weight * lower.weight);
        const law = Math.exp(-(((2 * p.force * time) ** 2 * 2.25) / 2 + (p.force * time * time) ** 2 / 18));
        rows.push(check('spin coherence = beam overlap', Math.abs(overlap - law), 1e-5));
      }
      return rows;
    }
  }
};

// One row of the beam map: |Ψ(z)|² summed over spin, sampled across the view.
function record(state, time) {
  const density = spinorDensity(state.grid, state.spinor).total;
  const row = new Float32Array(200);
  const origin = state.grid.x[0];
  const cell = state.grid.x[1] - state.grid.x[0];
  for (let column = 0; column < 200; column += 1) {
    const index = Math.round((-50 + ((column + 0.5) / 200) * 100 - origin) / cell);
    row[column] = density[index] || 0;
  }
  state.rows.push(row);
  state.times.push(time);
}

function check(label, error, tolerance) {
  return [label, `${error <= tolerance ? '✓' : '✗'} ${error.toExponential(1)}`, error <= tolerance ? 'ok' : 'bad'];
}

const app = { mode: 'qubit', parameters: {}, state: null, time: 0, playing: false, events: null, last: null, detections: 0 };

function buildControls() {
  const mode = MODES[app.mode];
  const host = $('controls');
  host.replaceChildren();
  app.parameters = {};
  for (const control of mode.controls) {
    app.parameters[control.id] = control.value;
    const label = document.createElement('label');
    label.className = 'field';
    const value = document.createElement('span');
    const input = document.createElement('input');
    input.type = 'range';
    input.min = control.min;
    input.max = control.max;
    input.step = control.step;
    input.value = control.value;
    input.id = `ctl-${control.id}`;
    const show = () => { value.textContent = control.forever && Number(input.value) >= FOREVER ? '∞' : Number(input.value).toFixed(control.step < 1 ? 2 : 0); };
    input.addEventListener('input', () => {
      app.parameters[control.id] = Number(input.value);
      show();
      restart();
    });
    show();
    label.append(control.label, value, input);
    host.append(label);
  }
  $('explain').textContent = mode.explain;
}

function restart() {
  const mode = MODES[app.mode];
  app.state = mode.setup(app.parameters);
  app.time = 0;
  app.events = null;
  app.detections = 0;
  if (app.mode === 'sg') record(app.state, 0);
  renderer.resetScale();
  draw();
}

function draw() {
  const mode = MODES[app.mode];
  const frame = mode.frame(app.state, app.parameters);
  frame.events = app.events;
  renderer.render(frame, { strict: $('strict').checked, duration: mode.duration });
  const list = $('readouts');
  list.replaceChildren();
  for (const [label, value, state] of mode.readouts(app.state, app.parameters, app.time)) {
    const term = document.createElement('dt');
    term.textContent = label;
    const detail = document.createElement('dd');
    detail.textContent = value;
    if (state) detail.className = state;
    list.append(term, detail);
  }
}

function tick(now) {
  if (app.playing) {
    const mode = MODES[app.mode];
    const elapsed = app.last === null ? 0.016 : Math.min(0.05, (now - app.last) / 1000);
    // Two time units per second of wall clock, in whole solver steps.
    let steps = Math.max(1, Math.round((elapsed * 2) / app.state.dt));
    while (steps > 0 && app.time < mode.duration - app.state.dt / 2) {
      mode.step(app.state, app.time);
      app.time += app.state.dt;
      steps -= 1;
    }
    if (app.time >= mode.duration - app.state.dt / 2) setPlaying(false);
    draw();
  }
  app.last = now;
  requestAnimationFrame(tick);
}

function setPlaying(value) {
  app.playing = value;
  $('play').textContent = value ? 'Pause' : 'Play';
  $('play').setAttribute('aria-pressed', String(value));
}

// 200 detections sampled from the density on the screen (Born rule) with the seed in the box.
function detect() {
  const frame = MODES[app.mode].frame(app.state, app.parameters);
  const cumulative = new Float64Array(frame.density.length);
  let total = 0;
  for (let index = 0; index < cumulative.length; index += 1) {
    total += Math.max(0, frame.density[index]);
    cumulative[index] = total;
  }
  app.detections += 1;
  const rng = createRng(`${$('seed').value}/${app.mode}/${app.detections}`);
  const positions = new Float64Array(200);
  const cell = frame.x[1] - frame.x[0];
  let above = 0;
  for (let shot = 0; shot < 200; shot += 1) {
    const target = rng.uniform() * total;
    let low = 0;
    let high = cumulative.length - 1;
    while (low < high) {
      const middle = (low + high) >> 1;
      if (cumulative[middle] < target) low = middle + 1;
      else high = middle;
    }
    positions[shot] = frame.x[low] + (rng.uniform() - 0.5) * cell;
    if (positions[shot] > 0) above += 1;
  }
  app.events = positions;
  $('status').textContent = `200 detections: ${above} above the centre line, ${200 - above} below.`;
  draw();
}

function setMode(name) {
  app.mode = name;
  $('mode-qubit').setAttribute('aria-pressed', String(name === 'qubit'));
  $('mode-sg').setAttribute('aria-pressed', String(name === 'sg'));
  $('status').textContent = '';
  setPlaying(false);
  buildControls();
  restart();
}

// Buttons are wired by delegation on the document.
document.addEventListener('click', (event) => {
  const id = event.target && event.target.id;
  if (id === 'mode-qubit') setMode('qubit');
  else if (id === 'mode-sg') setMode('sg');
  else if (id === 'play') {
    if (!app.playing && app.time >= MODES[app.mode].duration - app.state.dt) restart();
    setPlaying(!app.playing);
  } else if (id === 'reset') {
    setPlaying(false);
    $('status').textContent = '';
    restart();
  } else if (id === 'detect') detect();
});
$('strict').addEventListener('change', draw);
window.addEventListener('resize', draw);

setMode('qubit');
requestAnimationFrame(tick);
