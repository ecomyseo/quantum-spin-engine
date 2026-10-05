// The Bloch renderer draws only what the frame carries. A spy canvas records every stroke, so the test
// can check where the state arrow ends and that strict mode really removes the interpretive layers.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createBlochRenderer, project, LABELS } from '../src/bloch-renderer.js';

function environment() {
  const calls = { count: 0, nan: 0, text: [] };
  const context = () => new Proxy({}, {
    get(target, property) {
      if (property in target) return target[property];
      if (property === 'createImageData') return (width, height) => ({ width, height, data: new Uint8ClampedArray(width * height * 4) });
      if (property === 'fillText') return (label) => { calls.count += 1; calls.text.push(String(label)); };
      return (...args) => {
        calls.count += 1;
        for (const value of args) if (typeof value === 'number' && Number.isNaN(value)) calls.nan += 1;
      };
    },
    set(target, property, value) { target[property] = value; return true; }
  });
  const element = (tag) => {
    const node = { tagName: tag, children: [], style: {}, attributes: {}, clientWidth: 960, clientHeight: 420, width: 300, height: 150 };
    node.setAttribute = (name, value) => { node.attributes[name] = String(value); };
    node.appendChild = (child) => { node.children.push(child); return child; };
    if (tag === 'canvas') {
      let cached = null;
      node.getContext = () => { cached = cached || context(); return cached; };
    }
    return node;
  };
  return { calls, documentRef: { createElement: element }, windowRef: { devicePixelRatio: 1 }, host: element('div') };
}

const axis = Float64Array.from({ length: 65 }, (_, index) => -3.2 + 0.1 * index);
const baseFrame = { x: axis, density: Float64Array.from(axis, (value) => Math.exp(-((value - 1.5) ** 2) / 0.23)), view: [-3.2, 3.2] };

test('the projection is a rigid rotation with the poles on the vertical', () => {
  const north = project([0, 0, 1]);
  assert.ok(Math.abs(north[0]) < 1e-12 && north[1] > 0.9);
  for (const point of [[1, 0, 0], [0, 1, 0], [0.6, -0.48, 0.64]]) assert.ok(Math.abs(Math.hypot(...project(point)) - Math.hypot(...point)) < 1e-12);
  const a = project([1, 0, 0]);
  const b = project([0, 1, 0]);
  assert.ok(Math.abs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) < 1e-12);
});

test('the state arrow ends where the Bloch vector says, with its length', () => {
  const env = environment();
  const renderer = createBlochRenderer(env.host, { documentRef: env.documentRef, windowRef: env.windowRef });
  const context = renderer.overlay.getContext('2d');
  const strokes = [];
  let last = null;
  context.lineTo = (x, y) => { last = [x, y]; };
  context.stroke = () => { strokes.push({ width: context.lineWidth, to: last }); };
  const tipOf = (vector) => {
    strokes.length = 0;
    renderer.render({ ...baseFrame, bloch: { vector } }, {});
    const arrow = strokes.filter((stroke) => stroke.width === 2.6);
    assert.equal(arrow.length, 1);
    return arrow[0].to;
  };
  const centre = tipOf([0, 0, 0]);
  const north = tipOf([0, 0, 1]);
  const half = tipOf([0, 0, 0.5]);
  assert.ok(north[1] < centre[1] && tipOf([0, 0, -1])[1] > centre[1]);
  assert.ok(tipOf([0, 1, 0])[0] > centre[0]);
  assert.ok(Math.abs((half[1] - centre[1]) - 0.5 * (north[1] - centre[1])) < 1e-9);
});

test('strict mode removes trail, samples and reference state; labels can be replaced', () => {
  const frame = {
    ...baseFrame,
    events: Float64Array.of(-1.5, 1.4, 1.6),
    screenLabels: ['up', 'down'],
    bloch: { vector: [0.3, 0.2, 0.5], trail: [[0, 0, 1], [0.1, 0, 0.99], [0.3, 0.2, 0.5]], axis: [0, 0, 2], measure: [1, 0, 0], samples: [[0, 1, 0]], ghost: [0, 0, 1] },
    circuit: { elements: [{ label: 'H', detail: '', flow: 1, blocked: 0 }, { label: 'SG', detail: 'up', flow: 0.5, blocked: 0.5 }], cursor: 0.4 },
    carpet: { columns: 8, rows: [new Float32Array(8).fill(0.2), new Float32Array(8).fill(0.1)], times: [0, 1] }
  };
  const env = environment();
  const renderer = createBlochRenderer(env.host, { documentRef: env.documentRef, windowRef: env.windowRef, labels: { caption: ['Bloch sphere: a map of states'], time: 'time →', absorbed: '% absorbed', decimal: '.' } });
  assert.equal(renderer.overlay.attributes.role, 'img');
  renderer.render(frame, { strict: false, duration: 2 });
  const normal = env.calls.count;
  renderer.render(frame, { strict: true, duration: 2 });
  assert.ok(env.calls.count - normal < normal, 'strict mode draws fewer primitives');
  assert.equal(env.calls.nan, 0);
  assert.ok(env.calls.text.includes('Bloch sphere: a map of states'));
  assert.ok(env.calls.text.includes('time →'));
  assert.ok(env.calls.text.some((label) => label.endsWith('% absorbed') && label.startsWith('50')));
  assert.ok(env.calls.text.some((label) => label.startsWith('|⟨σ⟩| = 0.6')), 'decimal point instead of comma');
  assert.ok(!env.calls.text.includes(LABELS.time));
});
