// Static checks of the demo and of the package layout: every file the demo and the package refer to
// exists, the engine never touches the browser, and nothing draws from unseeded randomness.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const read = (file) => readFileSync(path.join(ROOT, file), 'utf8');

test('the demo page loads files that exist and its script only imports from src/', () => {
  const html = read('demo/index.html');
  assert.ok(html.includes('<html lang="en">') && html.includes('name="viewport"'));
  for (const [, href] of html.matchAll(/<link rel="stylesheet" href="([^"]+)"/g)) assert.ok(existsSync(path.join(ROOT, 'demo', href)), href);
  for (const [, source] of html.matchAll(/<script type="module" src="([^"]+)"/g)) assert.ok(existsSync(path.join(ROOT, 'demo', source)), source);
  const script = read('demo/main.js');
  const imports = [...script.matchAll(/from '([^']+)'/g)].map((match) => match[1]);
  assert.ok(imports.length >= 6);
  for (const target of imports) {
    assert.ok(target.startsWith('../src/'), target);
    assert.ok(existsSync(path.join(ROOT, 'demo', target)), target);
  }
  // Every element the script looks up by id is in the page.
  for (const [, id] of script.matchAll(/\$\('([a-z-]+)'\)/g)) assert.ok(html.includes(`id="${id}"`), `missing #${id}`);
  for (const id of ['mode-qubit', 'mode-sg', 'play', 'reset', 'detect']) assert.ok(script.includes(`'${id}'`) && html.includes(`id="${id}"`), id);
});

test('every import inside src/ resolves and the solver modules never touch the browser', () => {
  const files = readdirSync(path.join(ROOT, 'src')).filter((name) => name.endsWith('.js'));
  assert.ok(files.length >= 11);
  for (const name of files) {
    const source = read(`src/${name}`);
    for (const [, target] of source.matchAll(/from '([^']+)'/g)) {
      assert.ok(target.startsWith('./'), `${name} imports ${target}`);
      assert.ok(existsSync(path.join(ROOT, 'src', target)), `${name} imports a missing file: ${target}`);
    }
    assert.ok(!/Math\.random/.test(source), `${name} uses unseeded randomness`);
    if (name !== 'bloch-renderer.js') assert.ok(!/\b(document|window|requestAnimationFrame|getContext)\b/.test(source), `${name} touches the browser`);
  }
  assert.ok(!/Math\.random/.test(read('demo/main.js')));
});

test('package.json points to files that exist and the README documents the public names', async () => {
  const manifest = JSON.parse(read('package.json'));
  for (const target of Object.values(manifest.exports)) assert.ok(existsSync(path.join(ROOT, target)), target);
  for (const entry of manifest.files) assert.ok(existsSync(path.join(ROOT, entry)), entry);
  const engine = await import('../src/index.js');
  const readme = read('README.md');
  for (const name of ['expm', 'hermitianEigenvalues', 'blochVector', 'rotation', 'liouvillian', 'superPropagator', 'applySuper', 'createJumpUnraveling', 'measureAlong', 'tomography', 'createSpinor', 'reducedSpin', 'createBlochRenderer']) {
    assert.ok(readme.includes(name), `the README does not mention ${name}`);
    if (name !== 'createBlochRenderer') assert.equal(typeof engine[name], 'function', name);
  }
  assert.ok(readme.includes(manifest.homepage));
  for (const [, image] of readme.matchAll(/\]\((media\/[^)]+)\)/g)) assert.ok(existsSync(path.join(ROOT, image)), image);
  assert.ok(read('LICENSE').startsWith('MIT License'));
  assert.ok(read('CHANGELOG.md').includes(manifest.version));
});
