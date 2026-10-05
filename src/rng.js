// Generador pseudoaleatorio reproducible por semilla textual.
// La misma semilla produce siempre la misma secuencia: es la base de toda medición repetible.

export function hashSeed(seed) {
  const text = String(seed);
  let hash = 2166136261;
  for (let index = 0; index < text.length; index += 1) {
    hash ^= text.charCodeAt(index);
    hash = Math.imul(hash, 16777619);
  }
  return hash >>> 0;
}

export function createRng(seed) {
  const label = String(seed);
  let state = hashSeed(label);
  let spare = null;

  function uniform() {
    state = (state + 0x6D2B79F5) | 0;
    let mixed = Math.imul(state ^ (state >>> 15), 1 | state);
    mixed = (mixed + Math.imul(mixed ^ (mixed >>> 7), 61 | mixed)) ^ mixed;
    return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
  }

  function gaussian() {
    if (spare !== null) {
      const value = spare;
      spare = null;
      return value;
    }
    let first = uniform();
    while (first <= 1e-300) first = uniform();
    const second = uniform();
    const radius = Math.sqrt(-2 * Math.log(first));
    const angle = 2 * Math.PI * second;
    spare = radius * Math.sin(angle);
    return radius * Math.cos(angle);
  }

  return {
    seed: label,
    uniform,
    gaussian,
    // Secuencia independiente y reproducible derivada de esta semilla.
    fork(name) { return createRng(`${label}/${name}`); }
  };
}
