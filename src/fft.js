// Transformada rápida de Fourier compleja, radix-2 y en sitio.
// Convenio: forward calcula X_k = Σ x_n · exp(−2πi·k·n/N); inverse divide por N.

export function createFft(size) {
  if (!Number.isInteger(size) || size < 2 || (size & (size - 1)) !== 0) {
    throw new RangeError('El tamaño de la transformada debe ser una potencia de dos.');
  }

  const half = size >> 1;
  const cosTable = new Float64Array(half);
  const sinTable = new Float64Array(half);
  for (let index = 0; index < half; index += 1) {
    const angle = (2 * Math.PI * index) / size;
    cosTable[index] = Math.cos(angle);
    sinTable[index] = Math.sin(angle);
  }

  const bits = Math.round(Math.log2(size));
  const reverse = new Uint32Array(size);
  for (let index = 0; index < size; index += 1) {
    let reversed = 0;
    for (let bit = 0; bit < bits; bit += 1) {
      reversed |= ((index >>> bit) & 1) << (bits - 1 - bit);
    }
    reverse[index] = reversed >>> 0;
  }

  function transform(re, im, inverse) {
    for (let index = 0; index < size; index += 1) {
      const target = reverse[index];
      if (target > index) {
        const realValue = re[index];
        re[index] = re[target];
        re[target] = realValue;
        const imaginaryValue = im[index];
        im[index] = im[target];
        im[target] = imaginaryValue;
      }
    }

    for (let span = 1; span < size; span <<= 1) {
      const tableStep = size / (span << 1);
      for (let start = 0; start < size; start += span << 1) {
        for (let offset = 0, tableIndex = 0; offset < span; offset += 1, tableIndex += tableStep) {
          const wr = cosTable[tableIndex];
          const wi = inverse ? sinTable[tableIndex] : -sinTable[tableIndex];
          const low = start + offset;
          const high = low + span;
          const tr = re[high] * wr - im[high] * wi;
          const ti = re[high] * wi + im[high] * wr;
          re[high] = re[low] - tr;
          im[high] = im[low] - ti;
          re[low] += tr;
          im[low] += ti;
        }
      }
    }

    if (inverse) {
      const scale = 1 / size;
      for (let index = 0; index < size; index += 1) {
        re[index] *= scale;
        im[index] *= scale;
      }
    }
  }

  return Object.freeze({
    size,
    forward(re, im) { transform(re, im, false); },
    inverse(re, im) { transform(re, im, true); }
  });
}
