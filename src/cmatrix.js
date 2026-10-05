// Matrices complejas densas y pequeñas (hasta unas decenas de filas), sin dependencias.
// Una matriz es { n, re, im } con los elementos por filas. Es la base del motor de sistemas de pocos niveles:
// operadores, matrices densidad, superoperadores de Lindblad y exponenciales de matrices.

export function zeros(n) {
  return { n, re: new Float64Array(n * n), im: new Float64Array(n * n) };
}

export function identity(n) {
  const result = zeros(n);
  for (let index = 0; index < n; index += 1) result.re[index * n + index] = 1;
  return result;
}

// Construye una matriz a partir de filas; cada elemento es un número real o un par [real, imaginaria].
export function matrix(rows) {
  const n = rows.length;
  const result = zeros(n);
  for (let row = 0; row < n; row += 1) {
    if (rows[row].length !== n) throw new RangeError('La matriz tiene que ser cuadrada.');
    for (let column = 0; column < n; column += 1) {
      const value = rows[row][column];
      result.re[row * n + column] = Array.isArray(value) ? value[0] : value;
      result.im[row * n + column] = Array.isArray(value) ? value[1] : 0;
    }
  }
  return result;
}

export function clone(a) {
  return { n: a.n, re: Float64Array.from(a.re), im: Float64Array.from(a.im) };
}

export function add(a, b) {
  const result = zeros(a.n);
  for (let index = 0; index < a.re.length; index += 1) {
    result.re[index] = a.re[index] + b.re[index];
    result.im[index] = a.im[index] + b.im[index];
  }
  return result;
}

export function subtract(a, b) {
  const result = zeros(a.n);
  for (let index = 0; index < a.re.length; index += 1) {
    result.re[index] = a.re[index] - b.re[index];
    result.im[index] = a.im[index] - b.im[index];
  }
  return result;
}

// Producto por el escalar complejo (real + i·imaginary).
export function scale(a, real, imaginary = 0) {
  const result = zeros(a.n);
  for (let index = 0; index < a.re.length; index += 1) {
    result.re[index] = a.re[index] * real - a.im[index] * imaginary;
    result.im[index] = a.re[index] * imaginary + a.im[index] * real;
  }
  return result;
}

export function multiply(a, b) {
  const n = a.n;
  const result = zeros(n);
  for (let row = 0; row < n; row += 1) {
    for (let inner = 0; inner < n; inner += 1) {
      const ar = a.re[row * n + inner];
      const ai = a.im[row * n + inner];
      if (ar === 0 && ai === 0) continue;
      for (let column = 0; column < n; column += 1) {
        const br = b.re[inner * n + column];
        const bi = b.im[inner * n + column];
        result.re[row * n + column] += ar * br - ai * bi;
        result.im[row * n + column] += ar * bi + ai * br;
      }
    }
  }
  return result;
}

export function dagger(a) {
  const n = a.n;
  const result = zeros(n);
  for (let row = 0; row < n; row += 1) {
    for (let column = 0; column < n; column += 1) {
      result.re[column * n + row] = a.re[row * n + column];
      result.im[column * n + row] = -a.im[row * n + column];
    }
  }
  return result;
}

export function conjugate(a) {
  return { n: a.n, re: Float64Array.from(a.re), im: Float64Array.from(a.im, (value) => -value) };
}

export function trace(a) {
  let real = 0;
  let imaginary = 0;
  for (let index = 0; index < a.n; index += 1) {
    real += a.re[index * a.n + index];
    imaginary += a.im[index * a.n + index];
  }
  return [real, imaginary];
}

// Producto de Kronecker a ⊗ b.
export function kron(a, b) {
  const n = a.n * b.n;
  const result = zeros(n);
  for (let ar = 0; ar < a.n; ar += 1) {
    for (let ac = 0; ac < a.n; ac += 1) {
      const xr = a.re[ar * a.n + ac];
      const xi = a.im[ar * a.n + ac];
      if (xr === 0 && xi === 0) continue;
      for (let br = 0; br < b.n; br += 1) {
        for (let bc = 0; bc < b.n; bc += 1) {
          const yr = b.re[br * b.n + bc];
          const yi = b.im[br * b.n + bc];
          const target = (ar * b.n + br) * n + (ac * b.n + bc);
          result.re[target] = xr * yr - xi * yi;
          result.im[target] = xr * yi + xi * yr;
        }
      }
    }
  }
  return result;
}

export function commutator(a, b) {
  return subtract(multiply(a, b), multiply(b, a));
}

// Norma 1: la mayor suma de módulos por columna.
export function normOne(a) {
  let largest = 0;
  for (let column = 0; column < a.n; column += 1) {
    let total = 0;
    for (let row = 0; row < a.n; row += 1) total += Math.hypot(a.re[row * a.n + column], a.im[row * a.n + column]);
    largest = Math.max(largest, total);
  }
  return largest;
}

// Mayor diferencia, elemento a elemento, entre dos matrices.
export function distance(a, b) {
  let largest = 0;
  for (let index = 0; index < a.re.length; index += 1) largest = Math.max(largest, Math.hypot(a.re[index] - b.re[index], a.im[index] - b.im[index]));
  return largest;
}

// Exponencial de una matriz cualquiera: se reduce a la mitad hasta que su norma es pequeña, se suma la serie
// de Taylor y se eleva al cuadrado las veces que se redujo.
export function expm(a) {
  const norm = normOne(a);
  const squarings = norm > 0.5 ? Math.ceil(Math.log2(norm / 0.5)) : 0;
  const scaled = scale(a, 1 / 2 ** squarings);
  let result = identity(a.n);
  let term = identity(a.n);
  for (let order = 1; order <= 18; order += 1) {
    term = scale(multiply(term, scaled), 1 / order);
    result = add(result, term);
  }
  for (let count = 0; count < squarings; count += 1) result = multiply(result, result);
  return result;
}

// Aplica la matriz a un vector { re, im }.
export function applyToVector(a, vector) {
  const n = a.n;
  const result = { re: new Float64Array(n), im: new Float64Array(n) };
  for (let row = 0; row < n; row += 1) {
    let real = 0;
    let imaginary = 0;
    for (let column = 0; column < n; column += 1) {
      const mr = a.re[row * n + column];
      const mi = a.im[row * n + column];
      real += mr * vector.re[column] - mi * vector.im[column];
      imaginary += mr * vector.im[column] + mi * vector.re[column];
    }
    result.re[row] = real;
    result.im[row] = imaginary;
  }
  return result;
}

// Autovalores de una matriz hermítica por rotaciones de Jacobi complejas. Devuelve los valores ordenados.
export function hermitianEigenvalues(a) {
  const n = a.n;
  const work = clone(a);
  for (let sweep = 0; sweep < 60; sweep += 1) {
    let off = 0;
    for (let p = 0; p < n; p += 1) for (let q = p + 1; q < n; q += 1) off += work.re[p * n + q] ** 2 + work.im[p * n + q] ** 2;
    if (off < 1e-30) break;
    for (let p = 0; p < n; p += 1) {
      for (let q = p + 1; q < n; q += 1) {
        const xr = work.re[p * n + q];
        const xi = work.im[p * n + q];
        const modulus = Math.hypot(xr, xi);
        if (modulus < 1e-300) continue;
        // Fase que hace real el elemento (p,q) y ángulo que lo anula.
        const pr = xr / modulus;
        const pi = xi / modulus;
        const app = work.re[p * n + p];
        const aqq = work.re[q * n + q];
        const theta = 0.5 * Math.atan2(2 * modulus, app - aqq);
        const c = Math.cos(theta);
        const s = Math.sin(theta);
        // Columnas: A ← A·G, con G = [[c, −s·e^{iφ}], [s·e^{−iφ}, c]] en el plano (p,q).
        for (let k = 0; k < n; k += 1) {
          const kpR = work.re[k * n + p];
          const kpI = work.im[k * n + p];
          const kqR = work.re[k * n + q];
          const kqI = work.im[k * n + q];
          // (k,q)·e^{−iφ}
          const rotR = kqR * pr + kqI * pi;
          const rotI = kqI * pr - kqR * pi;
          work.re[k * n + p] = c * kpR + s * rotR;
          work.im[k * n + p] = c * kpI + s * rotI;
          // −s·e^{iφ}·(k,p) + c·(k,q)
          work.re[k * n + q] = -s * (kpR * pr - kpI * pi) + c * kqR;
          work.im[k * n + q] = -s * (kpR * pi + kpI * pr) + c * kqI;
        }
        // Filas: A ← G†·A.
        for (let k = 0; k < n; k += 1) {
          const pkR = work.re[p * n + k];
          const pkI = work.im[p * n + k];
          const qkR = work.re[q * n + k];
          const qkI = work.im[q * n + k];
          // e^{iφ}·(q,k)
          const rotR = qkR * pr - qkI * pi;
          const rotI = qkR * pi + qkI * pr;
          work.re[p * n + k] = c * pkR + s * rotR;
          work.im[p * n + k] = c * pkI + s * rotI;
          // −s·e^{−iφ}·(p,k) + c·(q,k)
          work.re[q * n + k] = -s * (pkR * pr + pkI * pi) + c * qkR;
          work.im[q * n + k] = -s * (pkI * pr - pkR * pi) + c * qkI;
        }
      }
    }
  }
  const values = [];
  for (let index = 0; index < n; index += 1) values.push(work.re[index * n + index]);
  return values.sort((left, right) => left - right);
}
