// Visor de espín, Stern–Gerlach y qubits sobre Canvas 2D.
// Dibuja exclusivamente lo que entrega el solver en frame(); no calcula física ni añade adornos.
//   · esfera de Bloch: el vector ⟨σ⟩ del estado, su rastro, el eje de giro y el eje de medida. Es un mapa
//     del espacio de estados, no un objeto que gira: por eso no lleva textura ni animación propia;
//   · mapa de haces (si el solver lo trae): |ψ(z,t)|² con el tiempo en horizontal y z en vertical;
//   · pantalla del analizador: perfil de probabilidad sobre el eje de medida y detecciones muestreadas;
//   · tira de circuito: puertas, pulsos o analizadores en orden, con la fracción del haz que sigue.

// Colores y formato de números del visor.
const PALETTE = Object.freeze({
  ink: '#f4efe6',
  muted: '#b7b0a5',
  quiet: '#817a71',
  gold: '#e8b56d',
  detector: '#ffd88a',
  reference: '#ff9d7a'
});

function formatTick(value) {
  const rounded = Math.abs(value) >= 10 ? value.toFixed(0) : Number(value.toFixed(2)).toString();
  return rounded.replace('.', ',').replace('-', '\u2212');
}

const CAMERA = Object.freeze({ azimuth: (28 * Math.PI) / 180, elevation: (17 * Math.PI) / 180 });
const BEAM = 'rgba(127, 214, 255, 0.9)';
// Rótulos por defecto. Quien use el visor fuera de este sitio puede pasar los suyos con «labels».
export const LABELS = Object.freeze({
  north: '|0⟩ · ↑z',
  south: '|1⟩ · ↓z',
  caption: ['Esfera de Bloch: mapa de estados, no un objeto que gira'],
  captionNarrow: ['Mapa de estados,', 'no un objeto que gira'],
  time: 'tiempo →',
  height: 'z ↑',
  absorbed: '% absorbido',
  decimal: ','
});

// Proyección ortográfica de un punto de la esfera: devuelve [derecha, arriba, hacia el observador].
export function project([x, y, z]) {
  const toward = x * Math.cos(CAMERA.azimuth) + y * Math.sin(CAMERA.azimuth);
  const right = -x * Math.sin(CAMERA.azimuth) + y * Math.cos(CAMERA.azimuth);
  return [right, z * Math.cos(CAMERA.elevation) - toward * Math.sin(CAMERA.elevation), toward * Math.cos(CAMERA.elevation) + z * Math.sin(CAMERA.elevation)];
}

// Desplazamiento horizontal reproducible para cada detección: no es un dato, solo evita que se tapen.
function spreadOf(index) {
  const value = Math.sin((index + 1) * 12.9898) * 43758.5453;
  return value - Math.floor(value);
}

export function createSpinRenderer(container, { documentRef = document, windowRef = window, labels = {} } = {}) {
  const t = { ...LABELS, ...labels };
  const number = (value) => (t.decimal === ',' ? formatTick(value) : formatTick(value).replace(',', '.').replace('−', '-'));
  const stage = documentRef.createElement('div');
  stage.className = 'fq-stage';
  const overlay = documentRef.createElement('canvas');
  overlay.className = 'fq-overlay';
  overlay.setAttribute('role', 'img');
  overlay.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;display:block';
  stage.appendChild(overlay);
  container.appendChild(stage);
  const context = overlay.getContext('2d');
  let mapCanvas = null;

  function resize() {
    const ratio = Math.min(windowRef.devicePixelRatio || 1, 2);
    const width = Math.max(320, Math.round((stage.clientWidth || 960) * ratio));
    const height = Math.max(200, Math.round((stage.clientHeight || 400) * ratio));
    if (overlay.width !== width || overlay.height !== height) {
      overlay.width = width;
      overlay.height = height;
    }
    return { width, height, ratio };
  }

  function text(label, x, y, { color = PALETTE.muted, size = 11, align = 'left', ratio = 1, weight = 500 } = {}) {
    context.fillStyle = color;
    context.font = `${weight} ${Math.round(size * ratio)}px Inter, system-ui, sans-serif`;
    context.textAlign = align;
    context.textBaseline = 'middle';
    context.fillText(label, x, y);
  }

  function line(points, { color, width = 1, dash = null }) {
    if (points.length < 2) return;
    context.beginPath();
    context.moveTo(points[0][0], points[0][1]);
    for (let index = 1; index < points.length; index += 1) context.lineTo(points[index][0], points[index][1]);
    context.strokeStyle = color;
    context.lineWidth = width;
    context.setLineDash(dash || []);
    context.stroke();
    context.setLineDash([]);
  }

  function drawSphere(bloch, box, ratio, strict) {
    // En pantallas estrechas el rótulo va en dos líneas para que no se corte.
    const caption = box.width >= 300 * ratio ? t.caption : t.captionNarrow;
    // El pie (módulo del vector y rótulo) ocupa una franja abajo y los rótulos de los polos sobresalen de
    // la esfera. En un lienzo bajo —un teléfono tumbado— la esfera se encoge y se centra en el hueco que
    // queda, para que el polo sur no se pise con el pie ni el polo norte se salga por arriba.
    const footer = (33 + 12 * caption.length) * ratio;
    const top = 9 * ratio;
    const reach = 1.34 * project([0, 0, 1])[1];
    const usable = box.height - footer - top;
    const radius = Math.min(Math.min(box.width, box.height) * 0.36, usable / (2 * reach));
    const cx = box.x + box.width / 2;
    const cy = box.y + Math.min(box.height * 0.47, top + usable / 2);
    const at = (point) => {
      const [right, up, toward] = project(point);
      return [cx + radius * right, cy - radius * up, toward];
    };
    // Contorno y círculos máximos: la mitad trasera, más tenue.
    context.beginPath();
    context.arc(cx, cy, radius, 0, 2 * Math.PI);
    context.fillStyle = 'rgba(232, 181, 109, 0.045)';
    context.fill();
    context.strokeStyle = 'rgba(244, 239, 230, 0.38)';
    context.lineWidth = 1.2 * ratio;
    context.stroke();
    const circles = [(angle) => [Math.cos(angle), Math.sin(angle), 0], (angle) => [Math.cos(angle), 0, Math.sin(angle)], (angle) => [0, Math.cos(angle), Math.sin(angle)]];
    for (const circle of circles) {
      for (const front of [false, true]) {
        let run = [];
        for (let step = 0; step <= 96; step += 1) {
          const [px, py, depth] = at(circle((step / 96) * 2 * Math.PI));
          if (depth >= 0 === front) run.push([px, py]);
          else {
            line(run, { color: front ? 'rgba(244, 239, 230, 0.26)' : 'rgba(244, 239, 230, 0.09)', width: ratio });
            run = [];
          }
        }
        line(run, { color: front ? 'rgba(244, 239, 230, 0.26)' : 'rgba(244, 239, 230, 0.09)', width: ratio });
      }
    }
    // Ejes y rótulos de los polos.
    for (const [axis, label] of [[[1, 0, 0], 'x'], [[0, 1, 0], 'y'], [[0, 0, 1], 'z']]) {
      const from = at(axis.map((value) => -value));
      const to = at(axis);
      line([from, to], { color: 'rgba(244, 239, 230, 0.16)', width: ratio });
      const tip = at(axis.map((value) => value * 1.16));
      text(label, tip[0], tip[1], { color: PALETTE.quiet, size: 11, align: 'center', ratio });
    }
    const north = at([0, 0, 1.34]);
    const south = at([0, 0, -1.34]);
    text(t.north, north[0], north[1], { color: PALETTE.muted, size: 11, align: 'center', ratio, weight: 700 });
    text(t.south, south[0], south[1], { color: PALETTE.muted, size: 11, align: 'center', ratio, weight: 700 });

    // Eje de giro (campo, pulso o puerta) y eje de medida.
    if (bloch.axis) {
      const norm = Math.hypot(...bloch.axis) || 1;
      const unit = bloch.axis.map((value) => value / norm);
      line([at(unit.map((value) => -1.08 * value)), at(unit.map((value) => 1.08 * value))], { color: 'rgba(232, 181, 109, 0.75)', width: 1.4 * ratio, dash: [6 * ratio, 5 * ratio] });
    }
    if (bloch.measure) {
      const tipUp = at(bloch.measure);
      const tipDown = at(bloch.measure.map((value) => -value));
      line([tipDown, tipUp], { color: 'rgba(127, 214, 255, 0.5)', width: 1.2 * ratio });
      for (const tip of [tipUp, tipDown]) {
        context.beginPath();
        context.arc(tip[0], tip[1], 3.2 * ratio, 0, 2 * Math.PI);
        context.strokeStyle = BEAM;
        context.lineWidth = 1.2 * ratio;
        context.stroke();
      }
    }
    if (!strict) {
      // Rastro del vector y, si las hay, realizaciones individuales.
      const trail = bloch.trail || [];
      for (let index = 1; index < trail.length; index += 1) {
        const from = at(trail[index - 1]);
        const to = at(trail[index]);
        const age = index / trail.length;
        line([from, to], { color: `rgba(232, 181, 109, ${(0.12 + 0.6 * age) * (to[2] >= 0 ? 1 : 0.4)})`, width: 1.6 * ratio });
      }
      for (const sample of bloch.samples || []) {
        const [px, py, depth] = at(sample);
        context.beginPath();
        context.arc(px, py, 2.2 * ratio, 0, 2 * Math.PI);
        context.fillStyle = depth >= 0 ? 'rgba(255, 216, 138, 0.8)' : 'rgba(255, 216, 138, 0.3)';
        context.fill();
      }
      if (bloch.ghost) line([[cx, cy], at(bloch.ghost)], { color: 'rgba(244, 239, 230, 0.4)', width: 1.4 * ratio, dash: [3 * ratio, 4 * ratio] });
    }
    // Vector del estado: su longitud es |⟨σ⟩|; menor que 1 si el estado es una mezcla.
    const tip = at(bloch.vector);
    line([[cx, cy], tip], { color: PALETTE.gold, width: 2.6 * ratio });
    context.beginPath();
    context.arc(tip[0], tip[1], 5.5 * ratio, 0, 2 * Math.PI);
    context.fillStyle = tip[2] >= 0 ? '#ffe1a6' : PALETTE.gold;
    context.fill();
    const length = Math.hypot(...bloch.vector);
    text(`|⟨σ⟩| = ${number(Math.round(length * 1000) / 1000)}`, cx, box.y + box.height - (19 + 12 * caption.length) * ratio, { color: PALETTE.ink, size: 12, align: 'center', ratio, weight: 700 });
    caption.forEach((line, index) => text(line, cx, box.y + box.height - (12 * (caption.length - index) + 1) * ratio, { color: PALETTE.quiet, size: 10, align: 'center', ratio }));
  }

  function drawMap(carpet, box, ratio, duration) {
    const rows = carpet.rows;
    const columns = carpet.columns;
    if (!mapCanvas) mapCanvas = documentRef.createElement('canvas');
    mapCanvas.width = rows.length;
    mapCanvas.height = columns;
    const mapContext = mapCanvas.getContext('2d');
    const image = mapContext.createImageData(rows.length, columns);
    let peak = 1e-12;
    for (const row of rows) for (let column = 0; column < columns; column += 1) if (row[column] > peak) peak = row[column];
    for (let time = 0; time < rows.length; time += 1) {
      for (let column = 0; column < columns; column += 1) {
        // z crece hacia arriba; brillo con raíz para que los haces ya ensanchados sigan viéndose.
        const level = Math.sqrt(Math.max(0, rows[time][column]) / peak);
        const offset = ((columns - 1 - column) * rows.length + time) * 4;
        image.data[offset] = 20 + 235 * level;
        image.data[offset + 1] = 26 + 190 * level;
        image.data[offset + 2] = 44 + 120 * level * level;
        image.data[offset + 3] = 255;
      }
    }
    mapContext.putImageData(image, 0, 0);
    const reached = duration > 0 ? Math.min(1, (carpet.times.at(-1) || 0) / duration) : 1;
    context.fillStyle = '#05070c';
    context.fillRect(box.x, box.y, box.width, box.height);
    context.imageSmoothingEnabled = true;
    if (reached > 0 && rows.length > 1) context.drawImage(mapCanvas, box.x, box.y, Math.max(1, box.width * reached), box.height);
    context.strokeStyle = 'rgba(244, 239, 230, 0.18)';
    context.lineWidth = ratio;
    context.strokeRect(box.x, box.y, box.width, box.height);
    text(t.time, box.x + box.width - 6 * ratio, box.y + box.height - 10 * ratio, { color: PALETTE.quiet, size: 10, align: 'right', ratio });
    text(t.height, box.x + 6 * ratio, box.y + 11 * ratio, { color: PALETTE.quiet, size: 10, ratio });
  }

  function drawScreen(frame, box, ratio, scaleState) {
    const [from, to] = frame.view;
    const count = frame.x.length;
    const origin = frame.x[0];
    const cell = frame.x[1] - frame.x[0];
    const yOf = (value) => box.y + box.height * (1 - (value - from) / (to - from));
    let peak = 1e-12;
    for (let index = 0; index < count; index += 1) if (frame.density[index] > peak) peak = frame.density[index];
    // La escala solo crece: evita que el perfil «respire» por un cambio de máximo entre fotogramas.
    scaleState.peak = Math.max(scaleState.peak * 0.985, peak);
    const profile = box.width * 0.46;
    const baseX = box.x + 8 * ratio;
    context.fillStyle = 'rgba(3, 5, 11, 0.85)';
    context.fillRect(box.x, box.y, box.width, box.height);
    context.strokeStyle = 'rgba(244, 239, 230, 0.18)';
    context.lineWidth = ratio;
    context.strokeRect(box.x, box.y, box.width, box.height);
    const draw = (values, fill, stroke) => {
      context.beginPath();
      context.moveTo(baseX, yOf(from));
      for (let pixel = 0; pixel <= box.height; pixel += 1) {
        const value = to - ((to - from) * pixel) / box.height;
        const index = Math.min(count - 1, Math.max(0, Math.round((value - origin) / cell)));
        context.lineTo(baseX + profile * Math.min(1, Math.max(0, values[index]) / scaleState.peak), box.y + pixel);
      }
      context.lineTo(baseX, box.y + box.height);
      context.closePath();
      if (fill) {
        context.fillStyle = fill;
        context.fill();
      }
      context.strokeStyle = stroke;
      context.lineWidth = 1.4 * ratio;
      context.stroke();
    };
    draw(frame.density, 'rgba(232, 181, 109, 0.28)', PALETTE.gold);
    for (const [order, ghost] of (frame.ghosts || []).entries()) draw(ghost.values, null, order === 0 ? 'rgba(255, 216, 138, 0.75)' : BEAM);
    // Línea de z = 0: separa las dos manchas.
    const middle = yOf(0);
    if (middle > box.y && middle < box.y + box.height) line([[box.x, middle], [box.x + box.width, middle]], { color: 'rgba(244, 239, 230, 0.14)', width: ratio, dash: [3 * ratio, 4 * ratio] });
    const events = frame.events;
    if (events && events.length) {
      const field = box.width - profile - 20 * ratio;
      context.fillStyle = PALETTE.detector;
      const limit = Math.min(events.length, 4000);
      for (let index = events.length - limit; index < events.length; index += 1) {
        const y = yOf(events[index]);
        if (y < box.y || y > box.y + box.height) continue;
        context.fillRect(baseX + profile + 10 * ratio + field * spreadOf(index), y, 1.6 * ratio, 1.6 * ratio);
      }
    }
    const labels = frame.screenLabels;
    if (labels) {
      text(labels[0], box.x + box.width - 6 * ratio, box.y + 11 * ratio, { color: PALETTE.muted, size: 10, align: 'right', ratio, weight: 700 });
      text(labels[1], box.x + box.width - 6 * ratio, box.y + box.height - 11 * ratio, { color: PALETTE.muted, size: 10, align: 'right', ratio, weight: 700 });
    } else {
      text(`z = ${number(to)}`, box.x + box.width - 6 * ratio, box.y + 11 * ratio, { color: PALETTE.quiet, size: 10, align: 'right', ratio });
      text(`z = ${number(from)}`, box.x + box.width - 6 * ratio, box.y + box.height - 11 * ratio, { color: PALETTE.quiet, size: 10, align: 'right', ratio });
    }
  }

  function drawCircuit(circuit, box, ratio) {
    const elements = circuit.elements;
    const slot = box.width / elements.length;
    const wireY = box.y + box.height * 0.42;
    const boxWidth = Math.min(slot * 0.62, 110 * ratio);
    const boxHeight = box.height * 0.56;
    let incoming = 1;
    elements.forEach((element, index) => {
      const centre = box.x + slot * (index + 0.5);
      const left = centre - boxWidth / 2;
      // Cable de entrada: su opacidad es la fracción del haz que llega hasta aquí.
      line([[box.x + slot * index, wireY], [left, wireY]], { color: `rgba(232, 181, 109, ${0.12 + 0.8 * incoming})`, width: (1 + 2.4 * incoming) * ratio });
      const progress = Math.min(1, Math.max(0, circuit.cursor - index));
      context.fillStyle = 'rgba(17, 18, 22, 0.95)';
      context.fillRect(left, wireY - boxHeight / 2, boxWidth, boxHeight);
      if (progress > 0) {
        context.fillStyle = 'rgba(232, 181, 109, 0.22)';
        context.fillRect(left, wireY - boxHeight / 2, boxWidth * progress, boxHeight);
      }
      context.strokeStyle = progress > 0 && progress < 1 ? PALETTE.gold : 'rgba(244, 239, 230, 0.34)';
      context.lineWidth = 1.3 * ratio;
      context.strokeRect(left, wireY - boxHeight / 2, boxWidth, boxHeight);
      text(element.label, centre, wireY - (element.detail ? 5 * ratio : 0), { color: PALETTE.ink, size: 12, align: 'center', ratio, weight: 700 });
      if (element.detail) text(element.detail, centre, wireY + 9 * ratio, { color: PALETTE.muted, size: 9.5, align: 'center', ratio });
      if (element.blocked > 0.001) {
        // Rama absorbida por el analizador.
        const drop = [[centre + boxWidth * 0.3, wireY + boxHeight / 2], [centre + boxWidth * 0.3, box.y + box.height - 2 * ratio]];
        line(drop, { color: `rgba(255, 157, 122, ${0.2 + 0.7 * element.blocked})`, width: (1 + 2 * element.blocked) * ratio });
        text(`${number(Math.round(element.blocked * 1000) / 10)} ${t.absorbed}`, centre + boxWidth * 0.3 + 5 * ratio, box.y + box.height - 6 * ratio, { color: PALETTE.reference, size: 9.5, ratio });
      }
      incoming = element.flow;
      line([[left + boxWidth, wireY], [box.x + slot * (index + 1), wireY]], { color: `rgba(232, 181, 109, ${0.12 + 0.8 * incoming})`, width: (1 + 2.4 * incoming) * ratio });
    });
  }

  const scaleState = { peak: 0 };

  function render(frame, { strict = false, duration = 0 } = {}) {
    const { width, height, ratio } = resize();
    context.setTransform(1, 0, 0, 1, 0, 0);
    context.fillStyle = '#03050b';
    context.fillRect(0, 0, width, height);
    const pad = 12 * ratio;
    const narrow = width < 560 * ratio;
    const sphereWidth = frame.bloch ? (narrow ? width * 0.5 : Math.min(width * 0.4, height * 1.02)) : 0;
    const right = { x: sphereWidth + pad * 0.5, y: pad, width: width - sphereWidth - pad * 1.5, height: height - 2 * pad };
    if (frame.bloch) drawSphere(frame.bloch, { x: 0, y: 0, width: sphereWidth, height }, ratio, strict);
    if (frame.circuit) {
      const strip = Math.min(74 * ratio, right.height * 0.26);
      drawCircuit(frame.circuit, { x: right.x, y: right.y, width: right.width, height: strip }, ratio);
      right.y += strip + pad * 0.5;
      right.height -= strip + pad * 0.5;
    }
    const hasMap = Boolean(frame.carpet && frame.carpet.rows.length);
    const screenWidth = hasMap ? Math.max(right.width * 0.3, 96 * ratio) : right.width;
    if (hasMap) drawMap(frame.carpet, { x: right.x, y: right.y, width: right.width - screenWidth - pad * 0.5, height: right.height }, ratio, duration);
    drawScreen(frame, { x: right.x + right.width - screenWidth, y: right.y, width: screenWidth, height: right.height }, ratio, scaleState);
  }

  return {
    overlay,
    kind: 'canvas2d',
    // El mapa de haces va dentro del escenario: el runner no tiene que repetirlo en un panel aparte.
    usesCarpet: true,
    resetScale() { scaleState.peak = 0; },
    render
  };
}

// Nombre con el que se publica el visor fuera de Fundamental Universe.
export const createBlochRenderer = createSpinRenderer;
