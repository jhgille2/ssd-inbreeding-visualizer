'use strict';
/* Single-Seed Descent visualizer.
 * Renders chromosome mosaics from docs/ssd_data.json (or a fresh in-browser
 * simulation of the same Haldane meiosis model) as dependency-free SVG. */

const NS = 'http://www.w3.org/2000/svg';
const COL_A = '#0072b2';   // Parent A
const COL_B = '#e69f00';   // Parent B
const COL_XO = '#b03a2e';  // crossover marks
const INK = '#1f2933', MUTED = '#5b6570';

const state = { data: null, genIdx: 0, showXO: true, showMeiosis: true };

/* ---------- SVG helpers ---------- */
function S(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function text(parent, x, y, str, attrs) {
  const t = S('text', Object.assign({ x, y }, attrs || {}), parent);
  t.textContent = str;
  return t;
}
function title(parent, str) {
  const t = document.createElementNS(NS, 'title');
  t.textContent = str;
  parent.appendChild(t);
}

/* ---------- data utilities ---------- */
// Compress a per-marker array into gapless [x0,x1) runs in cM.
function runsOf(arr, pos) {
  const runs = [];
  let s = 0;
  for (let i = 1; i <= arr.length; i++) {
    if (i === arr.length || arr[i] !== arr[s]) {
      runs.push({ x0: pos[s], x1: i < arr.length ? pos[i] : pos[pos.length - 1], v: arr[s] });
      s = i;
    }
  }
  return runs;
}
const originColor = v => (v === 0 ? COL_A : COL_B);
const originName = v => (v === 0 ? 'Parent A' : 'Parent B');

/* ---------- chromosome bar drawing ---------- */
function axis(g, X, y, plotW) {
  S('line', { x1: X(0), y1: y, x2: X(100), y2: y, stroke: INK, 'stroke-width': 1 }, g);
  for (let c = 0; c <= 100; c += 10) {
    S('line', { x1: X(c), y1: y, x2: X(c), y2: y + 5, stroke: INK, 'stroke-width': 1 }, g);
    text(g, X(c), y + 20, String(c), { 'text-anchor': 'middle', 'font-size': 11, fill: MUTED });
  }
  text(g, X(100) + 6, y + 4, 'cM', { 'font-size': 11, fill: MUTED });
}

function chromBar(g, X, y, h, runs, xos, opts) {
  opts = opts || {};
  runs.forEach((r, i) => {
    const rect = S('rect', {
      x: X(r.x0), y, width: Math.max(0.5, X(r.x1) - X(r.x0)), height: h,
      fill: originColor(r.v), stroke: 'none', rx: 2
    }, g);
    title(rect, `${originName(r.v)} — ${r.x0.toFixed(1)} to ${r.x1.toFixed(1)} cM`);
  });
  S('rect', { x: X(0), y, width: X(100) - X(0), height: h, fill: 'none', stroke: INK, 'stroke-width': 1.2, rx: 2 }, g);
  if (opts.showXO !== false) {
    (xos || []).forEach(xo => {
      const tick = S('line', {
        x1: X(xo), y1: y - 11, x2: X(xo), y2: y - 1,
        stroke: COL_XO, 'stroke-width': 2.4
      }, g);
      title(tick, `Crossover at ${xo.toFixed(1)} cM`);
      const xm = text(g, X(xo), y - 14, '×', {
        'text-anchor': 'middle', 'font-size': 13, fill: COL_XO, 'font-weight': 'bold'
      });
      title(xm, `Crossover at ${xo.toFixed(1)} cM`);
    });
  }
}

/* ---------- main figure ---------- */
function renderHomologs() {
  const d = state.data, pos = d.marker_pos_cm;
  const gen = d.generations[state.genIdx];
  const h1 = gen.homologs[0], h2 = gen.homologs[1];

  const W = 960, ML = 128, MR = 30, plotW = W - ML - MR;
  const X = cm => ML + (cm / 100) * plotW;
  const H = 268;

  const box = document.getElementById('homolog-viz');
  box.innerHTML = '';
  const svg = S('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img',
    'aria-label': `Chromosome pair of the ${gen.label} plant` }, box);

  const rows = [
    { label: 'Homolog 1', y: 30, runs: runsOf(h1.origin, pos), xos: h1.crossovers_cm },
    { label: 'Homolog 2', y: 96, runs: runsOf(h2.origin, pos), xos: h2.crossovers_cm },
  ];
  rows.forEach(r => {
    text(svg, ML - 12, r.y + 24, r.label, { 'text-anchor': 'end', 'font-size': 13, fill: INK });
    chromBar(svg, X, r.y, 36, r.runs, r.xos, { showXO: state.showXO });
  });

  // Genotype bar: code = h1*2 + h2 -> 0: A/A, 3: B/B, 1/2: heterozygous (phased).
  const gy = 178, gh = 36;
  text(svg, ML - 12, gy + 24, 'Genotype', { 'text-anchor': 'end', 'font-size': 13, fill: INK });
  const code = h1.origin.map((v, i) => v * 2 + h2.origin[i]);
  runsOf(code, pos).forEach(r => {
    const w = Math.max(0.5, X(r.x1) - X(r.x0));
    if (r.v === 0 || r.v === 3) {
      const rect = S('rect', { x: X(r.x0), y: gy, width: w, height: gh,
        fill: r.v === 0 ? COL_A : COL_B, rx: 2 }, svg);
      title(rect, `Homozygous ${r.v === 0 ? 'Parent A' : 'Parent B'} — ${r.x0.toFixed(1)} to ${r.x1.toFixed(1)} cM`);
    } else {
      // Heterozygous: code 1 = (h1: A, h2: B), code 2 = (h1: B, h2: A).
      const topC = originColor(r.v === 1 ? 0 : 1);
      const botC = originColor(r.v === 1 ? 1 : 0);
      const rt = S('rect', { x: X(r.x0), y: gy, width: w, height: gh / 2, fill: topC }, svg);
      const rb = S('rect', { x: X(r.x0), y: gy + gh / 2, width: w, height: gh / 2, fill: botC }, svg);
      title(rt, `Heterozygous — ${r.x0.toFixed(1)} to ${r.x1.toFixed(1)} cM (upper: Homolog 1)`);
      title(rb, `Heterozygous — ${r.x0.toFixed(1)} to ${r.x1.toFixed(1)} cM (lower: Homolog 2)`);
    }
  });
  S('rect', { x: X(0), y: gy, width: X(100) - X(0), height: gh, fill: 'none', stroke: INK, 'stroke-width': 1.2, rx: 2 }, svg);

  axis(svg, X, 236, plotW);
  document.getElementById('fig-title').textContent =
    `The ${gen.label} plant — one chromosome pair`;
}

/* ---------- stats + caption ---------- */
function renderStats() {
  const gen = state.data.generations[state.genIdx];
  const h1 = gen.homologs[0].origin, h2 = gen.homologs[1].origin;
  let hom = 0;
  for (let i = 0; i < h1.length; i++) if (h1[i] === h2[i]) hom++;
  const obsHom = hom / h1.length;
  const expHom = 1 - gen.expected_heterozygosity;
  const nXO = gen.homologs[0].crossovers_cm.length + gen.homologs[1].crossovers_cm.length;
  const pct = v => (v * 100).toFixed(1) + '%';

  const stats = [
    ['Observed homozygosity', pct(obsHom), `${hom} of ${h1.length} markers`],
    ['Expected homozygosity', pct(expHom), `1 − (1/2)^${gen.selfing_generations}`],
    ['Heterozygous markers', String(h1.length - hom), 'still segregating'],
    ['Crossovers in formative meioses', String(nXO), 'visible breakpoints above'],
  ];
  const grid = document.getElementById('stats-grid');
  grid.innerHTML = '';
  stats.forEach(([k, v, s]) => {
    const div = document.createElement('div');
    div.className = 'stat';
    div.innerHTML = `<div class="k">${k}</div><div class="v">${v}</div><div class="s">${s}</div>`;
    grid.appendChild(div);
  });

  const cap = document.getElementById('gen-caption');
  const L = gen.label;
  if (L === 'F1') {
    cap.innerHTML = `<strong>F1.</strong> The F1 carries one intact chromosome from each parent: ` +
      `every marker is heterozygous and no recombination has occurred yet. Recombination first ` +
      `appears in the gametes this plant produces.`;
  } else {
    const prev = state.data.generations[state.genIdx - 1].label;
    cap.innerHTML = `<strong>${L}.</strong> One seed from the ${prev} plant was grown and selfed. ` +
      `This plant received one gamete from each of two independent meioses in its parent; ` +
      `${nXO} crossover${nXO === 1 ? '' : 's'} left visible breakpoints. It is ${pct(obsHom)} ` +
      `homozygous (expectation ${pct(expHom)}): each selfing fixes, on average, half of the ` +
      `heterozygosity that remains.`;
  }
}

/* ---------- animated meiosis ---------- */
const MEIOSIS = { hh: 0, runId: 0, ui: null };
const lerp = (a, b, k) => a + (b - a) * k;
const clamp01 = v => Math.min(1, Math.max(0, v));
const easeInOut = k => (k < 0.5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2);
const hasRAF = typeof requestAnimationFrame === 'function';
const reducedMotion = () =>
  typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;

const M_W = 960, M_ML = 170, M_MR = 30, M_PLOT = M_W - M_ML - M_MR, M_H = 372, M_RH = 24;
const Y_SEP = [86, 120, 196, 230];    // homologs apart
const Y_BUN = [122, 148, 178, 204];  // paired bundle
const Y_SPR = [58, 132, 206, 280];   // segregated gametes
const X_MID = 175;                   // crossover marks between bundled rows 1 and 2
const Xof = cm => M_ML + (cm / 100) * M_PLOT;

// One meiosis in the parent plant, as four chromatids. Rows 0/3 are the parental
// (non-interacting) sister chromatids; rows 1/2 are the non-sister pair that
// crosses over. rows[1] is the recombinant starting on homolog 1, rows[2] the
// reciprocal recombinant; transmittedRow is the one kept by single-seed descent.
function buildMeiosisModel(gi, hh) {
  const d = state.data, pos = d.marker_pos_cm;
  const parent = d.generations[gi - 1], gen = d.generations[gi];
  const P1 = parent.homologs[0].origin, P2 = parent.homologs[1].origin;
  const T = gen.homologs[hh].origin;
  const xo = gen.homologs[hh].crossovers_cm.slice().sort((a, b) => a - b);
  const rec = s => {
    const g = new Array(pos.length);
    for (let i = 0; i < pos.length; i++) {
      let k = 0;
      while (k < xo.length && xo[k] <= pos[i]) k++;
      g[i] = ((k + s) % 2 === 0) ? P1[i] : P2[i];
    }
    return g;
  };
  const r0 = rec(0), r1 = rec(1);
  const same = (a, b) => { for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false; return true; };
  const start = same(r0, T) ? 0 : 1;
  return {
    P1, P2, pos, xo, start, parentLabel: parent.label, genLabel: gen.label, hh,
    rows: [P1, r0, r1, P2], transmittedRow: 1 + start,
  };
}

// Intervals (in cM) where a chromatid's origin changes after the exchange.
function exchangedIntervals(before, after, pos) {
  const ivs = [];
  let s = -1;
  for (let i = 0; i < before.length; i++) {
    const ch = before[i] !== after[i];
    if (ch && s < 0) s = i;
    if (!ch && s >= 0) { ivs.push({ x0: pos[s], x1: pos[i], from: before[s], to: after[s] }); s = -1; }
  }
  if (s >= 0) ivs.push({ x0: pos[s], x1: pos[pos.length - 1], from: before[s], to: after[s] });
  return ivs;
}

function staticBar(g, runs, h) {
  runs.forEach(r => {
    S('rect', {
      x: Xof(r.x0), y: 0, width: Math.max(0.5, Xof(r.x1) - Xof(r.x0)),
      height: h, fill: originColor(r.v), rx: 2,
    }, g);
  });
  S('rect', {
    x: Xof(0), y: 0, width: Xof(100) - Xof(0), height: h,
    fill: 'none', stroke: INK, 'stroke-width': 1, rx: 2,
  }, g);
}

function renderMeiosis() {
  const sec = document.getElementById('meiosis-section');
  sec.style.display = state.showMeiosis ? '' : 'none';
  if (!state.showMeiosis) return;
  MEIOSIS.runId++; // cancel any running animation
  const box = document.getElementById('meiosis-viz');
  box.innerHTML = '';
  const gi = state.genIdx;

  if (gi === 0) {
    const card = document.createElement('div');
    card.className = 'meiosis-card';
    card.innerHTML = '<h3>No meiosis to replay</h3><p>The F1\u2019s two homologs are the intact ' +
      'parental chromosomes \u2014 one full copy from each parent. Crossovers first appear in the ' +
      'gametes produced <em>by</em> this plant (advance to F2 to watch one happen).</p>';
    box.appendChild(card);
    return;
  }

  const ctrl = document.createElement('div');
  ctrl.className = 'meiosis-controls';
  const seg = document.createElement('div');
  seg.className = 'seg';
  seg.setAttribute('role', 'group');
  seg.setAttribute('aria-label', 'Choose meiosis');
  [0, 1].forEach(h => {
    const b = document.createElement('button');
    b.textContent = `Meiosis ${h + 1} \u2192 Homolog ${h + 1}`;
    if (h === MEIOSIS.hh) b.className = 'active';
    b.addEventListener('click', () => { MEIOSIS.hh = h; renderMeiosis(); });
    seg.appendChild(b);
  });
  const replay = document.createElement('button');
  replay.id = 'replay-meiosis';
  replay.textContent = '\u21BB Replay meiosis';
  replay.addEventListener('click', () => renderMeiosis());
  ctrl.appendChild(seg);
  ctrl.appendChild(replay);
  box.appendChild(ctrl);

  const stage = document.createElement('div');
  stage.className = 'meiosis-stage';
  box.appendChild(stage);
  const svg = S('svg', {
    viewBox: `0 0 ${M_W} ${M_H}`, role: 'img', 'aria-label': 'Animated meiosis',
  }, stage);

  const model = buildMeiosisModel(gi, MEIOSIS.hh);
  const rows = [];
  for (let i = 0; i < 4; i++) {
    const g = S('g', { transform: `translate(0 ${Y_SEP[i]})` }, svg);
    staticBar(g, runsOf(model.rows[i], model.pos), M_RH);
    rows.push(g);
  }
  const homLabels = [
    text(svg, M_ML - 12, (Y_SEP[0] + Y_SEP[1]) / 2 + 8, 'Homolog 1',
      { 'text-anchor': 'end', 'font-size': 13, fill: INK }),
    text(svg, M_ML - 12, (Y_SEP[2] + Y_SEP[3]) / 2 + 8, 'Homolog 2',
      { 'text-anchor': 'end', 'font-size': 13, fill: INK }),
  ];
  text(svg, M_ML, 30, `One meiosis in the ${model.parentLabel} plant`,
    { 'font-size': 13, fill: MUTED });
  const top = S('g', {}, svg);
  const cap = document.createElement('p');
  cap.className = 'phase-caption';
  stage.appendChild(cap);

  MEIOSIS.ui = { svg, rows, homLabels, top, cap, model, xg: null, xmarks: [], movers: null, moverRects: [] };
  playMeiosis();
}

function playMeiosis() {
  const ui = MEIOSIS.ui, model = ui.model;
  const run = ++MEIOSIS.runId;
  const t = model.transmittedRow;
  const n = model.xo.length;
  const ivs = exchangedIntervals(model.P1, model.rows[1], model.pos);
  const setRows = ys => ys.forEach((y, i) =>
    ui.rows[i].setAttribute('transform', `translate(0 ${y})`));
  const setCap = html => { ui.cap.innerHTML = html; };

  const steps = [
    {
      dur: 900,
      cap: '<strong>1 \u00B7 Synapsis.</strong> The parent\u2019s homologs pair up. Each homolog ' +
        'is two identical sister chromatids; crossovers will form between non-sister chromatids.',
      frame: k => {
        setRows(Y_SEP.map((y, i) => lerp(y, Y_BUN[i], k)));
        ui.homLabels.forEach(el => el.setAttribute('opacity', 1 - k));
      },
    },
    {
      dur: 1100,
      cap: n
        ? `<strong>2 \u00B7 Crossing over.</strong> ${n} crossover${n === 1 ? '' : 's'} ` +
          `form${n === 1 ? 's' : ''} between non-sister chromatids (\u00D7 marks).`
        : '<strong>2 \u00B7 Crossing over.</strong> No crossovers in this meiosis ' +
          '\u2014 all four chromatids stay parental.',
      start: () => {
        ui.xg = S('g', { opacity: 0 }, ui.svg);
        ui.xmarks = model.xo.map(xo => {
          const g = S('g', { transform: `translate(${Xof(xo)} ${X_MID}) scale(0.4)` }, ui.xg);
          S('line', { x1: -7, y1: -7, x2: 7, y2: 7, stroke: COL_XO, 'stroke-width': 2.6 }, g);
          S('line', { x1: -7, y1: 7, x2: 7, y2: -7, stroke: COL_XO, 'stroke-width': 2.6 }, g);
          const tg = S('g', {}, g);
          title(tg, `Crossover at ${xo.toFixed(1)} cM`);
          return { g, x: Xof(xo) };
        });
      },
      frame: k => {
        if (!ui.xg) return;
        const s = k < 0.7 ? (k / 0.7) * 1.2 : 1.2 - 0.2 * ((k - 0.7) / 0.3);
        ui.xg.setAttribute('opacity', k);
        ui.xmarks.forEach(m =>
          m.g.setAttribute('transform', `translate(${m.x} ${X_MID}) scale(${s.toFixed(3)})`));
      },
    },
    {
      dur: 2400,
      cap: ivs.length
        ? '<strong>3 \u00B7 Exchange.</strong> At each crossover the non-sister chromatids ' +
          'break and rejoin \u2014 segments are physically swapped between them.'
        : '<strong>3 \u00B7 Exchange.</strong> With no crossovers there is nothing to exchange.',
      start: () => {
        ui.movers = S('g', {}, ui.top);
        ui.moverRects = [];
        ivs.forEach(iv => {
          const w = Math.max(0.5, Xof(iv.x1) - Xof(iv.x0));
          // Each moving rect starts exactly over the same-colored segment it
          // detaches from, then slides to the other chromatid.
          const rA = S('rect', {
            x: Xof(iv.x0), y: Y_BUN[2], width: w, height: M_RH,
            fill: originColor(iv.to), stroke: INK, 'stroke-width': 1, rx: 2,
          }, ui.movers);
          const rB = S('rect', {
            x: Xof(iv.x0), y: Y_BUN[1], width: w, height: M_RH,
            fill: originColor(iv.from), stroke: INK, 'stroke-width': 1, rx: 2,
          }, ui.movers);
          ui.moverRects.push({ rA, rB });
        });
      },
      frame: k => {
        ui.moverRects.forEach((m, idx) => {
          const kk = easeInOut(clamp01(k * 1.6 - idx * 0.25));
          m.rA.setAttribute('y', lerp(Y_BUN[2], Y_BUN[1], kk).toFixed(1));
          m.rB.setAttribute('y', lerp(Y_BUN[1], Y_BUN[2], kk).toFixed(1));
        });
      },
      end: () => {
        if (ui.movers) ui.movers.setAttribute('opacity', 0);
        [1, 2].forEach(i => {
          ui.rows[i].innerHTML = '';
          staticBar(ui.rows[i], runsOf(model.rows[i], model.pos), M_RH);
        });
      },
    },
    {
      dur: 1400,
      cap: `<strong>4 \u00B7 Segregation.</strong> The four chromatids separate into gametes. ` +
        `Single-seed descent keeps a single seed \u2014 only the highlighted gamete is transmitted, ` +
        `becoming Homolog ${model.hh + 1} of the ${model.genLabel} plant.`,
      frame: k => {
        setRows(Y_BUN.map((y, i) => lerp(y, Y_SPR[i], k)));
        if (ui.xg) ui.xg.setAttribute('opacity', 1 - k);
      },
      end: () => {
        const kinds = n
          ? ['parental', 'recombinant', 'recombinant', 'parental']
          : ['parental', 'parental', 'parental', 'parental'];
        kinds.forEach((kind, i) => {
          text(ui.svg, M_ML - 12, Y_SPR[i] + 17,
            i === t ? `${kind} \u2605 transmitted` : kind,
            {
              'text-anchor': 'end', 'font-size': 12,
              fill: i === t ? INK : MUTED, 'font-weight': i === t ? 'bold' : 'normal',
            });
          if (i !== t) ui.rows[i].setAttribute('opacity', 0.35);
        });
        S('rect', {
          x: Xof(0) - 7, y: Y_SPR[t] - 7, width: M_PLOT + 14, height: M_RH + 14,
          fill: 'none', stroke: '#c9a227', 'stroke-width': 3, rx: 8,
        }, ui.top);
      },
    },
  ];

  if (!hasRAF || reducedMotion()) {
    // No animation support (or reduced motion): jump straight to the final state.
    steps.forEach(st => { setCap(st.cap); if (st.start) st.start(); st.frame(1); if (st.end) st.end(); });
    return;
  }
  let i = 0;
  const next = () => {
    if (run !== MEIOSIS.runId) return;
    if (i >= steps.length) return;
    const st = steps[i++];
    setCap(st.cap);
    if (st.start) st.start();
    const t0 = performance.now();
    const tick = now => {
      if (run !== MEIOSIS.runId) return;
      const k = Math.min(1, (now - t0) / st.dur);
      st.frame(easeInOut(k));
      if (k < 1) requestAnimationFrame(tick);
      else { if (st.end) st.end(); setTimeout(next, 280); }
    };
    requestAnimationFrame(tick);
  };
  next();
}

/* ---------- trajectory ---------- */
function renderTrajectory() {
  const d = state.data;
  const W = 960, H = 320, ML = 64, MR = 30, MT = 24, MB = 44;
  const plotW = W - ML - MR, plotH = H - MT - MB;
  const X = i => ML + (i / (d.generations.length - 1)) * plotW;
  const Y = v => MT + (1 - v) * plotH;

  const box = document.getElementById('trajectory-viz');
  box.innerHTML = '';
  const svg = S('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img',
    'aria-label': 'Homozygosity across generations' }, box);

  [0, 0.25, 0.5, 0.75, 1].forEach(v => {
    S('line', { x1: ML, y1: Y(v), x2: W - MR, y2: Y(v), stroke: '#d9dee3' }, svg);
    text(svg, ML - 10, Y(v) + 4, Math.round(v * 100) + '%',
      { 'text-anchor': 'end', 'font-size': 11, fill: MUTED });
  });
  d.generations.forEach((g, i) => {
    text(svg, X(i), H - 16, g.label, { 'text-anchor': 'middle', 'font-size': 12, fill: INK });
  });

  const expPts = d.generations.map((g, i) => `${X(i)},${Y(1 - g.expected_heterozygosity)}`).join(' ');
  S('polyline', { points: expPts, fill: 'none', stroke: MUTED, 'stroke-width': 2, 'stroke-dasharray': '7 5' }, svg);
  const obs = d.generations.map(g => {
    const h1 = g.homologs[0].origin, h2 = g.homologs[1].origin;
    let n = 0;
    for (let i = 0; i < h1.length; i++) if (h1[i] === h2[i]) n++;
    return n / h1.length;
  });
  const obsPts = obs.map((v, i) => `${X(i)},${Y(v)}`).join(' ');
  S('polyline', { points: obsPts, fill: 'none', stroke: COL_A, 'stroke-width': 2.5 }, svg);
  obs.forEach((v, i) => {
    const c = S('circle', { cx: X(i), cy: Y(v), r: i === state.genIdx ? 8 : 5,
      fill: COL_A, stroke: '#fff', 'stroke-width': 2,
      style: i === state.genIdx ? '' : 'cursor:pointer' }, svg);
    if (i !== state.genIdx) {
      c.addEventListener('click', () => { state.genIdx = i; syncControls(); render(); });
    }
    title(c, `${d.generations[i].label}: observed ${(v * 100).toFixed(1)}%, expected ${((1 - d.generations[i].expected_heterozygosity) * 100).toFixed(1)}%`);
  });
  text(svg, W - MR, MT - 8, '— expected', { 'text-anchor': 'end', 'font-size': 12, fill: MUTED });
  text(svg, W - MR, MT + 10, '— observed', { 'text-anchor': 'end', 'font-size': 12, fill: COL_A, 'font-weight': 'bold' });
}

/* ---------- in-browser simulation (same model, fresh seed) ---------- */
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function rpois(rng, lam) {
  const L = Math.exp(-lam);
  let k = 0, p = 1;
  do { p *= rng(); k++; } while (p > L);
  return k - 1;
}
function meiosisJS(h1, h2, pos, rng) {
  const n = rpois(rng, (pos[pos.length - 1] - pos[0]) / 100);
  const xo = [];
  for (let i = 0; i < n; i++) xo.push(rng() * 100);
  xo.sort((a, b) => a - b);
  const start = rng() < 0.5 ? 0 : 1;
  const g = new Array(pos.length);
  for (let i = 0; i < pos.length; i++) {
    let s = 0;
    while (s < xo.length && xo[s] <= pos[i]) s++;
    g[i] = ((s + start) % 2 === 0) ? h1[i] : h2[i];
  }
  return { g, xo: xo.map(x => Math.round(x * 1000) / 1000) };
}
function simulateLine(seed) {
  const rng = mulberry32(seed);
  const N = 1000;
  const pos = [];
  for (let i = 0; i < N; i++) pos.push(Math.round((i * 100 / (N - 1)) * 10000) / 10000);
  const mk = (label, sg, m1, m2) => ({
    label, selfing_generations: sg,
    expected_heterozygosity: Math.pow(0.5, sg),
    homologs: [
      { origin: m1.g, crossovers_cm: m1.xo },
      { origin: m2.g, crossovers_cm: m2.xo },
    ],
  });
  const zeros = new Array(N).fill(0), ones = new Array(N).fill(1);
  const generations = [mk('F1', 0, { g: zeros, xo: [] }, { g: ones, xo: [] })];
  let plant = [zeros, ones];
  for (let i = 1; i <= 7; i++) {
    const g1 = meiosisJS(plant[0], plant[1], pos, rng);
    const g2 = meiosisJS(plant[0], plant[1], pos, rng);
    plant = [g1.g, g2.g];
    generations.push(mk('F' + (i + 1), i, g1, g2));
  }
  return {
    meta: {
      title: 'Single-seed descent: from cross to inbred line',
      organism: 'Glycine max (soybean); one representative chromosome is shown (soybean has 20 chromosome pairs, 2n = 40)',
      chromosome_length_cm: 100.0,
      n_markers: N,
      parents: { 0: 'Parent A', 1: 'Parent B' },
      model: 'Single-seed descent simulated in-browser: F1 selfed; one random progeny selfed each generation to F8. Meiosis: Haldane model (Poisson crossover count, uniform positions, no interference); one random chromatid transmitted per gamete.',
      generator: 'in-browser simulation (app.js)',
      seed,
    },
    marker_pos_cm: pos,
    generations,
  };
}

/* ---------- controls + boot ---------- */
function syncControls() {
  document.getElementById('gen-slider').value = state.genIdx;
  document.getElementById('gen-label').textContent = state.data.generations[state.genIdx].label;
}
function render() {
  renderHomologs();
  renderStats();
  renderMeiosis();
  renderTrajectory();
}
function setSeedLabel() {
  const el = document.getElementById('line-seed');
  el.textContent = state.data.meta.generator.startsWith('in-browser')
    ? `Line seed: ${state.data.meta.seed} (simulated in your browser)`
    : `Line seed: ${state.data.meta.seed} (default dataset)`;
}

async function init() {
  try {
    const res = await fetch('ssd_data.json');
    if (!res.ok) throw new Error('fetch failed');
    state.data = await res.json();
  } catch (e) {
    document.getElementById('homolog-viz').innerHTML =
      '<p style="padding:1rem">Could not load <code>ssd_data.json</code>. ' +
      'Serve the <code>docs/</code> directory over HTTP (e.g. <code>python3 -m http.server</code>) ' +
      'or open the hosted site.</p>';
    return;
  }
  setSeedLabel();

  const slider = document.getElementById('gen-slider');
  slider.max = state.data.generations.length - 1;
  slider.addEventListener('input', () => { state.genIdx = +slider.value; syncControls(); render(); });
  document.getElementById('prev-gen').addEventListener('click', () => {
    state.genIdx = Math.max(0, state.genIdx - 1); syncControls(); render();
  });
  document.getElementById('next-gen').addEventListener('click', () => {
    state.genIdx = Math.min(state.data.generations.length - 1, state.genIdx + 1);
    syncControls(); render();
  });
  document.getElementById('toggle-xo').addEventListener('change', e => {
    state.showXO = e.target.checked; render();
  });
  document.getElementById('toggle-meiosis').addEventListener('change', e => {
    state.showMeiosis = e.target.checked; render();
  });
  document.getElementById('new-line-btn').addEventListener('click', () => {
    const seed = Math.floor(Math.random() * 1e9);
    state.data = simulateLine(seed);
    state.genIdx = 0;
    setSeedLabel();
    syncControls(); render();
  });

  syncControls();
  render();
}
document.addEventListener('DOMContentLoaded', init);
