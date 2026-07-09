// DGGS aspect-ratio explorer — reads the static data built by build_data.py
// (out/manifest.json, out/histograms.json, out/globe/{sys}_r{res}_*.{f32,u32,json}).
// Two views: dynamic overlaid histograms (Observable Plot) and two synced
// orthographic globes (ajglobe), cells colored by aspect ratio.
import { Orb, lnglatToQuat, quatToLngLat } from './vendor/ajglobe.min.js';
import { DNC_GREY, viridis, lut, fetchBin, sortedFinite, quantileT } from './_shared.js';

const fmt = (x, d = 4) => (x == null ? '—' : x.toFixed(d));
const $ = (sel) => document.querySelector(sel);

let M;        // manifest
let HIST;     // histograms.json {edges, nbins, amax, data}
const state = {
  series: new Map(),   // key "sys|res" -> {sys, res} (insertion order = draw order)
  bins: 0,             // set from the slider's HTML value in initHistControls
  yMode: 'density',
  yScale: 'linear',
};

// ---------- tabs ----------
function initTabs() {
  const tabs = [['#tab-hist', '#view-hist'], ['#tab-globe', '#view-globe']];
  for (const [btn, view] of tabs) {
    $(btn).addEventListener('click', () => {
      for (const [b, v] of tabs) {
        $(b).classList.toggle('active', b === btn);
        $(v).classList.toggle('active', v === view);
      }
      if (view === '#view-globe') globe.kick();
    });
  }
}

// ================= Histograms =================
function buildSeriesPicker() {
  const root = $('#seriesPicker');
  root.innerHTML = '';
  for (const sys of M.systems) {
    const block = document.createElement('div');
    block.className = 'sys-block';
    const name = document.createElement('div');
    name.className = 'sys-name';
    name.innerHTML = `<span class="swatch" style="background:${M.colors[sys]}"></span>${M.labels[sys]}`;
    block.appendChild(name);
    const chips = document.createElement('div');
    chips.className = 'chips';
    for (const res of M.hist_res[sys]) {
      const key = `${sys}|${res}`;
      const chip = document.createElement('span');
      chip.className = 'chip' + (res === M.target_res[sys] ? ' target' : '');
      chip.textContent = res;
      chip.style.color = M.colors[sys];
      chip.dataset.key = key;
      chip.addEventListener('click', () => toggleSeries(sys, res, chip));
      chips.appendChild(chip);
    }
    block.appendChild(chips);
    root.appendChild(block);
  }
}

function toggleSeries(sys, res, chip) {
  const key = `${sys}|${res}`;
  if (state.series.has(key)) {
    state.series.delete(key);
    chip.classList.remove('on');
    chip.style.background = '';
    chip.style.color = M.colors[sys];
  } else {
    state.series.set(key, { sys, res });
    chip.classList.add('on');
    chip.style.background = M.colors[sys];
    chip.style.color = '#fff';
  }
  renderHist();
}

// Color a series: base = system color; if several resolutions of one system are
// selected, lighten by rank so they read apart while staying "that system".
function seriesColors() {
  const bySys = {};
  for (const { sys } of state.series.values()) (bySys[sys] ??= 0, bySys[sys]++);
  const seen = {};
  const out = new Map();
  for (const [key, { sys }] of state.series) {
    const n = bySys[sys];
    const rank = (seen[sys] = (seen[sys] ?? -1) + 1);
    const base = M.colors[sys];
    const c = n > 1 ? d3.interpolateRgb(base, '#dfe7ef')(0.62 * rank / (n - 1)) : base;
    out.set(key, c);
  }
  return out;
}

// Re-aggregate the fixed 256-bin counts into ~`state.bins` display bins by
// merging adjacent fine bins (keeps edges aligned to the stored grid). Returns
// step-outline points for Plot.line.
function seriesLineData(key, color) {
  const { sys, res } = state.series.get(key);
  const rec = HIST.data[sys]?.[res];
  if (!rec) return [];
  const edges = HIST.edges, nfine = HIST.nbins;
  const factor = Math.max(1, Math.round(nfine / state.bins));
  const pts = [];
  for (let i = 0; i < nfine; i += factor) {
    const j = Math.min(i + factor, nfine);
    let count = 0;
    for (let k = i; k < j; k++) count += rec.counts[k];
    const x0 = edges[i], x1 = edges[j];
    let y = count;
    if (state.yMode === 'density') y = count / (rec.n * (x1 - x0));
    if (state.yScale === 'log' && y <= 0) continue;   // log breaks the line at gaps
    // two points per bin -> flat top + vertical step under a linear curve
    pts.push({ key, color, label: seriesLabel(sys, res), x: x0, y });
    pts.push({ key, color, label: seriesLabel(sys, res), x: x1, y });
  }
  return pts;
}

function seriesLabel(sys, res) {
  return `${M.labels[sys]} ${M.res_prefix[sys]}${res}`;
}

function renderHist() {
  const colors = seriesColors();
  const keys = [...state.series.keys()];
  const all = keys.flatMap((k) => seriesLineData(k, colors.get(k)));
  const plotEl = $('#histPlot');

  if (!all.length) {
    plotEl.innerHTML = '<p class="hint" style="padding:40px 0;text-align:center">'
      + 'Pick one or more system · resolution series on the left.</p>';
    $('#statsTable').innerHTML = '';
    return;
  }

  // Dynamic x-domain from the selected series (zoom to their span).
  const xmax = d3.max(keys, (k) => {
    const { sys, res } = state.series.get(k);
    return Math.min(HIST.data[sys][res].max, HIST.amax);
  });

  const plot = Plot.plot({
    width: plotEl.clientWidth || 760,
    height: 460,
    marginLeft: 58, marginBottom: 42, marginRight: 16, marginTop: 14,
    style: { background: 'transparent', color: '#9aa7b4', fontSize: '11px' },
    x: { domain: [1, xmax], label: 'aspect ratio →', grid: true,
         labelAnchor: 'right', tickFormat: '~f' },
    y: {
      type: state.yScale,
      label: state.yMode === 'density' ? '↑ density' : '↑ count',
      grid: true,
      clamp: true,
    },
    color: { type: 'identity' },   // each datum carries its literal hex in `color`
    marks: [
      // Baseline only in linear scale; log has no zero and a tiny floor would
      // stretch the whole axis down to it.
      ...(state.yScale === 'linear' ? [Plot.ruleY([0], { stroke: '#2a323d' })] : []),
      Plot.line(all, {
        x: 'x', y: 'y', z: 'key', stroke: 'color',
        strokeWidth: 1.6, curve: 'linear',
      }),
      Plot.tip(all, Plot.pointerX({
        x: 'x', y: 'y', stroke: 'color',
        title: (d) => `${d.label}\nAR ${fmt(d.x, 3)}\n${state.yMode} ${d.y.toExponential(2)}`,
        fontFamily: 'monospace',
      })),
    ],
  });
  plotEl.innerHTML = '';
  plotEl.appendChild(plot);
  renderStats(colors);
}

function renderStats(colors) {
  const rows = [...state.series].map(([key, { sys, res }]) => {
    const r = HIST.data[sys][res];
    const dnc = r.dnc ? `<td class="dnc">${r.dnc.toLocaleString()}</td>` : '<td>0</td>';
    return `<tr>
      <td><span class="sw" style="background:${colors.get(key)}"></span>${seriesLabel(sys, res)}</td>
      <td>${r.n.toLocaleString()}</td>${dnc}
      <td>${fmt(r.min)}</td><td>${fmt(r.median)}</td><td>${fmt(r.p99)}</td><td>${fmt(r.max)}</td>
    </tr>`;
  });
  $('#statsTable').innerHTML =
    `<tr><th>series</th><th>n</th><th>DNC</th><th>min</th><th>median</th><th>p99</th><th>max</th></tr>`
    + rows.join('');
}

function initHistControls() {
  const slider = $('#binSlider'), binVal = $('#binVal');
  const sync = () => { binVal.textContent = state.bins; };
  // Coalesce drag events to one re-render per frame — each render rebuilds the
  // whole plot SVG.
  let raf = 0;
  slider.addEventListener('input', () => {
    state.bins = +slider.value; sync();
    if (!raf) raf = requestAnimationFrame(() => { raf = 0; renderHist(); });
  });
  state.bins = +slider.value; sync();

  for (const [grp, prop] of [['#yMode', 'yMode'], ['#yScale', 'yScale']]) {
    $(grp).addEventListener('click', (e) => {
      const b = e.target.closest('button'); if (!b) return;
      state[prop] = b.dataset.v;
      $(grp).querySelectorAll('button').forEach((x) => x.classList.toggle('on', x === b));
      renderHist();
    });
  }
  $('#histNote').textContent =
    `Stored as ${HIST.nbins} fixed bins over 1–${fmt(HIST.amax, 2)}; `
    + `anything beyond folds into overflow (working resolutions never reach it).`;
}

// ================= Globe (ajglobe Orb, WebGL2) =================
const globe = (() => {
  const gamma = 0.4;
  const css = (c) => `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;

  const cache = new Map();           // key -> Promise<[pos, starts, ar, ids]>
  let domainMode = 'adaptive';       // default: each globe scaled to its own AR range
  const panels = {};                 // id -> {orb, layer, label, ar, ids, qt, arMax, epoch}
  let spinning = false, spinReq = null, lastT = 0;

  // Two coloring modes:
  //   shared    — absolute AR over [1, global max] with a sub-linear (gamma)
  //               stretch, so the two globes are directly comparable.
  //   per-globe — color by the cell's quantile within its own globe, spreading
  //               colors evenly however the distribution is shaped (maximizes
  //               visible structure; not cross-comparable).
  function fillFn(p) {
    if (domainMode === 'shared') {
      const hi = M.globe_ar_max;
      return (i) => {
        const a = p.ar[i];
        if (!Number.isFinite(a)) return DNC_GREY;
        return lut(Math.pow((Math.min(a, hi) - 1) / (hi - 1), gamma));
      };
    }
    return (i) => {
      const t = p.qt[i];
      return Number.isNaN(t) ? DNC_GREY : lut(t);
    };
  }

  function buildCard(id, defaultSys, defaultRes) {
    const card = $(`#card-${id}`);
    const sysSel = document.createElement('select');
    const resSel = document.createElement('select');
    for (const s of M.systems) sysSel.add(new Option(M.labels[s], s));
    sysSel.value = defaultSys;
    const head = document.createElement('div');
    head.className = 'globe-head';
    head.append(sysSel, resSel);
    const label = document.createElement('span');
    label.className = 'label';
    head.appendChild(label);

    const holder = document.createElement('div');
    holder.className = 'globe-holder';
    const canvas = document.createElement('canvas');
    holder.appendChild(canvas);

    card.innerHTML = '';
    card.append(head, holder);

    const orb = new Orb(canvas, { background: '#0b0e13', sphere: '#11151c' });
    orb.lookAt(0, 20);
    orb.borders({ color: '#0b0e13', width: 1 });   // Natural Earth outlines (fetched from a CDN at runtime)
    panels[id] = { orb, label, layer: null, epoch: 0 };

    // Hover: GPU pick -> highlight tint + the shared tooltip, tracked to the cursor.
    const tip = $('#tooltip');
    orb.on('hover', (e) => {
      const p = panels[id];
      orb.highlight(e.index ?? -1, p.layer);
      if (e.index == null || !p.ids) { tip.style.opacity = 0; return; }
      const a = p.ar[e.index];
      tip.innerHTML = `${p.ids[e.index]}<br>AR ${Number.isFinite(a) ? fmt(a, 4) : 'DNC'}`;
      const r = canvas.getBoundingClientRect();
      tip.style.left = `${r.left + e.x + 14}px`;
      tip.style.top = `${r.top + e.y + 14}px`;
      tip.style.opacity = 1;
    });
    canvas.addEventListener('pointerleave', () => { tip.style.opacity = 0; });

    const fillRes = () => {
      resSel.innerHTML = '';
      for (const r of M.globe_res[sysSel.value]) {
        resSel.add(new Option(`${M.res_prefix[sysSel.value]}${r}`, r));
      }
    };
    const reload = () => loadPanel(id, sysSel.value, +resSel.value);
    sysSel.addEventListener('change', () => { fillRes(); reload(); });
    resSel.addEventListener('change', reload);
    fillRes();
    resSel.value = defaultRes;

    loadPanel(id, defaultSys, defaultRes);
  }

  async function loadPanel(id, sys, res) {
    const key = `${sys}_r${res}`;
    if (!cache.has(key)) {
      cache.set(key, Promise.all([
        fetchBin(`out/globe/${key}_pos.f32`, Float32Array),
        fetchBin(`out/globe/${key}_idx.u32`, Uint32Array),
        fetchBin(`out/globe/${key}_ar.f32`, Float32Array),
        fetch(`out/globe/${key}_ids.json`).then((r) => r.json()),
      ]));
    }
    const p = panels[id];
    const e = ++p.epoch;                 // ignore a stale load after a quick re-select
    const [pos, starts, ar, ids] = await cache.get(key);
    if (e !== p.epoch) return;
    const sorted = sortedFinite(ar);
    p.ar = ar;
    p.ids = ids;
    p.qt = quantileT(ar, sorted);
    p.arMax = sorted.length ? sorted[sorted.length - 1] : M.globe_ar_max;
    p.label.textContent = `${ids.length.toLocaleString()} cells · max AR ${fmt(p.arMax, 3)}`;
    if (p.layer) p.layer.remove();
    p.layer = p.orb.polygons({ lnglat: pos, starts, fill: fillFn(p) });
  }

  // Auto-spin drives globe A; the viewchange sync carries globe B along.
  function spinLoop(t) {
    if (!spinning) return;
    if (lastT) {
      const orb = panels.A.orb;
      const { q, zoom } = orb.getView();
      const v = quatToLngLat(q);
      orb.setView({ q: lnglatToQuat(v.lng + (t - lastT) * 0.01, v.lat, v.roll), zoom });
    }
    lastT = t;
    spinReq = requestAnimationFrame(spinLoop);
  }

  function updateLegend() {
    // Shared: bar color follows the gamma stretch, so a swatch's position reads
    // as an absolute AR (labelled 1.0 .. global max). Per-globe: color is by
    // quantile, so the bar is a plain low→high key and each globe's own max
    // lives in its card label.
    const n = 48;
    const stops = Array.from({ length: n }, (_, i) => {
      const f = i / (n - 1);
      const t = domainMode === 'shared' ? Math.pow(f, gamma) : f;
      return `${css(viridis(t))} ${100 * f}%`;
    }).join(',');
    $('#legendGrad').style.background = `linear-gradient(90deg, ${stops})`;
    $('#legLo').textContent = domainMode === 'shared' ? '1.0' : 'low';
    $('#legHi').textContent = domainMode === 'shared' ? fmt(M.globe_ar_max, 2) : 'high';
  }

  let built = false;
  function kick() {
    if (built) return;
    built = true;
    updateLegend();
    // A and B start on two different systems for an immediate comparison.
    buildCard('A', 'h3', M.globe_res.h3.at(-1));
    buildCard('B', 'ivea7h', M.globe_res.ivea7h.at(-1));
    // Rotation+zoom sync: copy either globe's view to the other on every change.
    // setView() is idempotent, so the echoed viewchange no-ops — no guard flag.
    const a = panels.A.orb, b = panels.B.orb;
    a.on('viewchange', () => b.setView(a.getView()));
    b.on('viewchange', () => a.setView(b.getView()));
    $('#spin').addEventListener('change', (e) => {
      spinning = e.target.checked; lastT = 0;
      if (spinning) spinReq = requestAnimationFrame(spinLoop);
      else if (spinReq) cancelAnimationFrame(spinReq);
    });
    $('#domainMode').addEventListener('click', (e) => {
      const btn = e.target.closest('button'); if (!btn) return;
      domainMode = btn.dataset.v;
      $('#domainMode').querySelectorAll('button').forEach((x) => x.classList.toggle('on', x === btn));
      // Restyle in place — rewrites the per-feature style texture, no re-tessellation.
      for (const id in panels) {
        const p = panels[id];
        if (p.layer) p.layer.update({ fill: fillFn(p) });
      }
      updateLegend();
    });
  }

  return { kick };
})();

// ================= boot =================
async function main() {
  [M, HIST] = await Promise.all([
    fetch('out/manifest.json').then((r) => r.json()),
    fetch('out/histograms.json').then((r) => r.json()),
  ]);
  $('#subtitle').textContent =
    `${M.systems.length} systems · skar gap_tol ${M.gap_tol.toExponential()} · AR over 1.0`;
  buildSeriesPicker();
  initHistControls();
  // Open on the cross-system comparison at each system's working resolution.
  for (const sys of M.systems) {
    const res = M.target_res[sys];
    const chip = document.querySelector(`.chip[data-key="${sys}|${res}"]`);
    if (chip) toggleSeries(sys, res, chip);
  }
  initTabs();
}
main();
