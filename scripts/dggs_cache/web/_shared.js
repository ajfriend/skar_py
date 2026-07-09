// Helpers shared by the two viewer pages (app.js and globe_full.html) so the
// color pipeline can't drift between them. App-level code by design — the
// ajglobe library owns no color scales or data loading.

export const DNC_GREY = [68, 68, 68, 255];

// viridis via control stops — no per-cell color-string parsing.
const STOPS = [[68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142],
               [33, 144, 141], [39, 173, 129], [92, 200, 99], [253, 231, 37]];
export function viridis(t) {
  t = Math.max(0, Math.min(1, t)) * (STOPS.length - 1);
  const i = Math.min(STOPS.length - 2, t | 0), f = t - i, a = STOPS[i], b = STOPS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, 255];
}

// 256-entry LUT of preallocated RGBA arrays: fill callbacks run once per cell
// (up to 1.18M), so color lookups should be indexed, not allocated.
export const VIRIDIS = Array.from({ length: 256 }, (_, i) => viridis(i / 255));
export const lut = (t) => VIRIDIS[Math.min(255, (t * 255) | 0)];

export async function fetchBin(path, Ctor) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path} ${r.status}`);
  return new Ctor(await r.arrayBuffer());
}

// Ascending finite values of an AR array (NaN = did-not-converge dropped).
export function sortedFinite(ar) {
  return Float64Array.from(ar).filter(Number.isFinite).sort();
}

// Per-cell quantile position in [0, 1] (NaN where the cell DNC'd), computed
// once per load so fill callbacks are lookups. Rank is bisect-right into the
// sorted values — equal ARs get equal color, ties never split.
export function quantileT(ar, sorted = sortedFinite(ar)) {
  const t = new Float32Array(ar.length).fill(NaN);
  const m = sorted.length;
  for (let i = 0; i < ar.length; i++) {
    const a = ar[i];
    if (!Number.isFinite(a)) continue;
    let lo = 0, hi = m;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] <= a) lo = mid + 1; else hi = mid; }
    t[i] = m > 1 ? Math.min(1, lo / (m - 1)) : 0.5;
  }
  return t;
}
