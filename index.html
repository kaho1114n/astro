<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>アストロラーベ | ネイタルチャート計算機</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600&family=Inter:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.tailwindcss.com"></script>
<style>
  body { margin: 0; background: #14142B; }
</style>
</head>
<body>
<div id="root"></div>

<!-- React / ReactDOM (UMD build) -->
<script src="https://unpkg.com/react@18/umd/react.production.min.js" crossorigin></script>
<script src="https://unpkg.com/react-dom@18/umd/react-dom.production.min.js" crossorigin></script>
<!-- Babel standalone: ブラウザ上でJSXを変換するために使用 -->
<script src="https://unpkg.com/@babel/standalone/babel.min.js"></script>

<script type="text/babel" data-presets="react">
const { useState, useEffect, useMemo } = React;


/* ============================================================
   Astronomy module
   Base method: Paul Schlyter's low-precision Keplerian elements
   (accuracy ~1 arcmin for inner planets, ~1-2 arcmin for outer
   planets and the Moon). This version adds the main perturbation
   terms for the Moon, Jupiter, Saturn and Uranus, and uses the
   dedicated Fourier fit for Pluto, following the same reference
   (stjarnhimlen.se/comp/ppcomp.html), which improves accuracy
   beyond the plain two-body Kepler solution used previously.
   ============================================================ */

const D2R = Math.PI / 180;
const R2D = 180 / Math.PI;
const rev = (a) => a - Math.floor(a / 360) * 360;
const toRad = (a) => a * D2R;
const toDeg = (a) => a * R2D;
const clamp1 = (x) => Math.max(-1, Math.min(1, x));

const ELEMENTS = {
  sun: { N: [0, 0], i: [0, 0], w: [282.9404, 4.70935e-5], a: [1.0, 0], e: [0.016709, -1.151e-9], M: [356.047, 0.9856002585] },
  moon: { N: [125.1228, -0.0529538083], i: [5.1454, 0], w: [318.0634, 0.1643573223], a: [60.2666, 0], e: [0.0549, 0], M: [115.3654, 13.0649929509] },
  mercury: { N: [48.3313, 3.24587e-5], i: [7.0047, 5.0e-8], w: [29.1241, 1.01444e-5], a: [0.387098, 0], e: [0.205635, 5.59e-10], M: [168.6562, 4.0923344368] },
  venus: { N: [76.6799, 2.4659e-5], i: [3.3946, 2.75e-8], w: [54.891, 1.38374e-5], a: [0.72333, 0], e: [0.006773, -1.302e-9], M: [48.0052, 1.6021302244] },
  mars: { N: [49.5574, 2.11081e-5], i: [1.8497, -1.78e-8], w: [286.5016, 2.92961e-5], a: [1.523688, 0], e: [0.093405, 2.516e-9], M: [18.6021, 0.5240207766] },
  jupiter: { N: [100.4542, 2.76854e-5], i: [1.303, -1.557e-7], w: [273.8777, 1.64505e-5], a: [5.20256, 0], e: [0.048498, 4.469e-9], M: [19.895, 0.0830853001] },
  saturn: { N: [113.6634, 2.3898e-5], i: [2.4886, -1.081e-7], w: [339.3939, 2.97661e-5], a: [9.55475, 0], e: [0.055546, -9.499e-9], M: [316.967, 0.0334442282] },
  uranus: { N: [74.0005, 1.3978e-5], i: [0.7733, 1.9e-8], w: [96.6612, 3.0565e-5], a: [19.18171, -1.55e-8], e: [0.047318, 7.45e-9], M: [142.5905, 0.011725806] },
  neptune: { N: [131.7806, 3.0173e-5], i: [1.77, -2.55e-7], w: [272.8461, -6.027e-6], a: [30.05826, 3.313e-8], e: [0.008606, 2.15e-9], M: [260.2471, 0.005995147] },
};

function evalEl(pair, d) {
  return pair[0] + pair[1] * d;
}

function keplerXYZ(el, d) {
  const N = rev(evalEl(el.N, d));
  const i = evalEl(el.i, d);
  const w = rev(evalEl(el.w, d));
  const a = evalEl(el.a, d);
  const e = evalEl(el.e, d);
  const M = rev(evalEl(el.M, d));

  let E = M + R2D * e * Math.sin(toRad(M)) * (1 + e * Math.cos(toRad(M)));
  for (let k = 0; k < 12; k++) {
    const dE = (E - R2D * e * Math.sin(toRad(E)) - M) / (1 - e * Math.cos(toRad(E)));
    E -= dE;
    if (Math.abs(dE) < 1e-7) break;
  }

  const xv = a * (Math.cos(toRad(E)) - e);
  const yv = a * Math.sqrt(1 - e * e) * Math.sin(toRad(E));
  const r = Math.sqrt(xv * xv + yv * yv);
  const v = toDeg(Math.atan2(yv, xv));

  const vw = toRad(v + w);
  const Nr = toRad(N);
  const ir = toRad(i);

  const xh = r * (Math.cos(Nr) * Math.cos(vw) - Math.sin(Nr) * Math.sin(vw) * Math.cos(ir));
  const yh = r * (Math.sin(Nr) * Math.cos(vw) + Math.cos(Nr) * Math.sin(vw) * Math.cos(ir));
  const zh = r * (Math.sin(vw) * Math.sin(ir));
  return { xh, yh, zh, r };
}

// 黄経に補正量(度)を加え、緯度・距離は保ったまま直交座標を再構成する
function applyLonCorrection(xyz, deltaLonDeg, deltaLatDeg = 0) {
  const r = Math.sqrt(xyz.xh * xyz.xh + xyz.yh * xyz.yh + xyz.zh * xyz.zh);
  const lon = toDeg(Math.atan2(xyz.yh, xyz.xh)) + deltaLonDeg;
  const lat = toDeg(Math.atan2(xyz.zh, Math.sqrt(xyz.xh * xyz.xh + xyz.yh * xyz.yh))) + deltaLatDeg;
  const lonR = toRad(lon), latR = toRad(lat);
  return {
    xh: r * Math.cos(latR) * Math.cos(lonR),
    yh: r * Math.cos(latR) * Math.sin(lonR),
    zh: r * Math.sin(latR),
  };
}

function julianDay(year, month, day, hourUT) {
  let y = year, m = month;
  if (m <= 2) { y -= 1; m += 12; }
  const A = Math.floor(y / 100);
  const B = 2 - A + Math.floor(A / 4);
  return Math.floor(365.25 * (y + 4716)) + Math.floor(30.6001 * (m + 1)) + day + B - 1524.5 + hourUT / 24;
}

function plutoXYZ(d) {
  const S = 50.03 + 0.033459652 * d;
  const P = 238.95 + 0.003968789 * d;
  const s = (n) => Math.sin(toRad(n * P));
  const c = (n) => Math.cos(toRad(n * P));

  const lonecl =
    238.9508 + 0.00400703 * d
    - 19.799 * s(1) + 19.848 * c(1)
    + 0.897 * s(2) - 4.956 * c(2)
    + 0.61 * s(3) + 1.211 * c(3)
    - 0.341 * s(4) - 0.19 * c(4)
    + 0.128 * s(5) - 0.034 * c(5)
    - 0.038 * s(6) + 0.031 * c(6)
    + 0.02 * Math.sin(toRad(S - P)) - 0.01 * Math.cos(toRad(S - P));

  const latecl =
    -3.9082
    - 5.453 * s(1) - 14.975 * c(1)
    + 3.527 * s(2) + 1.673 * c(2)
    - 1.051 * s(3) + 0.328 * c(3)
    + 0.179 * s(4) - 0.292 * c(4)
    + 0.019 * s(5) + 0.1 * c(5)
    - 0.031 * s(6) - 0.026 * c(6)
    + 0.011 * Math.cos(toRad(S - P));

  const r =
    40.72
    + 6.68 * s(1) + 6.9 * c(1)
    - 1.18 * s(2) - 0.03 * c(2)
    + 0.15 * s(3) - 0.14 * c(3);

  const lonR = toRad(lonecl), latR = toRad(latecl);
  return {
    xh: r * Math.cos(latR) * Math.cos(lonR),
    yh: r * Math.cos(latR) * Math.sin(lonR),
    zh: r * Math.sin(latR),
  };
}

function computeBodies(d) {
  const sunXYZ = keplerXYZ(ELEMENTS.sun, d);
  const sunLon = rev(toDeg(Math.atan2(sunXYZ.yh, sunXYZ.xh)));
  const Ms = rev(evalEl(ELEMENTS.sun.M, d));
  const ws = rev(evalEl(ELEMENTS.sun.w, d));

  // --- Moon: base position + main perturbation terms ---
  let moonXYZ = keplerXYZ(ELEMENTS.moon, d);
  const Mm = rev(evalEl(ELEMENTS.moon.M, d));
  const Nm = rev(evalEl(ELEMENTS.moon.N, d));
  const wm = rev(evalEl(ELEMENTS.moon.w, d));
  const Ls = rev(Ms + ws);
  const Lm = rev(Mm + wm + Nm);
  const moonD = rev(Lm - Ls);
  const moonF = rev(Lm - Nm);
  const sinD = (x) => Math.sin(toRad(x));

  const moonDLon =
    -1.274 * sinD(Mm - 2 * moonD)
    + 0.658 * sinD(2 * moonD)
    - 0.186 * sinD(Ms)
    - 0.059 * sinD(2 * Mm - 2 * moonD)
    - 0.057 * sinD(Mm - 2 * moonD + Ms)
    + 0.053 * sinD(Mm + 2 * moonD)
    + 0.046 * sinD(2 * moonD - Ms)
    + 0.041 * sinD(Mm - Ms)
    - 0.035 * sinD(moonD)
    - 0.031 * sinD(Mm + Ms)
    - 0.015 * sinD(2 * moonF - 2 * moonD)
    + 0.011 * sinD(Mm - 4 * moonD);

  const moonDLat =
    -0.173 * sinD(moonF - 2 * moonD)
    - 0.055 * sinD(Mm - moonF - 2 * moonD)
    - 0.046 * sinD(Mm + moonF - 2 * moonD)
    + 0.033 * sinD(moonF + 2 * moonD)
    + 0.017 * sinD(2 * Mm + moonF);

  moonXYZ = applyLonCorrection(moonXYZ, moonDLon, moonDLat);
  const moonLon = rev(toDeg(Math.atan2(moonXYZ.yh, moonXYZ.xh)));

  // --- Mercury, Venus, Mars: no significant perturbation (per reference) ---
  const mercuryXYZ = keplerXYZ(ELEMENTS.mercury, d);
  const venusXYZ = keplerXYZ(ELEMENTS.venus, d);
  const marsXYZ = keplerXYZ(ELEMENTS.mars, d);

  // --- Jupiter, Saturn, Uranus: base + mutual perturbation terms ---
  const Mj = rev(evalEl(ELEMENTS.jupiter.M, d));
  const Msa = rev(evalEl(ELEMENTS.saturn.M, d));
  const Mu = rev(evalEl(ELEMENTS.uranus.M, d));
  const sn = (x) => Math.sin(toRad(x));
  const cs = (x) => Math.cos(toRad(x));

  let jupiterXYZ = keplerXYZ(ELEMENTS.jupiter, d);
  const jupDLon =
    -0.332 * sn(2 * Mj - 5 * Msa - 67.6)
    - 0.056 * sn(2 * Mj - 2 * Msa + 21)
    + 0.042 * sn(3 * Mj - 5 * Msa + 21)
    - 0.036 * sn(Mj - 2 * Msa)
    + 0.022 * cs(Mj - Msa)
    + 0.023 * sn(2 * Mj - 3 * Msa + 52)
    - 0.016 * sn(Mj - 5 * Msa - 69);
  jupiterXYZ = applyLonCorrection(jupiterXYZ, jupDLon);

  let saturnXYZ = keplerXYZ(ELEMENTS.saturn, d);
  const satDLon =
    0.812 * sn(2 * Mj - 5 * Msa - 67.6)
    - 0.229 * cs(2 * Mj - 4 * Msa - 2)
    + 0.119 * sn(Mj - 2 * Msa - 3)
    + 0.046 * sn(2 * Mj - 6 * Msa - 69)
    + 0.014 * sn(Mj - 3 * Msa + 32);
  const satDLat =
    -0.02 * cs(2 * Mj - 4 * Msa - 2)
    + 0.018 * sn(2 * Mj - 6 * Msa - 49);
  saturnXYZ = applyLonCorrection(saturnXYZ, satDLon, satDLat);

  let uranusXYZ = keplerXYZ(ELEMENTS.uranus, d);
  const uraDLon =
    0.04 * sn(Msa - 2 * Mu + 6)
    + 0.035 * sn(Msa - 3 * Mu + 33)
    - 0.015 * sn(Mj - Mu + 20);
  uranusXYZ = applyLonCorrection(uranusXYZ, uraDLon);

  const neptuneXYZ = keplerXYZ(ELEMENTS.neptune, d);
  const plutoXYZv = plutoXYZ(d);

  const helio = {
    mercury: mercuryXYZ, venus: venusXYZ, mars: marsXYZ,
    jupiter: jupiterXYZ, saturn: saturnXYZ, uranus: uranusXYZ,
    neptune: neptuneXYZ, pluto: plutoXYZv,
  };

  const geoLon = {};
  Object.keys(helio).forEach((k) => {
    const gx = helio[k].xh + sunXYZ.xh;
    const gy = helio[k].yh + sunXYZ.yh;
    geoLon[k] = rev(toDeg(Math.atan2(gy, gx)));
  });

  return { sunLon, moonLon, geoLon, sunXYZ, helio };
}

function computeRetrograde(d, geoLon) {
  const prev = computeBodies(d - 1);
  const retro = {};
  Object.keys(geoLon).forEach((k) => {
    if (k === "sun" || k === "moon") { retro[k] = false; return; }
    let delta = geoLon[k] - prev.geoLon[k];
    if (delta > 180) delta -= 360;
    if (delta < -180) delta += 360;
    retro[k] = delta < 0;
  });
  return retro;
}

/* ============================================================
   Houses: Placidus (iterative) with Equal-house fallback
   ============================================================ */

// プラシダス方式: 各カスプが「その天体自身の半昼弧(または半夜弧)」を
// 3等分する時間位置にあるという定義を、黄経λについて不動点反復で解く。
// H0(赤緯)=昼弧の半分=acos(-tanφ・tanδ)、夜弧の半分=acos(+tanφ・tanδ)。
// 11室・12室はMCから、2室・3室はICから、それぞれ算出する。
function solveCuspLongitude(raRef, frac, side, lat, eps) {
  let lambda = raRef;
  for (let iter = 0; iter < 40; iter++) {
    const dec = Math.asin(Math.sin(eps) * Math.sin(toRad(lambda)));
    const tanPhi = Math.tan(toRad(lat));
    const tanDec = Math.tan(dec);
    const cosH0 = clamp1(side === "diurnal" ? -tanPhi * tanDec : tanPhi * tanDec);
    const H0 = toDeg(Math.acos(cosH0));
    const sign = side === "diurnal" ? 1 : -1;
    const raTarget = rev(raRef + sign * frac * H0);
    const lambdaNew = rev(toDeg(Math.atan2(Math.sin(toRad(raTarget)), Math.cos(toRad(raTarget)) * Math.cos(eps))));
    if (Math.abs(rev(lambdaNew - lambda + 180) - 180) < 1e-7) { lambda = lambdaNew; break; }
    lambda = lambdaNew;
  }
  return rev(lambda);
}

function computePlacidusCusps(ascLon, mcLon, RAMC, lat, obliquityDeg) {
  const eps = toRad(obliquityDeg);
  if (Math.abs(lat) >= 66) return null; // 高緯度(北緯/南緯66度以上)では破綻しやすいため不採用

  try {
    const RAIC = rev(RAMC + 180);
    const c11 = solveCuspLongitude(RAMC, 1 / 3, "diurnal", lat, eps);
    const c12 = solveCuspLongitude(RAMC, 2 / 3, "diurnal", lat, eps);
    const c3 = solveCuspLongitude(RAIC, 1 / 3, "nocturnal", lat, eps);
    const c2 = solveCuspLongitude(RAIC, 2 / 3, "nocturnal", lat, eps);

    const cusps = new Array(12);
    cusps[0] = rev(ascLon);        // house 1 (ASC)
    cusps[1] = rev(c2);            // house 2
    cusps[2] = rev(c3);            // house 3
    cusps[3] = rev(mcLon + 180);   // house 4 (IC)
    cusps[4] = rev(c11 + 180);     // house 5 (opposite house 11)
    cusps[5] = rev(c12 + 180);     // house 6 (opposite house 12)
    cusps[6] = rev(ascLon + 180);  // house 7 (DESC)
    cusps[7] = rev(c2 + 180);      // house 8 (opposite house 2)
    cusps[8] = rev(c3 + 180);      // house 9 (opposite house 3)
    cusps[9] = rev(mcLon);         // house 10 (MC)
    cusps[10] = rev(c11);          // house 11
    cusps[11] = rev(c12);          // house 12

    // 妥当性チェック: 12室分の合計が360度で、どの室も170度未満であること
    let total = 0, maxSpan = 0;
    for (let i = 0; i < 12; i++) {
      const span = rev(cusps[(i + 1) % 12] - cusps[i]);
      total += span;
      maxSpan = Math.max(maxSpan, span);
    }
    if (cusps.some((v) => Number.isNaN(v)) || Math.abs(total - 360) > 0.5 || maxSpan >= 170) return null;
    return cusps;
  } catch (e) {
    return null;
  }
}

function computeEqualCusps(ascLon) {
  return Array.from({ length: 12 }, (_, k) => rev(ascLon + k * 30));
}

function houseOfLon(lon, cusps) {
  const l = rev(lon);
  for (let i = 0; i < 12; i++) {
    const start = cusps[i];
    const span = rev(cusps[(i + 1) % 12] - start);
    const pos = rev(l - start);
    if (pos < (span === 0 ? 360 : span)) return i + 1;
  }
  return 12;
}

/* ============================================================
   Full chart computation
   ============================================================ */

function computeChart({ year, month, day, hourUT, lat, lonEast, houseSystem = "placidus" }) {
  const jd = julianDay(year, month, day, hourUT);
  const d = jd - 2451543.5;

  const bodies = computeBodies(d);
  const retro = computeRetrograde(d, bodies.geoLon);

  const obliquity = 23.4393 - 3.563e-7 * d;
  const eps = toRad(obliquity);

  const GMST0 = rev(bodies.sunLon + 180);
  const GMST = rev(GMST0 + hourUT * 15);
  const LST = rev(GMST + lonEast);
  const RAMC = LST;

  const mcLon = rev(toDeg(Math.atan2(Math.sin(toRad(RAMC)), Math.cos(toRad(RAMC)) * Math.cos(eps))));
  // 注: atan2の両引数の符号を反転させることで180度回転させ、
  // (ディセンダントではなく)正しいアセンダントを得ている
  const ascLon = rev(
    toDeg(
      Math.atan2(
        Math.cos(toRad(RAMC)),
        -(Math.sin(eps) * Math.tan(toRad(lat)) + Math.cos(eps) * Math.sin(toRad(RAMC)))
      )
    )
  );

  let cusps = null;
  let usedSystem = houseSystem;
  if (houseSystem === "placidus") {
    cusps = computePlacidusCusps(ascLon, mcLon, RAMC, lat, obliquity);
    if (!cusps) { cusps = computeEqualCusps(ascLon); usedSystem = "equal_fallback"; }
  } else {
    cusps = computeEqualCusps(ascLon);
  }

  return {
    sunLon: bodies.sunLon,
    moonLon: bodies.moonLon,
    geoLon: bodies.geoLon,
    retro,
    ascLon,
    mcLon,
    obliquity,
    jd,
    cusps,
    usedSystem,
  };
}

/* ============================================================
   Zodiac and planet reference data
   ============================================================ */

const SIGNS = [
  { name: "牡羊座", glyph: "\u2648" }, { name: "牡牛座", glyph: "\u2649" },
  { name: "双子座", glyph: "\u264A" }, { name: "蟹座", glyph: "\u264B" },
  { name: "獅子座", glyph: "\u264C" }, { name: "乙女座", glyph: "\u264D" },
  { name: "天秤座", glyph: "\u264E" }, { name: "蠍座", glyph: "\u264F" },
  { name: "射手座", glyph: "\u2650" }, { name: "山羊座", glyph: "\u2651" },
  { name: "水瓶座", glyph: "\u2652" }, { name: "魚座", glyph: "\u2653" },
];

const PLANETS_ORDER = [
  { key: "sun", name: "太陽", glyph: "\u2609", color: "#D9B44A" },
  { key: "moon", name: "月", glyph: "\u263D", color: "#C9C9E0" },
  { key: "mercury", name: "水星", glyph: "\u263F", color: "#9FBF8F" },
  { key: "venus", name: "金星", glyph: "\u2640", color: "#E5A0A8" },
  { key: "mars", name: "火星", glyph: "\u2642", color: "#C9645C" },
  { key: "jupiter", name: "木星", glyph: "\u2643", color: "#D9975A" },
  { key: "saturn", name: "土星", glyph: "\u2644", color: "#8A87A8" },
  { key: "uranus", name: "天王星", glyph: "\u2645", color: "#7FB8C4" },
  { key: "neptune", name: "海王星", glyph: "\u2646", color: "#7A8FD9" },
  { key: "pluto", name: "冥王星", glyph: "\u2647", color: "#A87AC9" },
];

function signOf(lon) {
  const idx = Math.floor(rev(lon) / 30);
  return { ...SIGNS[idx], deg: rev(lon) % 30, idx };
}

function fmtDeg(d) {
  const deg = Math.floor(d);
  const min = Math.floor((d - deg) * 60);
  return deg + "\u00B0" + String(min).padStart(2, "0") + "'";
}

/* ============================================================
   Personality text (sun / moon / ascendant)
   ============================================================ */

const SUN_TRAITS = [
  "情熱的で行動力があり、物事に真っ先に飛び込むタイプです。自分の意志をはっきり持ち、負けず嫌いな一面もあります。",
  "着実で粘り強く、心地よさや美しいものを大切にする性質です。一度決めたことは簡単には揺らがない安定感があります。",
  "好奇心旺盛で頭の回転が速く、会話や情報のやり取りを楽しむタイプです。複数の物事を同時に進めるのが得意です。",
  "情に厚く、身近な人を守り育むことに喜びを感じるタイプです。家庭的で、感受性の豊かさが強みになります。",
  "自分らしさを堂々と表現し、周囲を明るく照らす存在感を持つタイプです。誇りを持って物事に取り組みます。",
  "観察力に優れ、細部まで丁寧に整えることを得意とするタイプです。誠実に努力を積み重ねる力があります。",
  "調和とバランス感覚に優れ、人と人との間を上手につなぐタイプです。美意識が高く、公平さを大切にします。",
  "物事の本質を見抜く探究心と、一度決めたら貫く芯の強さを持つタイプです。深い絆を大切にします。",
  "自由と広い世界への探究心を持ち、楽観的に前へ進むタイプです。哲学や旅を通じて視野を広げます。",
  "現実的で責任感が強く、コツコツと目標を積み上げるタイプです。長期的な視点で物事を計画します。",
  "独自の視点を持ち、既存の枠にとらわれない発想力があるタイプです。仲間との対等な関係を大切にします。",
  "共感力と想像力に富み、目に見えないものを感じ取る感受性を持つタイプです。優しさと柔軟さが持ち味です。",
];

const MOON_TRAITS = [
  "感情の動きが早く、思ったことをすぐ表に出す率直さがあります。安心には新しい刺激が必要です。",
  "感情の変化はゆっくりで、安定した環境やお気に入りのものに囲まれると落ち着きます。",
  "気持ちを言葉にして整理するタイプです。会話や情報のやり取りが心の栄養になります。",
  "感受性が非常に豊かで、身近な人との情緒的なつながりに安心を感じます。",
  "感情表現がドラマチックで、認められることや温かい注目が心の支えになります。",
  "感情を分析し整理してから表現するタイプです。役に立てているという実感が安心につながります。",
  "対立や不調和を避けたい気持ちが強く、周囲との調和が心の安定を左右します。",
  "感情の起伏が深く、表には出にくいものの内面では強い情熱を秘めています。",
  "感情に縛られず自由でいたい気持ちが強く、笑いや冒険心が心を軽くします。",
  "感情を表に出すよりコントロールする傾向があり、達成感が安心につながります。",
  "感情よりも客観性を重視する傾向があり、仲間との共感や自由な空気が心地よく感じられます。",
  "感受性が繊細で、周囲の空気や他者の感情を強く受け取りやすい性質があります。",
];

const ASC_TRAITS = [
  "初対面ではエネルギッシュでストレートな印象を与えます。行動が早く見えるでしょう。",
  "落ち着きがあり、穏やかで安定した第一印象を持たれやすいタイプです。",
  "話し上手で軽やかな雰囲気を持ち、知的で社交的な印象を与えます。",
  "親しみやすく、面倒見が良さそうな柔らかい第一印象を持たれます。",
  "存在感があり、堂々とした華やかな雰囲気を周囲に与えます。",
  "きちんとしていて誠実、控えめながらも信頼できそうな印象を持たれます。",
  "洗練された物腰で、感じの良い社交的な第一印象を与えます。",
  "ミステリアスで芯の強さを感じさせる、印象に残るオーラを持っています。",
  "明るくオープンで、自由な雰囲気を持つ第一印象を与えます。",
  "落ち着きと責任感を感じさせる、大人びた第一印象を持たれます。",
  "個性的でユニーク、他と違う視点を持っていそうな印象を与えます。",
  "柔らかく夢見がちな雰囲気で、優しそうな第一印象を持たれます。",
];

function buildPersonality(natal) {
  const sun = signOf(natal.sunLon);
  const moon = signOf(natal.moonLon);
  const asc = signOf(natal.ascLon);
  return {
    sun: { sign: sun, text: SUN_TRAITS[sun.idx] },
    moon: { sign: moon, text: MOON_TRAITS[moon.idx] },
    asc: { sign: asc, text: ASC_TRAITS[asc.idx] },
  };
}

/* ============================================================
   Per-planet detailed interpretation (planet meaning x sign style x house area)
   10天体すべてについて「その天体は何を象徴するか」「サインごとの表れ方の傾向」
   「どのハウス=生活領域で働きやすいか」を組み合わせて解釈文を生成する。
   120通り(10惑星x12サイン)を手書きする代わりに、天体の意味・サインの様式・
   ハウスの領域という3つの部品を掛け合わせる方式(占星術の伝統的な合成読みに近い)。
   ============================================================ */

const PLANET_THEMES = {
  sun: "「自分らしさ」や「何のために頑張るか」",
  moon: "ホッとする瞬間や、気持ちの動き方",
  mercury: "考えグセや、話し方・伝え方",
  venus: "好きなものの好みや、人との仲良くなり方",
  mars: "やる気の出し方や、動き出すときの勢い",
  jupiter: "運が向いてきたときの伸ばし方や、チャンスのつかみ方",
  saturn: "苦手なことへの向き合い方や、コツコツ頑張る力",
  uranus: "変化がほしくなるときの、自分らしさの出し方",
  neptune: "夢見がちな部分や、感じやすさ",
  pluto: "とことんハマる力や、譲れないこだわり",
};

// サインごとの「タイプ」。どの天体にも共通して使える汎用パーツ(日常的な言い回しで)
const SIGN_STYLE = [
  "まっすぐ突っ走るタイプ。勢いで動く",
  "マイペースにじっくり。心地よさを優先する",
  "フットワーク軽く、いろいろ試してみる",
  "情に厚く、身近な人を大事にする",
  "堂々と自分を出す。目立つのも好き",
  "細かいところまで気を配る、実用重視",
  "空気を読んでバランスを取る",
  "一点集中でとことん深掘りする",
  "自由に、興味の赴くまま広げていく",
  "地に足つけて、コツコツ積み上げる",
  "人と違う視点で、マイペースに独自路線",
  "感覚重視で、なんとなく察してしまう",
];

// ハウスごとの「よく出る場面」を日常語で
const HOUSE_THEME = {
  1: "自分そのものや、見た目の第一印象",
  2: "お金や持ち物、大事にしているもの",
  3: "普段のやり取りや、ちょっとした学び",
  4: "家族や、ホッとできる場所",
  5: "恋愛や趣味、好きなことを楽しむ時間",
  6: "仕事のやり方や、毎日の習慣・健康",
  7: "パートナーや周りの人との関係",
  8: "深い付き合いや、人と分け合うもの",
  9: "学びや遠出、新しい価値観との出会い",
  10: "仕事や、まわりからどう見られるか",
  11: "友達付き合いや、これからの夢",
  12: "一人の時間や、なんとなく感じること",
};

function buildPlanetInterpretation(row) {
  const style = SIGN_STYLE[row.sign.idx];
  const area = HOUSE_THEME[row.house] || HOUSE_THEME[1];
  const theme = PLANET_THEMES[row.key] || "";
  const text = row.name + "は" + theme + "を表す天体です。" + row.sign.name + "は「" + style + "」タイプなので、"
    + area + "の場面で、それが特に出やすいでしょう。";
  return { key: row.key, name: row.name, glyph: row.glyph, sign: row.sign, house: row.house, text };
}

/* ============================================================
   Transit comparison (current sky vs natal chart)
   ============================================================ */

const ASPECTS = [
  { key: "conjunction", angle: 0, orb: 6, label: "コンジャンクション(合)" },
  { key: "sextile", angle: 60, orb: 4, label: "セクスタイル" },
  { key: "square", angle: 90, orb: 5, label: "スクエア" },
  { key: "trine", angle: 120, orb: 5, label: "トライン" },
  { key: "opposition", angle: 180, orb: 6, label: "オポジション" },
];

const TRANSIT_PLANET_LIST = ["sun", "moon", "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto"];

const TRANSIT_MEANING = {
  sun: "自己表現や活力",
  moon: "感情や気分の流れ",
  mercury: "思考や会話、情報のやり取り",
  venus: "恋愛やお金、美的センス",
  mars: "行動力や情熱、勢い",
  jupiter: "幸運や拡大、チャンス",
  saturn: "試練や責任、現実的な課題",
  uranus: "変化や刺激、自由への欲求",
  neptune: "直感や想像力、曖昧さ",
  pluto: "変容や深層心理、徹底したこだわり",
};

const NATAL_TARGET_MEANING = {
  sun: "自分らしさや生き方",
  moon: "感情や心の安定",
  asc: "第一印象や行動スタイル",
  mc: "仕事や社会的な立場",
};

const NATAL_TARGET_LABEL = { sun: "太陽", moon: "月", asc: "アセンダント", mc: "MC" };
const PLANET_LABEL = Object.fromEntries(PLANETS_ORDER.map((p) => [p.key, p.name]));

const ASPECT_TONE = {
  conjunction: "重なり、エネルギーが強まる",
  sextile: "良い形で結びつき、チャンスをつかみやすい",
  trine: "調和し、物事がスムーズに進みやすい",
  square: "ぶつかり合い、葛藤や試練を感じやすい",
  opposition: "向き合い、バランスを見直す気づきが訪れやすい",
};

const PLANET_FACTOR = {
  sun: 0.3, moon: 0.2, mercury: 0.3, venus: 1.2, mars: 0.2,
  jupiter: 1.8, saturn: -1.6, uranus: -0.6, neptune: -0.5, pluto: -0.9,
};

const TARGET_WEIGHT = { sun: 1.3, moon: 1.1, asc: 1.2, mc: 1.0 };

// アスペクトの基本スコア: トライン/セクスタイルは惑星の強さに応じて常にプラス、
// スクエア/オポジションは常にマイナス、合(コンジャンクション)は惑星自身の性質に従う。
// (以前は合以外にも惑星の符号をそのまま掛けていたため、土星や天王星などの
//  「試練寄り」惑星のトライン/セクスタイルまでマイナスになってしまう不整合があったため統一)
function aspectBaseScore(planetKey, aspectKey) {
  const intensity = Math.abs(PLANET_FACTOR[planetKey]);
  if (aspectKey === "conjunction") return PLANET_FACTOR[planetKey];
  if (aspectKey === "trine" || aspectKey === "sextile") return 1 + intensity * 0.5;
  return -(1 + intensity * 0.5);
}

/* ============================================================
   Synastry (compatibility between two natal charts)
   ============================================================ */

// 相性診断では、恋愛・対人関係で重視されやすい太陽・月・金星・火星を
// やや重めに、天王星〜冥王星(世代天体)はやや軽めに評価する
const SYNASTRY_WEIGHT = {
  sun: 1.3, moon: 1.3, mercury: 1.0, venus: 1.4, mars: 1.2,
  jupiter: 1.0, saturn: 1.0, uranus: 0.8, neptune: 0.8, pluto: 0.9,
};

const REL_ASPECT_TONE = {
  conjunction: "強く重なり合い、お互いに大きな影響を与え合う組み合わせです",
  sextile: "程よい刺激があり、自然と協力し合いやすい組み合わせです",
  trine: "無理なく噛み合う、心地よい組み合わせです",
  square: "刺激的な反面、摩擦も生まれやすく、成長のきっかけになりやすい組み合わせです",
  opposition: "強く惹かれ合う一方で、価値観の違いに気づかされやすい組み合わせです",
};

function lonsOf(natal) {
  return { sun: natal.sunLon, moon: natal.moonLon, ...natal.geoLon };
}

// 合(コンジャンクション)は両惑星の性質の平均、トライン/セクスタイルは常にプラス、
// スクエア/オポジションは常にマイナス(単独の惑星版aspectBaseScoreと同じ考え方の2惑星版)
function synastryBaseScore(planetA, planetB, aspectKey) {
  const va = PLANET_FACTOR[planetA], vb = PLANET_FACTOR[planetB];
  const intensity = (Math.abs(va) + Math.abs(vb)) / 2;
  if (aspectKey === "conjunction") return (va + vb) / 2;
  if (aspectKey === "trine" || aspectKey === "sextile") return 1 + intensity * 0.5;
  return -(1 + intensity * 0.5);
}

function computeSynastry(natalA, natalB) {
  const lonsA = lonsOf(natalA);
  const lonsB = lonsOf(natalB);

  const matches = [];
  TRANSIT_PLANET_LIST.forEach((pa) => {
    TRANSIT_PLANET_LIST.forEach((pb) => {
      const diff = angleDiff(lonsA[pa], lonsB[pb]);
      ASPECTS.forEach((asp) => {
        const gap = Math.abs(diff - asp.angle);
        if (gap <= asp.orb) {
          const closeness = asp.orb - gap;
          const weight = (SYNASTRY_WEIGHT[pa] + SYNASTRY_WEIGHT[pb]) / 2;
          const score = synastryBaseScore(pa, pb, asp.key) * weight;
          matches.push({ planetA: pa, planetB: pb, aspectKey: asp.key, aspectLabel: asp.label, closeness, score });
        }
      });
    });
  });

  matches.sort((a, b) => b.closeness - a.closeness);
  const totalScore = matches.reduce((sum, m) => sum + m.score, 0);
  return { matches, totalScore };
}

function synastryBucket(score) {
  if (score >= 10) return { label: "とても好相性", color: "#D9B44A", desc: "自然と息が合い、一緒にいて安心感のある組み合わせです。" };
  if (score >= 3) return { label: "好相性", color: "#9FBF8F", desc: "お互いの良さを引き出しやすい、バランスの取れた組み合わせです。" };
  if (score > -3) return { label: "個性がぶつかる相性", color: "#8A87A8", desc: "似ている部分と違う部分がはっきりしていて、知るほど発見がある組み合わせです。" };
  if (score > -10) return { label: "刺激の多い相性", color: "#D9975A", desc: "惹かれ合う反面、摩擦も起きやすいので、対話を大切にすると関係が深まります。" };
  return { label: "課題の多い相性", color: "#C9645C", desc: "価値観の違いが大きく出やすい組み合わせです。違いを尊重する姿勢が鍵になります。" };
}

function angleDiff(a, b) {
  let diff = Math.abs(rev(a) - rev(b));
  if (diff > 180) diff = 360 - diff;
  return diff;
}

function computeAspects(natal, transit) {
  const transitLons = { sun: transit.sunLon, moon: transit.moonLon, ...transit.geoLon };
  const natalTargets = { sun: natal.sunLon, moon: natal.moonLon, asc: natal.ascLon, mc: natal.mcLon };

  const matches = [];
  TRANSIT_PLANET_LIST.forEach((p) => {
    Object.keys(natalTargets).forEach((t) => {
      const diff = angleDiff(transitLons[p], natalTargets[t]);
      ASPECTS.forEach((asp) => {
        const gap = Math.abs(diff - asp.angle);
        if (gap <= asp.orb) {
          const closeness = asp.orb - gap;
          const matchScore = aspectBaseScore(p, asp.key) * TARGET_WEIGHT[t];
          matches.push({ planet: p, target: t, aspectKey: asp.key, aspectLabel: asp.label, closeness, score: matchScore });
        }
      });
    });
  });

  matches.sort((a, b) => b.closeness - a.closeness);
  const totalScore = matches.reduce((sum, m) => sum + m.score, 0);
  return { matches, totalScore };
}

function fortuneBucket(score) {
  if (score >= 4) return { label: "絶好調", color: "#D9B44A", desc: "追い風が吹いています。積極的に動くと良い結果につながりやすい時期です。" };
  if (score >= 1) return { label: "順調", color: "#9FBF8F", desc: "落ち着いて物事を進めやすい、安定した時期です。" };
  if (score > -1) return { label: "静かな時期", color: "#8A87A8", desc: "大きな変化は少なく、現状維持や準備に向いている時期です。" };
  if (score > -4) return { label: "要注意", color: "#D9975A", desc: "衝突や誤解が起きやすいので、慎重な判断を心がけましょう。" };
  return { label: "低調", color: "#C9645C", desc: "無理をせず、休養や見直しを優先すると良い時期です。" };
}

/* ============================================================
   Transit timeline (about 1 month before/after today)
   ============================================================ */

// 動きが遅く、数週間〜数ヶ月単位でオーブに入り続ける天体のみを対象にする
// (太陽・月・水星・金星は動きが速すぎて「期間」として示す意味が薄いため対象外)
const SLOW_TRANSIT_PLANETS = ["mars", "jupiter", "saturn", "uranus", "neptune", "pluto"];

function fmtMD(date) {
  return (date.getMonth() + 1) + "/" + date.getDate();
}

// baseDateを中心に daysBefore〜daysAfter の範囲を1日刻みで走査し、
// その日ごとの「今日の運勢」と同じ計算式(全10天体 x 太陽・月・ASC・MC)でスコアを出す。
// 折れ線グラフ用のデータ系列を作るための関数。
// stepDays間隔のオフセット列を作る。0(今日)と両端(-daysBefore, +daysAfter)は
// 間引かれても必ず含める(前後非対称な範囲でも「今日」を確実に表示するため)
function buildOffsets(daysBefore, daysAfter, stepDays) {
  const offsets = [];
  for (let o = 0; o >= -daysBefore; o -= stepDays) offsets.push(o);
  offsets.reverse();
  if (daysBefore > 0 && offsets[0] !== -daysBefore) offsets.unshift(-daysBefore);
  for (let o = stepDays; o <= daysAfter; o += stepDays) offsets.push(o);
  if (daysAfter > 0 && offsets[offsets.length - 1] !== daysAfter) offsets.push(daysAfter);
  return offsets;
}

function computeDailyScoreSeries(natal, daysBefore, daysAfter, baseDate, stepDays = 1) {
  const natalTargets = { sun: natal.sunLon, moon: natal.moonLon, asc: natal.ascLon, mc: natal.mcLon };
  const series = [];
  buildOffsets(daysBefore, daysAfter, stepDays).forEach((off) => {
    const dt = new Date(baseDate.getTime() + off * 86400000);
    const jd = julianDay(dt.getUTCFullYear(), dt.getUTCMonth() + 1, dt.getUTCDate(), dt.getUTCHours() + dt.getUTCMinutes() / 60);
    const d = jd - 2451543.5;
    const bodies = computeBodies(d);
    const lons = { sun: bodies.sunLon, moon: bodies.moonLon, ...bodies.geoLon };
    let score = 0;
    TRANSIT_PLANET_LIST.forEach((p) => {
      Object.keys(natalTargets).forEach((t) => {
        const diff = angleDiff(lons[p], natalTargets[t]);
        ASPECTS.forEach((asp) => {
          const gap = Math.abs(diff - asp.angle);
          if (gap <= asp.orb) score += aspectBaseScore(p, asp.key) * TARGET_WEIGHT[t];
        });
      });
    });
    series.push({ offset: off, date: dt, score });
  });
  return series;
}

// baseDateを中心に daysBefore〜daysAfter の範囲を1日刻みで走査し、
// スロームーバーが各ネイタルポイントとアスペクトを組む期間(開始/正確な日/終了)を洗い出す
function buildTransitTimeline(natal, lat, lonEast, daysBefore, daysAfter, baseDate) {
  const natalTargets = { sun: natal.sunLon, moon: natal.moonLon, asc: natal.ascLon, mc: natal.mcLon };

  const days = [];
  for (let off = -daysBefore; off <= daysAfter; off++) {
    const dt = new Date(baseDate.getTime() + off * 86400000);
    const jd = julianDay(dt.getUTCFullYear(), dt.getUTCMonth() + 1, dt.getUTCDate(), dt.getUTCHours() + dt.getUTCMinutes() / 60);
    const d = jd - 2451543.5;
    const bodies = computeBodies(d);
    days.push({ date: dt, lons: bodies.geoLon });
  }

  const events = [];
  SLOW_TRANSIT_PLANETS.forEach((p) => {
    Object.keys(natalTargets).forEach((t) => {
      ASPECTS.forEach((asp) => {
        let windowStartIdx = null, minGap = Infinity, minGapIdx = null;
        days.forEach((day, idx) => {
          const diff = angleDiff(day.lons[p], natalTargets[t]);
          const gap = Math.abs(diff - asp.angle);
          const inOrb = gap <= asp.orb;
          if (inOrb) {
            if (windowStartIdx === null) windowStartIdx = idx;
            if (gap < minGap) { minGap = gap; minGapIdx = idx; }
          }
          const isLast = idx === days.length - 1;
          if ((!inOrb || isLast) && windowStartIdx !== null) {
            const windowEndIdx = inOrb && isLast ? idx : idx - 1;
            // 範囲の端でオーブに入ったまま/出たままの場合、「正確な日」がこの走査範囲の
            // 外側にある可能性があるため、その旨を openStart / openEnd / peakUncertain で示す
            const openStart = windowStartIdx === 0;
            const openEnd = windowEndIdx === days.length - 1;
            const peakUncertain = (openStart && minGapIdx === windowStartIdx) || (openEnd && minGapIdx === windowEndIdx);
            events.push({
              planet: p, target: t, aspectKey: asp.key, aspectLabel: asp.label,
              start: days[windowStartIdx].date, end: days[windowEndIdx].date, exact: days[minGapIdx].date,
              startOffset: windowStartIdx - daysBefore, endOffset: windowEndIdx - daysBefore, exactOffset: minGapIdx - daysBefore,
              openStart, openEnd, peakUncertain,
              score: aspectBaseScore(p, asp.key),
            });
            windowStartIdx = null; minGap = Infinity; minGapIdx = null;
          }
        });
      });
    });
  });

  events.sort((a, b) => a.startOffset - b.startOffset || a.exactOffset - b.exactOffset);
  return events;
}

/* ============================================================
   Chart wheel (SVG)
   ============================================================ */

function ChartWheel({ result }) {
  const cx = 300, cy = 300;
  const rOuter = 282, rSignRing = 240, rHouseRing = 190, rPlanetBase = 158, rInner = 70;
  const ascLon = result.ascLon;

  const relAngle = (lon) => rev(lon - ascLon);
  const toXY = (lon, radius) => {
    const theta = toRad(180 + relAngle(lon));
    return { x: cx + radius * Math.cos(theta), y: cy - radius * Math.sin(theta) };
  };

  const signBoundaries = Array.from({ length: 12 }, (_, k) => k * 30);
  const houseCusps = result.cusps;

  const bodies = [
    { key: "asc", name: "ASC", glyph: "AS", color: "#D9B44A", lon: result.ascLon },
    { key: "mc", name: "MC", glyph: "MC", color: "#D9B44A", lon: result.mcLon },
    { key: "sun", name: "太陽", glyph: "\u2609", color: "#D9B44A", lon: result.sunLon },
    { key: "moon", name: "月", glyph: "\u263D", color: "#C9C9E0", lon: result.moonLon },
    ...PLANETS_ORDER.filter((p) => p.key !== "sun" && p.key !== "moon").map((p) => ({
      ...p,
      lon: result.geoLon[p.key],
    })),
  ].sort((a, b) => rev(a.lon - ascLon) - rev(b.lon - ascLon));

  let lastAngle = -999, tier = 0;
  const placed = bodies.map((b) => {
    const ang = relAngle(b.lon);
    if (Math.abs(ang - lastAngle) < 9) tier = tier === 0 ? 1 : 0;
    else tier = 0;
    lastAngle = ang;
    const radius = b.key === "asc" || b.key === "mc" ? rHouseRing + 14 : rPlanetBase - tier * 26;
    return { ...b, ...toXY(b.lon, radius) };
  });

  return (
    <svg viewBox="0 0 600 600" className="w-full h-full">
      <defs>
        <radialGradient id="bgGlow" cx="50%" cy="50%" r="60%">
          <stop offset="0%" stopColor="#22224a" />
          <stop offset="100%" stopColor="#14142B" />
        </radialGradient>
      </defs>
      <circle cx={cx} cy={cy} r={rOuter} fill="url(#bgGlow)" stroke="#C9A227" strokeWidth="1" />
      <circle cx={cx} cy={cy} r={rSignRing} fill="none" stroke="#4a4a72" strokeWidth="1" />
      <circle cx={cx} cy={cy} r={rHouseRing} fill="none" stroke="#3a3a5e" strokeWidth="1" />
      <circle cx={cx} cy={cy} r={rInner} fill="none" stroke="#3a3a5e" strokeWidth="1" />

      {signBoundaries.map((deg, k) => {
        const p1 = toXY(deg, rSignRing);
        const p2 = toXY(deg, rOuter);
        const mid = toXY(deg + 15, (rSignRing + rOuter) / 2);
        return (
          <g key={"sign-" + k}>
            <line x1={p1.x} y1={p1.y} x2={p2.x} y2={p2.y} stroke="#C9A227" strokeWidth="0.6" opacity="0.6" />
            <text x={mid.x} y={mid.y} fill="#C9A227" fontSize="18" textAnchor="middle" dominantBaseline="middle" fontFamily="Cormorant Garamond, serif">
              {SIGNS[k].glyph}
            </text>
          </g>
        );
      })}

      {houseCusps.map((deg, k) => {
        const p1 = toXY(deg, rInner);
        const p2 = toXY(deg, rSignRing);
        const nextDeg = houseCusps[(k + 1) % 12];
        const midDeg = rev(deg + rev(nextDeg - deg) / 2);
        const mid = toXY(midDeg, rInner + 20);
        const major = k === 0 || k === 3 || k === 6 || k === 9;
        return (
          <g key={"house-" + k}>
            <line x1={p1.x} y1={p1.y} x2={p2.x} y2={p2.y} stroke={major ? "#C9A227" : "#5a5a86"} strokeWidth={major ? 1.4 : 0.5} />
            <text x={mid.x} y={mid.y} fill="#8A87A8" fontSize="11" textAnchor="middle" dominantBaseline="middle" fontFamily="Inter, sans-serif">
              {k + 1}
            </text>
          </g>
        );
      })}

      {placed.map((b) => {
        const radius = Math.hypot(b.x - cx, b.y - cy);
        const tick1 = toXY(b.lon, rSignRing);
        const tick2 = toXY(b.lon, radius + 12);
        return (
          <g key={b.key}>
            <line x1={tick1.x} y1={tick1.y} x2={tick2.x} y2={tick2.y} stroke={b.color} strokeWidth="0.5" opacity="0.5" />
            <circle cx={b.x} cy={b.y} r="13" fill="#14142B" stroke={b.color} strokeWidth="1" />
            <text x={b.x} y={b.y} fill={b.color} fontSize={b.key === "asc" || b.key === "mc" ? "10" : "14"} textAnchor="middle" dominantBaseline="middle" fontFamily="Inter, serif">
              {b.glyph}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

/* ============================================================
   Main app
   ============================================================ */

// 永続化ヘルパー: Claudeのアーティファクト環境ではwindow.storageを、
// それ以外の通常ブラウザ(このHTMLを直接開いた場合など)ではlocalStorageを使う。
// どちらも使えない/失敗する場合は静かに諦める(保存できないだけで、アプリ自体は動く)。
const persist = {
  async get(key) {
    if (typeof window !== "undefined" && window.storage) {
      try {
        const r = await window.storage.get(key, false);
        return r ? r.value : null;
      } catch (e) { return null; }
    }
    try { return window.localStorage.getItem(key); } catch (e) { return null; }
  },
  async set(key, value) {
    if (typeof window !== "undefined" && window.storage) {
      try { await window.storage.set(key, value, false); return; } catch (e) { /* fallthrough */ }
    }
    try { window.localStorage.setItem(key, value); } catch (e) { /* 保存できなくても致命的ではない */ }
  },
  async delete(key) {
    if (typeof window !== "undefined" && window.storage) {
      try { await window.storage.delete(key, false); return; } catch (e) { /* fallthrough */ }
    }
    try { window.localStorage.removeItem(key); } catch (e) { /* 保存できなくても致命的ではない */ }
  },
};

function App() {
  useEffect(() => {
    const link = document.createElement("link");
    link.href = "https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600&family=Inter:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap";
    link.rel = "stylesheet";
    document.head.appendChild(link);
    return () => document.head.removeChild(link);
  }, []);

  const [form, setForm] = useState({
    date: "1994-06-21",
    time: "12:00",
    tz: 9,
    lat: 35.6812,
    lon: 139.7671,
    place: "東京",
    houseSystem: "placidus",
  });
  const [showPersonB, setShowPersonB] = useState(false);
  const [formB, setFormB] = useState({
    date: "1994-06-21",
    time: "12:00",
    tz: 9,
    lat: 35.6812,
    lon: 139.7671,
    place: "東京",
    houseSystem: "placidus",
  });
  const [people, setPeople] = useState([]); // 登録済みの人一覧 [{id, name, date, time, tz, lat, lon, place, houseSystem}]
  const [natal, setNatal] = useState(null);
  const [natalB, setNatalB] = useState(null);
  const [synastry, setSynastry] = useState(null);
  const [transit, setTransit] = useState(null);
  const [transitTime, setTransitTime] = useState(null);
  const [rangeMode, setRangeMode] = useState("month"); // "month" | "pastYear" | "nextYear"
  const [topTab, setTopTab] = useState("fortune"); // "fortune" | "synastry"
  const [showPlanetDetails, setShowPlanetDetails] = useState(false);
  const [synastryTab, setSynastryTab] = useState("summary"); // "summary" | "chart" | "trend" | "aspects"
  const [error, setError] = useState("");

  const parseBirth = (f) => {
    const [y, m, day] = f.date.split("-").map(Number);
    const [hh, mm] = f.time.split(":").map(Number);
    if (!y || !m || !day || Number.isNaN(hh) || Number.isNaN(mm)) {
      throw new Error("生年月日・時刻の形式を確認してください");
    }
    const hourUT = hh + mm / 60 - Number(f.tz);
    return computeChart({
      year: y, month: m, day, hourUT,
      lat: Number(f.lat), lonEast: Number(f.lon),
      houseSystem: f.houseSystem,
    });
  };

  const saveInputs = async (f, fB, showB) => {
    await persist.set("natal-chart:formA", JSON.stringify(f));
    await persist.set("natal-chart:formB", JSON.stringify(fB));
    await persist.set("natal-chart:showPersonB", JSON.stringify(showB));
  };

  const persistPeople = async (list) => {
    await persist.set("natal-chart:people", JSON.stringify(list));
  };

  // 現在のフォーム内容を名前付きで登録する
  const registerPerson = (name, formData) => {
    const entry = {
      id: Date.now().toString(36) + Math.random().toString(36).slice(2, 7),
      name,
      date: formData.date, time: formData.time, tz: formData.tz,
      lat: formData.lat, lon: formData.lon, place: formData.place,
      houseSystem: formData.houseSystem,
    };
    const next = [...people, entry];
    setPeople(next);
    persistPeople(next);
  };

  const deletePerson = (id) => {
    const next = people.filter((p) => p.id !== id);
    setPeople(next);
    persistPeople(next);
  };

  // 登録済みの人のデータを、指定したフォーム(あなた/お相手)へ反映する
  const applyPersonTo = (id, setter) => {
    const p = people.find((x) => x.id === id);
    if (!p) return;
    setter({
      date: p.date, time: p.time, tz: p.tz,
      lat: p.lat, lon: p.lon, place: p.place,
      houseSystem: p.houseSystem || "placidus",
    });
  };

  // f/fB/showB を明示的に受け取ることで、起動時の読み込み直後にも
  // (Reactのstate反映を待たずに)そのまま計算できるようにしている
  const runCalc = (f, fB, showB) => {
    setError("");
    try {
      const natalResult = parseBirth(f);

      const now = new Date();
      const transitResult = computeChart({
        year: now.getUTCFullYear(),
        month: now.getUTCMonth() + 1,
        day: now.getUTCDate(),
        hourUT: now.getUTCHours() + now.getUTCMinutes() / 60 + now.getUTCSeconds() / 3600,
        lat: Number(f.lat), lonEast: Number(f.lon),
        houseSystem: f.houseSystem,
      });

      let natalBResult = null, synastryResult = null;
      if (showB) {
        natalBResult = parseBirth(fB);
        synastryResult = computeSynastry(natalResult, natalBResult);
      }

      setNatal(natalResult);
      setNatalB(natalBResult);
      setSynastry(synastryResult);
      setTransit(transitResult);
      setTransitTime(now);
      saveInputs(f, fB, showB);
    } catch (e) {
      setError(e.message || "計算に失敗しました");
      setNatal(null);
      setNatalB(null);
      setSynastry(null);
      setTransit(null);
    }
  };

  const handleCalc = () => runCalc(form, formB, showPersonB);

  // 起動時: 保存済みの入力があれば読み込んでから計算する。なければ初期値のまま計算する。
  useEffect(() => {
    let cancelled = false;
    (async () => {
      let loadedForm = form;
      let loadedFormB = formB;
      let loadedShowB = showPersonB;
      try {
        const a = await persist.get("natal-chart:formA");
        if (a) loadedForm = JSON.parse(a);
      } catch (e) { /* 保存データなし */ }
      try {
        const b = await persist.get("natal-chart:formB");
        if (b) loadedFormB = JSON.parse(b);
      } catch (e) { /* 保存データなし */ }
      try {
        const s = await persist.get("natal-chart:showPersonB");
        if (s) loadedShowB = JSON.parse(s);
      } catch (e) { /* 保存データなし */ }
      try {
        const p = await persist.get("natal-chart:people");
        if (p) setPeople(JSON.parse(p));
      } catch (e) { /* 保存データなし */ }
      if (cancelled) return;
      setForm(loadedForm);
      setFormB(loadedFormB);
      setShowPersonB(loadedShowB);
      runCalc(loadedForm, loadedFormB, loadedShowB);
    })();
    return () => { cancelled = true; };
  }, []); // eslint-disable-line

  // お相手の情報が無効化/未計算になったら、相性タブに取り残されないようにする
  useEffect(() => {
    if (!(natalB && synastry) && topTab === "synastry") setTopTab("fortune");
  }, [natalB, synastry, topTab]);

  // 表示期間: 1ヶ月は前後対称、1年は「過去だけ」「今後だけ」を選べるようにする。
  // 1年系は3日おきに間引いて計算する(グラフを見やすく・軽くするため)
  const RANGE_OPTIONS = {
    month: { label: "前後1ヶ月", before: 30, after: 30, step: 1 },
    pastYear: { label: "過去1年", before: 365, after: 0, step: 3 },
    nextYear: { label: "今後1年", before: 0, after: 365, step: 3 },
  };
  const range = RANGE_OPTIONS[rangeMode];
  const rangeBefore = range.before;
  const rangeAfter = range.after;
  const seriesStep = range.step;
  const rangeLabel = range.label;

  const timeline = useMemo(
    () => (natal && transitTime ? buildTransitTimeline(natal, Number(form.lat), Number(form.lon), rangeBefore, rangeAfter, transitTime) : []),
    [natal, transitTime, rangeBefore, rangeAfter, form.lat, form.lon]
  );
  const seriesA = useMemo(
    () => (natal && transitTime ? computeDailyScoreSeries(natal, rangeBefore, rangeAfter, transitTime, seriesStep) : []),
    [natal, transitTime, rangeBefore, rangeAfter, seriesStep]
  );
  const seriesB = useMemo(
    () => (natalB && transitTime ? computeDailyScoreSeries(natalB, rangeBefore, rangeAfter, transitTime, seriesStep) : []),
    [natalB, transitTime, rangeBefore, rangeAfter, seriesStep]
  );

  const rows = useMemo(() => {
    if (!natal) return [];
    const list = [
      { key: "sun", name: "太陽", glyph: "\u2609", lon: natal.sunLon, retro: false },
      { key: "moon", name: "月", glyph: "\u263D", lon: natal.moonLon, retro: false },
      ...PLANETS_ORDER.filter((p) => p.key !== "sun" && p.key !== "moon").map((p) => ({
        key: p.key, name: p.name, glyph: p.glyph, lon: natal.geoLon[p.key], retro: natal.retro[p.key],
      })),
    ];
    return list.map((b) => {
      const s = signOf(b.lon);
      return { ...b, sign: s, house: houseOfLon(b.lon, natal.cusps) };
    });
  }, [natal]);

  const personality = useMemo(() => (natal ? buildPersonality(natal) : null), [natal]);
  const planetInterpretations = useMemo(() => (rows.length > 0 ? rows.map(buildPlanetInterpretation) : []), [rows]);
  const fortune = useMemo(() => (natal && transit ? computeAspects(natal, transit) : null), [natal, transit]);
  const bucket = fortune ? fortuneBucket(fortune.totalScore) : null;

  // 相性の時期別スコア = あなたとお相手、それぞれの日次スコアの平均
  const relationshipSeries = useMemo(() => {
    if (seriesA.length === 0 || seriesB.length === 0 || seriesA.length !== seriesB.length) return [];
    return seriesA.map((a, idx) => ({ offset: a.offset, date: a.date, score: (a.score + seriesB[idx].score) / 2 }));
  }, [seriesA, seriesB]);

  return (
    <div className="min-h-screen w-full" style={{ background: "#14142B", color: "#EDE7DA", fontFamily: "Inter, sans-serif" }}>
      <div className="max-w-5xl mx-auto px-6 py-10">
        <header className="mb-8 border-b pb-6" style={{ borderColor: "#2c2c52" }}>
          <p style={{ color: "#9C8FC9", letterSpacing: "0.3em", fontSize: "11px" }}>NATAL CHART CALCULATOR</p>
          <h1 style={{ fontFamily: "Cormorant Garamond, serif", color: "#EDE7DA", fontSize: "40px", fontWeight: 600, letterSpacing: "0.02em" }}>
            アストロラーベ
          </h1>
          <p style={{ color: "#8A87A8", fontSize: "13px", marginTop: "4px" }}>
            生年月日・時刻・出生地からネイタルチャートを計算し、性格の傾向と現在の運気の流れを表示します。
          </p>
        </header>

        <div className="grid md:grid-cols-[320px_1fr] gap-8">
          <div>
            <div className="space-y-4 p-5 rounded" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
              <PersonPicker
                people={people}
                label="登録済みの人から選ぶ(あなた)"
                onLoad={(id) => applyPersonTo(id, setForm)}
                onSave={(name) => registerPerson(name, form)}
                onDelete={deletePerson}
              />
              <Field label="生年月日">
                <input type="date" value={form.date} onChange={(e) => setForm({ ...form, date: e.target.value })} style={inputStyle} />
              </Field>
              <Field label="出生時刻(現地時間)">
                <input type="time" value={form.time} onChange={(e) => setForm({ ...form, time: e.target.value })} style={inputStyle} />
              </Field>
              <Field label="タイムゾーン(UTC+)">
                <input type="number" step="0.5" value={form.tz} onChange={(e) => setForm({ ...form, tz: e.target.value })} style={inputStyle} />
              </Field>
              <Field label="出生地(参考ラベル)">
                <input type="text" value={form.place} onChange={(e) => setForm({ ...form, place: e.target.value })} style={inputStyle} />
              </Field>
              <div className="grid grid-cols-2 gap-3">
                <Field label="緯度">
                  <input type="number" step="0.0001" value={form.lat} onChange={(e) => setForm({ ...form, lat: e.target.value })} style={inputStyle} />
                </Field>
                <Field label="経度">
                  <input type="number" step="0.0001" value={form.lon} onChange={(e) => setForm({ ...form, lon: e.target.value })} style={inputStyle} />
                </Field>
              </div>
              <Field label="ハウスシステム">
                <select value={form.houseSystem} onChange={(e) => setForm({ ...form, houseSystem: e.target.value })} style={inputStyle}>
                  <option value="placidus">プラシダス</option>
                  <option value="equal">イコールハウス</option>
                </select>
              </Field>
              <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.5 }}>
                緯度・経度はGoogleマップ等で出生地を検索し取得してください(東経・北緯はプラス、西経・南緯はマイナス)。
              </p>

              <label className="flex items-center gap-2" style={{ fontSize: "12px", color: "#9C8FC9", cursor: "pointer" }}>
                <input type="checkbox" checked={showPersonB} onChange={(e) => setShowPersonB(e.target.checked)} />
                お相手との相性も見る
              </label>

              {showPersonB && (
                <div className="space-y-4 p-4 rounded" style={{ background: "#14142B", border: "1px dashed #3a3a5e" }}>
                  <p style={{ color: "#D9B44A", fontSize: "12px", fontWeight: 600 }}>お相手の情報</p>
                  <PersonPicker
                    people={people}
                    label="登録済みの人から選ぶ(お相手)"
                    onLoad={(id) => applyPersonTo(id, setFormB)}
                    onSave={(name) => registerPerson(name, formB)}
                    onDelete={deletePerson}
                  />
                  <Field label="生年月日">
                    <input type="date" value={formB.date} onChange={(e) => setFormB({ ...formB, date: e.target.value })} style={inputStyle} />
                  </Field>
                  <Field label="出生時刻(現地時間)">
                    <input type="time" value={formB.time} onChange={(e) => setFormB({ ...formB, time: e.target.value })} style={inputStyle} />
                  </Field>
                  <Field label="タイムゾーン(UTC+)">
                    <input type="number" step="0.5" value={formB.tz} onChange={(e) => setFormB({ ...formB, tz: e.target.value })} style={inputStyle} />
                  </Field>
                  <Field label="出生地(参考ラベル)">
                    <input type="text" value={formB.place} onChange={(e) => setFormB({ ...formB, place: e.target.value })} style={inputStyle} />
                  </Field>
                  <div className="grid grid-cols-2 gap-3">
                    <Field label="緯度">
                      <input type="number" step="0.0001" value={formB.lat} onChange={(e) => setFormB({ ...formB, lat: e.target.value })} style={inputStyle} />
                    </Field>
                    <Field label="経度">
                      <input type="number" step="0.0001" value={formB.lon} onChange={(e) => setFormB({ ...formB, lon: e.target.value })} style={inputStyle} />
                    </Field>
                  </div>
                </div>
              )}

              <button onClick={handleCalc} style={buttonStyle}>チャートと運勢を計算する</button>
              <button
                onClick={async () => {
                  await persist.delete("natal-chart:formA");
                  await persist.delete("natal-chart:formB");
                  await persist.delete("natal-chart:showPersonB");
                  setForm({ date: "1994-06-21", time: "12:00", tz: 9, lat: 35.6812, lon: 139.7671, place: "東京", houseSystem: "placidus" });
                  setFormB({ date: "1994-06-21", time: "12:00", tz: 9, lat: 35.6812, lon: 139.7671, place: "東京", houseSystem: "placidus" });
                  setShowPersonB(false);
                }}
                style={{ width: "100%", background: "transparent", color: "#6f6c8f", fontSize: "11px", padding: "6px", border: "1px solid #3a3a5e", borderRadius: "4px", cursor: "pointer" }}
              >
                保存した入力をクリア
              </button>
              <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.5 }}>
                入力内容はこの端末のブラウザに自動で保存され、次回開いたときも引き継がれます。登録した人の一覧も同様に保存されます(このアプリを使う人だけが見られる情報で、他の人には共有されません)。
              </p>
              {error && <p style={{ color: "#C9645C", fontSize: "12px" }}>{error}</p>}
              {natal && natal.usedSystem === "equal_fallback" && (
                <p style={{ color: "#D9975A", fontSize: "11px", lineHeight: 1.5 }}>
                  緯度が高すぎるためプラシダス方式を計算できず、イコールハウスで代用しています。
                </p>
              )}
            </div>

            <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.6, marginTop: "14px" }}>
              天体位置はケプラー軌道要素に月・木星・土星・天王星の主要な摂動項、および冥王星専用のフーリエ近似式を加えた計算で、精度はおよそ1分角(0.02度)程度です。ハウスはプラシダス方式(反復計算)を採用し、高緯度地域など計算が破綻する場合は自動的にイコールハウスに切り替わります。性格解釈・運勢診断・相性診断はあくまで自己理解のための参考情報であり、断定的な予言ではありません。精密な鑑定には専門ソフトウェアや占星術師の併用を推奨します。
            </p>
          </div>

          <div>
            {natal && (
              <>
                <div className="rounded p-4 mb-6" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
                  <ChartWheel result={natal} />
                </div>

                <div className="grid sm:grid-cols-2 gap-3 mb-6">
                  <SummaryCard label="アセンダント (ASC)" value={signOf(natal.ascLon).name + " " + fmtDeg(signOf(natal.ascLon).deg)} />
                  <SummaryCard label="MC (天頂)" value={signOf(natal.mcLon).name + " " + fmtDeg(signOf(natal.mcLon).deg)} />
                </div>

                <div className="rounded overflow-hidden mb-8" style={{ border: "1px solid #2c2c52" }}>
                  <table className="w-full text-sm" style={{ borderCollapse: "collapse" }}>
                    <thead>
                      <tr style={{ background: "#1B1B3A", color: "#9C8FC9" }}>
                        <th style={thStyle}>天体</th>
                        <th style={thStyle}>星座</th>
                        <th style={thStyle}>度数</th>
                        <th style={thStyle}>ハウス</th>
                      </tr>
                    </thead>
                    <tbody style={{ fontFamily: "IBM Plex Mono, monospace" }}>
                      {rows.map((r, idx) => (
                        <tr key={r.name} style={{ background: idx % 2 === 0 ? "#17172f" : "#1B1B3A" }}>
                          <td style={tdStyle}>
                            <span style={{ marginRight: 8, color: "#D9B44A" }}>{r.glyph}</span>
                            {r.name} {r.retro && <span style={{ color: "#C9645C", fontSize: "11px" }}>R</span>}
                          </td>
                          <td style={tdStyle}><span style={{ marginRight: 6 }}>{r.sign.glyph}</span>{r.sign.name}</td>
                          <td style={tdStyle}>{fmtDeg(r.sign.deg)}</td>
                          <td style={tdStyle}>第{r.house}ハウス</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>

                {personality && (
                  <section className="mb-8">
                    <SectionTitle eyebrow="PERSONALITY" title="あなたはどんな人?" />
                    <div className="grid gap-3">
                      <TraitCard label={"太陽星座 " + personality.sun.sign.name} sub="核となる性質・生き方" text={personality.sun.text} />
                      <TraitCard label={"月星座 " + personality.moon.sign.name} sub="感情・心の反応" text={personality.moon.text} />
                      <TraitCard label={"アセンダント " + personality.asc.sign.name} sub="第一印象・外に見える姿" text={personality.asc.text} />
                    </div>

                    <button
                      onClick={() => setShowPlanetDetails((v) => !v)}
                      style={{
                        marginTop: "12px", width: "100%", background: "transparent", color: "#9C8FC9",
                        fontSize: "12px", padding: "8px", border: "1px dashed #3a3a5e", borderRadius: "4px", cursor: "pointer",
                      }}
                    >
                      {showPlanetDetails ? "▲ 天体ごとの詳しい意味を閉じる" : "▼ 10天体すべての詳しい意味を見る(どの位置に何があるとどう働くか)"}
                    </button>

                    {showPlanetDetails && (
                      <div className="mt-3">
                        <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.6, marginBottom: "10px" }}>
                          天体ごとに「何を表すか」があり、それが「どのサイン(タイプ)」に入っているか、「どのハウス(場面)」に入っているかで、出方が変わってくる、という考え方です。
                        </p>
                        <div className="grid gap-3">
                          {planetInterpretations.map((p) => (
                            <TraitCard
                              key={p.key}
                              label={p.glyph + " " + p.name + " " + p.sign.name + " ・ 第" + p.house + "ハウス"}
                              sub={PLANET_THEMES[p.key]}
                              text={p.text}
                            />
                          ))}
                        </div>
                      </div>
                    )}
                  </section>
                )}

                {(fortune || (natalB && synastry)) && (
                  <div className="flex rounded overflow-hidden mb-4" style={{ border: "1px solid #3a3a5e", width: "fit-content" }}>
                    <button
                      onClick={() => setTopTab("fortune")}
                      style={{
                        padding: "7px 16px", fontSize: "12px", cursor: "pointer",
                        background: topTab === "fortune" ? "#C9A227" : "transparent",
                        color: topTab === "fortune" ? "#14142B" : "#8A87A8",
                        fontWeight: topTab === "fortune" ? 600 : 400,
                      }}
                    >
                      運勢の流れ
                    </button>
                    {natalB && synastry && (
                      <button
                        onClick={() => setTopTab("synastry")}
                        style={{
                          padding: "7px 16px", fontSize: "12px", cursor: "pointer",
                          background: topTab === "synastry" ? "#C9A227" : "transparent",
                          color: topTab === "synastry" ? "#14142B" : "#8A87A8",
                          fontWeight: topTab === "synastry" ? 600 : 400,
                        }}
                      >
                        お相手との相性
                      </button>
                    )}
                  </div>
                )}

                {fortune && bucket && topTab === "fortune" && (
                  <section>
                    <div className="flex items-center justify-between flex-wrap gap-2" style={{ marginBottom: "4px" }}>
                      <SectionTitle eyebrow="FORTUNE" title="運勢の流れ" />
                      <div className="flex rounded overflow-hidden" style={{ border: "1px solid #3a3a5e", height: "fit-content" }}>
                        {[{ key: "month", label: "1ヶ月" }, { key: "pastYear", label: "過去1年" }, { key: "nextYear", label: "今後1年" }].map((opt) => (
                          <button
                            key={opt.key}
                            onClick={() => setRangeMode(opt.key)}
                            style={{
                              padding: "5px 12px",
                              fontSize: "11px",
                              background: rangeMode === opt.key ? "#C9A227" : "transparent",
                              color: rangeMode === opt.key ? "#14142B" : "#8A87A8",
                              fontWeight: rangeMode === opt.key ? 600 : 400,
                              cursor: "pointer",
                            }}
                          >
                            {opt.label}
                          </button>
                        ))}
                      </div>
                    </div>
                    <p style={{ color: "#6f6c8f", fontSize: "11px", marginBottom: "10px" }}>
                      {transitTime && (transitTime.getFullYear() + "年" + (transitTime.getMonth() + 1) + "月" + transitTime.getDate() + "日 時点の空と、あなたのネイタルチャートを比較しています。")}
                      {(rangeMode === "pastYear" || rangeMode === "nextYear") && "(1年表示は3日おきの間引きで大まかな流れを見せています)"}
                    </p>

                    {/* 今日のスナップショット */}
                    <div className="rounded p-5 mb-4" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
                      <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>今日</p>
                      <div className="flex items-center justify-between mb-2">
                        <span style={{ fontFamily: "Cormorant Garamond, serif", fontSize: "26px", color: bucket.color, fontWeight: 600 }}>
                          {bucket.label}
                        </span>
                        <span style={{ fontFamily: "IBM Plex Mono, monospace", fontSize: "12px", color: "#8A87A8" }}>
                          score {fortune.totalScore.toFixed(1)}
                        </span>
                      </div>
                      <Gauge score={fortune.totalScore} color={bucket.color} />
                      <p style={{ color: "#EDE7DA", fontSize: "13px", marginTop: "10px", lineHeight: 1.6 }}>{bucket.desc}</p>
                    </div>

                    {/* 前後1ヶ月の視覚タイムライン */}
                    {timeline.length > 0 && transitTime && (
                      <>
                        <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>
                          {rangeLabel}のタイムライン(緑=追い風、赤=試練寄り、●=最も強まる日)
                        </p>
                        <div className="mb-4">
                          <TimelineChart events={timeline} daysBefore={rangeBefore} daysAfter={rangeAfter} baseDate={transitTime} />
                        </div>
                      </>
                    )}

                    {/* 前後の運勢スコア推移(折れ線グラフ) */}
                    {seriesA.length > 0 && (
                      <>
                        <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>
                          運勢スコアの推移(折れ線グラフ)
                        </p>
                        <div className="rounded p-3 mb-4" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
                          <ScoreLineChart lines={[{ key: "scoreA", label: "運勢スコア", color: "#D9B44A", series: seriesA }]} />
                        </div>
                      </>
                    )}

                    <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.6, marginBottom: "10px" }}>
                      各項目のスコアは、天体どうしの角度(アスペクト)を点数化したものです。トライン/セクスタイルはプラス、スクエア/オポジションはマイナスに働きます。木星・金星は影響が大きく、土星・天王星・冥王星はマイナス方向に強く働きやすい一方、太陽・月・水星は影響がおだやかです。ぴったりの角度(オーブが小さい)ほど1件あたりの影響も強くなります。
                    </p>

                    {/* 今日のアスペクト内訳(速い天体を含む) */}
                    <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>今日のアスペクト</p>
                    <div className="space-y-2 mb-6">
                      {fortune.matches.slice(0, 6).map((m, idx) => (
                        <div key={idx} className="rounded p-3" style={{ background: "#17172f", border: "1px solid #2c2c52" }}>
                          <div className="flex items-baseline justify-between" style={{ marginBottom: "3px" }}>
                            <span style={{ fontSize: "12px", color: m.score >= 0 ? "#9FBF8F" : "#C9645C" }}>
                              トランジット{PLANET_LABEL[m.planet]} {m.aspectLabel} ネイタル{NATAL_TARGET_LABEL[m.target]}
                            </span>
                            <span style={{ fontSize: "11px", color: m.score >= 0 ? "#9FBF8F" : "#C9645C", fontFamily: "IBM Plex Mono, monospace" }}>
                              {m.score >= 0 ? "+" : ""}{m.score.toFixed(1)}
                            </span>
                          </div>
                          <p style={{ fontSize: "13px", color: "#EDE7DA", lineHeight: 1.6 }}>
                            {TRANSIT_MEANING[m.planet]}が、あなたの{NATAL_TARGET_MEANING[m.target]}と{ASPECT_TONE[m.aspectKey]}時期です。
                          </p>
                        </div>
                      ))}
                      {fortune.matches.length === 0 && (
                        <p style={{ color: "#8A87A8", fontSize: "13px" }}>現在、特に目立った天体の重なりはなく、穏やかな時期です。</p>
                      )}
                    </div>

                    {/* 前後1ヶ月の詳細リスト(開始日順) */}
                    {timeline.length > 0 && (
                      <>
                        <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>
                          {rangeLabel}の主な動き(動きの遅い火星〜冥王星のみ、開始日順)
                        </p>
                        <div className="space-y-2">
                          {timeline.map((e, idx) => {
                            const ongoing = e.startOffset <= 0 && e.endOffset >= 0;
                            const statusTag = ongoing ? "進行中" : e.startOffset > 0 ? "今後" : "終了";
                            const statusColor = ongoing ? "#D9B44A" : e.startOffset > 0 ? "#7FB8C4" : "#6f6c8f";
                            return (
                              <div key={idx} className="rounded p-3" style={{ background: "#17172f", border: "1px solid #2c2c52" }}>
                                <div className="flex items-baseline justify-between" style={{ marginBottom: "3px" }}>
                                  <span style={{ fontSize: "12px" }}>
                                    <span style={{ color: statusColor, fontSize: "10px", border: "1px solid " + statusColor, borderRadius: "3px", padding: "1px 4px", marginRight: "6px" }}>
                                      {statusTag}
                                    </span>
                                    <span style={{ color: e.score >= 0 ? "#9FBF8F" : "#C9645C" }}>
                                      トランジット{PLANET_LABEL[e.planet]} {e.aspectLabel} ネイタル{NATAL_TARGET_LABEL[e.target]}
                                    </span>
                                  </span>
                                  <span style={{ fontSize: "11px", color: e.score >= 0 ? "#9FBF8F" : "#C9645C", fontFamily: "IBM Plex Mono, monospace" }}>
                                    {e.score >= 0 ? "+" : ""}{e.score.toFixed(1)}
                                  </span>
                                </div>
                                <p style={{ fontSize: "12px", color: "#8A87A8", marginBottom: "4px", fontFamily: "IBM Plex Mono, monospace" }}>
                                  {(e.openStart ? "〜" : "") + fmtMD(e.start) + " 〜 " + fmtMD(e.end) + (e.openEnd ? "〜" : "")}
                                  {!e.peakUncertain && "  最も強まる: " + fmtMD(e.exact)}
                                </p>
                                <p style={{ fontSize: "13px", color: "#EDE7DA", lineHeight: 1.6 }}>
                                  {TRANSIT_MEANING[e.planet]}が、あなたの{NATAL_TARGET_MEANING[e.target]}と{ASPECT_TONE[e.aspectKey]}時期です。
                                </p>
                              </div>
                            );
                          })}
                        </div>
                      </>
                    )}
                  </section>
                )}

                {natalB && synastry && topTab === "synastry" && (() => {
                  const bucket2 = synastryBucket(synastry.totalScore);
                  const personA_asc = signOf(natal.ascLon);
                  const personB_sun = signOf(natalB.sunLon);
                  const personB_moon = signOf(natalB.moonLon);
                  const personB_asc = signOf(natalB.ascLon);
                  const SYNASTRY_TABS = [
                    { key: "summary", label: "サマリー" },
                    { key: "chart", label: "お相手のチャート" },
                    { key: "trend", label: "推移グラフ" },
                    { key: "aspects", label: "アスペクト一覧" },
                  ];
                  return (
                    <section className="mt-8">
                      <SectionTitle eyebrow="SYNASTRY" title="お相手との相性" />

                      <div className="flex rounded overflow-hidden mb-4 flex-wrap" style={{ border: "1px solid #3a3a5e", width: "fit-content" }}>
                        {SYNASTRY_TABS.map((t) => (
                          <button
                            key={t.key}
                            onClick={() => setSynastryTab(t.key)}
                            style={{
                              padding: "6px 12px",
                              fontSize: "11px",
                              background: synastryTab === t.key ? "#C9A227" : "transparent",
                              color: synastryTab === t.key ? "#14142B" : "#8A87A8",
                              fontWeight: synastryTab === t.key ? 600 : 400,
                              cursor: "pointer",
                              whiteSpace: "nowrap",
                            }}
                          >
                            {t.label}
                          </button>
                        ))}
                      </div>

                      {synastryTab === "summary" && (
                        <>
                          <div className="grid sm:grid-cols-3 gap-3 mb-4">
                            <SummaryCard label="お相手の太陽星座" value={personB_sun.name} />
                            <SummaryCard label="お相手の月星座" value={personB_moon.name} />
                            <SummaryCard label="お相手のASC" value={personB_asc.name} />
                          </div>

                          <div className="rounded p-5 mb-4" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
                            <div className="flex items-center justify-between mb-2">
                              <span style={{ fontFamily: "Cormorant Garamond, serif", fontSize: "26px", color: bucket2.color, fontWeight: 600 }}>
                                {bucket2.label}
                              </span>
                              <span style={{ fontFamily: "IBM Plex Mono, monospace", fontSize: "12px", color: "#8A87A8" }}>
                                score {synastry.totalScore.toFixed(1)}
                              </span>
                            </div>
                            <Gauge score={synastry.totalScore / 2.2} color={bucket2.color} />
                            <p style={{ color: "#EDE7DA", fontSize: "13px", marginTop: "10px", lineHeight: 1.6 }}>{bucket2.desc}</p>
                          </div>

                          <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.5 }}>
                            このスコアは2人の生まれ持った相性(不変)です。時期による運気の重なりは「推移グラフ」タブ、個別のアスペクトは「アスペクト一覧」タブをご覧ください。
                          </p>
                        </>
                      )}

                      {synastryTab === "chart" && (
                        <div className="rounded p-4" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
                          <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "8px", letterSpacing: "0.05em" }}>
                            お相手のネイタルチャート{formB.place ? "(" + formB.place + ")" : ""}
                          </p>
                          <ChartWheel result={natalB} />
                        </div>
                      )}

                      {synastryTab === "trend" && (
                        <>
                          {relationshipSeries.length > 0 ? (
                            <>
                              <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>
                                {rangeLabel}の運気の重なり(折れ線グラフ)
                              </p>
                              <div className="rounded p-3 mb-3" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
                                <ScoreLineChart
                                  lines={[
                                    { key: "you", label: "あなた", color: "#7FB8C4", series: seriesA },
                                    { key: "partner", label: "お相手", color: "#E5A0A8", series: seriesB },
                                    { key: "rel", label: "相性(平均)", color: "#D9B44A", series: relationshipSeries },
                                  ]}
                                  height={210}
                                />
                              </div>
                              <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.5 }}>
                                金色の線(相性)が高い日ほど2人にとって動きやすいタイミング、低い日は無理をしない方が良いタイミングの目安です。期間の切替は上の「運勢の流れ」セクションのトグルと連動しています。
                              </p>
                            </>
                          ) : (
                            <p style={{ color: "#8A87A8", fontSize: "13px" }}>推移データを計算中、またはデータがありません。</p>
                          )}
                        </>
                      )}

                      {synastryTab === "aspects" && (
                        <>
                          <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.6, marginBottom: "10px" }}>
                            あなたと お相手の10天体どうし(合計100通り)の角度を1件ずつ点数化し、合計したものが総合スコアです。トライン/セクスタイルはプラス、スクエア/オポジションはマイナスに働きます。太陽・月・金星・火星の組み合わせは影響が大きく、天王星〜冥王星の組み合わせはやや控えめに扱っています。以下は上位8件です。
                          </p>
                          <div className="space-y-2">
                            {synastry.matches.slice(0, 8).map((m, idx) => (
                              <div key={idx} className="rounded p-3" style={{ background: "#17172f", border: "1px solid #2c2c52" }}>
                                <div className="flex items-baseline justify-between" style={{ marginBottom: "3px" }}>
                                  <span style={{ fontSize: "12px", color: m.score >= 0 ? "#9FBF8F" : "#C9645C" }}>
                                    あなたの{PLANET_LABEL[m.planetA]} {m.aspectLabel} お相手の{PLANET_LABEL[m.planetB]}
                                  </span>
                                  <span style={{ fontSize: "11px", color: m.score >= 0 ? "#9FBF8F" : "#C9645C", fontFamily: "IBM Plex Mono, monospace" }}>
                                    {m.score >= 0 ? "+" : ""}{m.score.toFixed(1)}
                                  </span>
                                </div>
                                <p style={{ fontSize: "13px", color: "#EDE7DA", lineHeight: 1.6 }}>
                                  あなたの{TRANSIT_MEANING[m.planetA]}と、お相手の{TRANSIT_MEANING[m.planetB]}が{REL_ASPECT_TONE[m.aspectKey]}
                                </p>
                              </div>
                            ))}
                            {synastry.matches.length === 0 && (
                              <p style={{ color: "#8A87A8", fontSize: "13px" }}>目立ったアスペクトは見つかりませんでした。オーブの狭い、控えめな結びつきの組み合わせです。</p>
                            )}
                          </div>
                          <p style={{ color: "#6f6c8f", fontSize: "11px", lineHeight: 1.6, marginTop: "10px" }}>
                            {"参考: あなたのアセンダントは" + personA_asc.name + "です。"}
                          </p>
                        </>
                      )}
                    </section>
                  );
                })()}
              </>

            )}
          </div>
        </div>
      </div>
    </div>
  );
}

function Field({ label, children }) {
  return (
    <label className="block">
      <span style={{ display: "block", fontSize: "11px", color: "#9C8FC9", marginBottom: "4px", letterSpacing: "0.05em" }}>{label}</span>
      {children}
    </label>
  );
}

// 登録済みの人から選んでフォームに反映する/現在の内容を新規登録する/削除する、ためのUI部品
function PersonPicker({ people, label, onLoad, onSave, onDelete }) {
  const [name, setName] = useState("");
  return (
    <div className="p-3 rounded" style={{ background: "#14142B", border: "1px dashed #3a3a5e" }}>
      <p style={{ color: "#9C8FC9", fontSize: "11px", marginBottom: "6px", letterSpacing: "0.05em" }}>{label}</p>

      {people.length > 0 && (
        <select
          defaultValue=""
          onChange={(e) => {
            if (e.target.value) { onLoad(e.target.value); e.target.value = ""; }
          }}
          style={{ ...inputStyle, marginBottom: "8px" }}
        >
          <option value="">登録済みから選んで反映...</option>
          {people.map((p) => (
            <option key={p.id} value={p.id}>{p.name}({p.place || "場所未設定"})</option>
          ))}
        </select>
      )}

      <div className="flex gap-2">
        <input
          type="text"
          placeholder="名前を付けて登録(例: 母、友人A)"
          value={name}
          onChange={(e) => setName(e.target.value)}
          style={{ ...inputStyle, flex: 1 }}
        />
        <button
          onClick={() => { if (name.trim()) { onSave(name.trim()); setName(""); } }}
          style={{ padding: "0 12px", background: "#C9A227", color: "#14142B", borderRadius: "4px", fontSize: "12px", fontWeight: 600, cursor: "pointer", whiteSpace: "nowrap" }}
        >
          今の内容を登録
        </button>
      </div>

      {people.length > 0 && (
        <div className="flex flex-wrap gap-2" style={{ marginTop: "8px" }}>
          {people.map((p) => (
            <span
              key={p.id}
              style={{ fontSize: "11px", color: "#8A87A8", background: "#1B1B3A", padding: "3px 6px 3px 10px", borderRadius: "999px", display: "inline-flex", alignItems: "center", gap: "6px", border: "1px solid #2c2c52" }}
            >
              {p.name}
              <button
                onClick={() => onDelete(p.id)}
                title="削除"
                style={{ color: "#C9645C", cursor: "pointer", background: "none", border: "none", padding: "0 4px", fontSize: "13px", lineHeight: 1 }}
              >
                ×
              </button>
            </span>
          ))}
        </div>
      )}
    </div>
  );
}

function SummaryCard({ label, value }) {
  return (
    <div className="rounded p-4" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
      <p style={{ color: "#8A87A8", fontSize: "11px", letterSpacing: "0.05em" }}>{label}</p>
      <p style={{ color: "#D9B44A", fontFamily: "Cormorant Garamond, serif", fontSize: "22px", marginTop: "2px" }}>{value}</p>
    </div>
  );
}

function SectionTitle({ eyebrow, title }) {
  return (
    <div className="mb-4">
      <p style={{ color: "#9C8FC9", letterSpacing: "0.25em", fontSize: "10px" }}>{eyebrow}</p>
      <h2 style={{ fontFamily: "Cormorant Garamond, serif", fontSize: "24px", color: "#EDE7DA" }}>{title}</h2>
    </div>
  );
}

function TraitCard({ label, sub, text }) {
  return (
    <div className="rounded p-4" style={{ background: "#1B1B3A", border: "1px solid #2c2c52" }}>
      <div className="flex items-baseline justify-between mb-1">
        <span style={{ color: "#D9B44A", fontSize: "14px", fontWeight: 600 }}>{label}</span>
        <span style={{ color: "#6f6c8f", fontSize: "11px" }}>{sub}</span>
      </div>
      <p style={{ color: "#EDE7DA", fontSize: "13px", lineHeight: 1.7 }}>{text}</p>
    </div>
  );
}

function Gauge({ score, color }) {
  const pct = Math.max(0, Math.min(100, ((score + 6) / 12) * 100));
  return (
    <div style={{ background: "#0f0f26", borderRadius: "999px", height: "8px", overflow: "hidden" }}>
      <div style={{ width: pct + "%", height: "100%", background: color, transition: "width 0.4s ease" }} />
    </div>
  );
}

// 前後n日のアスペクトを、横軸=日付の帯グラフとして可視化する。
// 同じ時期に複数のアスペクトが重なる場合は自動的にレーン(段)を分けて重なりを避ける。
function TimelineChart({ events, daysBefore, daysAfter, baseDate }) {
  const total = daysBefore + daysAfter;
  const pct = (offset) => ((offset + daysBefore) / total) * 100;
  const rangeStart = new Date(baseDate.getTime() - daysBefore * 86400000);
  const rangeEnd = new Date(baseDate.getTime() + daysAfter * 86400000);

  const sorted = [...events].sort((a, b) => a.startOffset - b.startOffset);
  const laneEnds = [];
  const laned = sorted.map((e) => {
    let lane = laneEnds.findIndex((end) => end < e.startOffset);
    if (lane === -1) { lane = laneEnds.length; laneEnds.push(e.endOffset); }
    else laneEnds[lane] = e.endOffset;
    return { ...e, lane };
  });
  const laneCount = Math.max(laneEnds.length, 1);
  const laneHeight = 20;
  const topPad = 8;
  const bottomPad = 22;
  const height = laneCount * laneHeight + topPad + bottomPad;
  const todayPct = pct(0);

  return (
    <div style={{ position: "relative", height: height + "px", background: "#14142B", borderRadius: "6px", border: "1px solid #2c2c52", padding: "0 12px" }}>
      <div style={{ position: "absolute", left: todayPct + "%", top: 4, bottom: 18, width: "1px", background: "#C9A227", opacity: 0.7 }} />
      <div style={{ position: "absolute", left: todayPct + "%", bottom: 2, transform: "translateX(-50%)", fontSize: "10px", color: "#C9A227", whiteSpace: "nowrap" }}>今日</div>
      {daysBefore > 0 && <div style={{ position: "absolute", left: 12, bottom: 2, fontSize: "10px", color: "#6f6c8f" }}>{fmtMD(rangeStart)}</div>}
      {daysAfter > 0 && <div style={{ position: "absolute", right: 12, bottom: 2, fontSize: "10px", color: "#6f6c8f" }}>{fmtMD(rangeEnd)}</div>}
      {laned.map((e, idx) => {
        const left = pct(e.startOffset);
        const width = Math.max(pct(e.endOffset) - pct(e.startOffset), 1.5);
        const exactLeft = pct(e.exactOffset);
        const color = e.score >= 0 ? "#9FBF8F" : "#C9645C";
        const top = e.lane * laneHeight + topPad;
        return (
          <React.Fragment key={idx}>
            <div style={{ position: "absolute", left: left + "%", width: width + "%", top: top + 6 + "px", height: "5px", borderRadius: "3px", background: color, opacity: 0.5 }} />
            {!e.peakUncertain && (
              <div style={{ position: "absolute", left: exactLeft + "%", top: top + 2 + "px", width: "9px", height: "9px", marginLeft: "-4.5px", borderRadius: "50%", background: color, border: "1px solid #14142B" }} />
            )}
          </React.Fragment>
        );
      })}
    </div>
  );
}

// 日次スコア推移(computeDailyScoreSeriesの出力)を折れ線グラフで表示する。
// lines: [{ key, label, color, series }] を渡すと複数系列を重ねて描画できる
// (個別運勢は1本、相性は本人・お相手の2本、を想定)
function ScoreLineChart({ lines, height = 190 }) {
  const base = lines[0].series;
  const n = base.length;
  const data = base.map((pt, idx) => {
    const row = { offset: pt.offset, label: fmtMD(pt.date) };
    lines.forEach((l) => { row[l.key] = Number(l.series[idx].score.toFixed(2)); });
    return row;
  });
  const todayIdx = data.findIndex((d) => d.offset === 0);

  const allScores = lines.flatMap((l) => l.series.map((s) => s.score));
  let min = Math.min(0, ...allScores);
  let max = Math.max(0, ...allScores);
  if (min === max) { min -= 1; max += 1; }
  const pad = (max - min) * 0.12;
  min -= pad;
  max += pad;

  const W = 600;
  const marginLeft = 30, marginRight = 10, marginTop = 14, marginBottom = 20;
  const plotW = W - marginLeft - marginRight;
  const plotH = height - marginTop - marginBottom;

  const xAt = (i) => marginLeft + (n <= 1 ? 0 : (i / (n - 1)) * plotW);
  const yAt = (v) => marginTop + (1 - (v - min) / (max - min)) * plotH;

  const tickCount = Math.min(6, n);
  const tickIdxs = [...new Set(Array.from({ length: tickCount }, (_, k) => Math.round((k / (tickCount - 1 || 1)) * (n - 1))))];
  const yTicks = [min + pad, (min + max) / 2, max - pad];

  return (
    <div>
      <svg viewBox={"0 0 " + W + " " + height} style={{ width: "100%", height: height + "px" }}>
        {yTicks.map((v, idx) => (
          <React.Fragment key={idx}>
            <line x1={marginLeft} y1={yAt(v)} x2={W - marginRight} y2={yAt(v)} stroke="#2c2c52" strokeDasharray="3 3" />
            <text x={marginLeft - 5} y={yAt(v)} fill="#8A87A8" fontSize="9" textAnchor="end" dominantBaseline="middle">{v.toFixed(0)}</text>
          </React.Fragment>
        ))}
        <line x1={marginLeft} y1={yAt(0)} x2={W - marginRight} y2={yAt(0)} stroke="#4a4a72" strokeWidth="1" />
        {todayIdx >= 0 && (
          <>
            <line x1={xAt(todayIdx)} y1={marginTop} x2={xAt(todayIdx)} y2={marginTop + plotH} stroke="#C9A227" strokeDasharray="4 2" />
            <text x={xAt(todayIdx)} y={marginTop - 3} fill="#C9A227" fontSize="9" textAnchor="middle">今日</text>
          </>
        )}
        {tickIdxs.map((i) => (
          <text key={i} x={xAt(i)} y={height - 4} fill="#8A87A8" fontSize="9" textAnchor="middle">{data[i].label}</text>
        ))}
        {lines.map((l) => {
          const points = data.map((d, i) => xAt(i) + "," + yAt(d[l.key])).join(" ");
          return <polyline key={l.key} points={points} fill="none" stroke={l.color} strokeWidth="2" />;
        })}
        {lines.map((l) => (
          <React.Fragment key={l.key + "-pts"}>
            {data.map((d, i) => (
              <circle key={i} cx={xAt(i)} cy={yAt(d[l.key])} r="6" fill="transparent">
                <title>{d.label + " " + l.label + ": " + d[l.key]}</title>
              </circle>
            ))}
          </React.Fragment>
        ))}
      </svg>
      {lines.length > 1 && (
        <div className="flex gap-4 justify-center flex-wrap" style={{ fontSize: "11px", color: "#8A87A8", marginTop: "2px" }}>
          {lines.map((l) => (
            <span key={l.key} style={{ display: "inline-flex", alignItems: "center", gap: "4px" }}>
              <span style={{ width: "10px", height: "10px", borderRadius: "50%", background: l.color, display: "inline-block" }} />
              {l.label}
            </span>
          ))}
        </div>
      )}
    </div>
  );
}

const inputStyle = {
  width: "100%",
  background: "#14142B",
  border: "1px solid #3a3a5e",
  borderRadius: "4px",
  padding: "8px 10px",
  color: "#EDE7DA",
  fontSize: "13px",
  fontFamily: "IBM Plex Mono, monospace",
};

const buttonStyle = {
  width: "100%",
  background: "#C9A227",
  color: "#14142B",
  fontWeight: 600,
  padding: "10px",
  borderRadius: "4px",
  fontSize: "13px",
  letterSpacing: "0.05em",
  cursor: "pointer",
};

const thStyle = { textAlign: "left", padding: "8px 12px", fontSize: "11px", letterSpacing: "0.05em", fontWeight: 500 };
const tdStyle = { padding: "8px 12px", borderTop: "1px solid #2c2c52" };


const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);
</script>
</body>
</html>
