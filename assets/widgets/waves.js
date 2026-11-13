/* =========================================================
   Ch3 地震波與波傳理論：互動元件
   ========================================================= */
document.addEventListener('sd:ready', () => {
  const TAU = Math.PI * 2, f = SD.f;
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

  function clearG(g) { while (g.firstChild) g.removeChild(g.firstChild); }
  function svgText(g, x, y, s, cls, anchor) {
    const t = svgEl('text', { x, y, class: cls || 'txt', 'text-anchor': anchor || 'start' }, g);
    svgRich(t, s);
    return t;
  }

  /* ---- Rayleigh 半空間速度比 K(ν)：解 K⁶ - 8K⁴ + (24 - 16α²)K² - 16(1 - α²) = 0 ---- */
  function rayleighK(nu) {
    const a2 = (1 - 2 * nu) / (2 * (1 - nu));
    const fn = x => x * x * x - 8 * x * x + (24 - 16 * a2) * x - 16 * (1 - a2);
    let lo = 0.5, hi = 0.99;
    for (let i = 0; i < 60; i++) {
      const m = 0.5 * (lo + hi);
      if (fn(m) < 0) lo = m; else hi = m;
    }
    return Math.sqrt(0.5 * (lo + hi));
  }

  /* ---- 共振柱頻率方程式 β tan β = r 的第一正根 ---- */
  function solveBeta1(r) {
    if (!isFinite(r) || r > 1e6) return Math.PI / 2;
    if (r <= 0) return 0;
    let a = 1e-9, b = Math.PI / 2 - 1e-9;
    for (let i = 0; i < 60; i++) {
      const m = 0.5 * (a + b);
      if (m * Math.tan(m) < r) a = m; else b = m;
    }
    return 0.5 * (a + b);
  }

  /* ---- Kausel & Roësset (1981) 層狀半空間 Rayleigh 波求解器 ---- */
  const KR = (() => {
    const GRAV = 9.81;
    const _hb = new Float64Array(8);
    function deriveModel(layers) {
      const n = layers.length - 1;
      const M = { n, h: [], vs: [], vp: [], rho: [], nu: [], gamma: [] };
      for (const L of layers) {
        M.h.push(L.h); M.vs.push(L.vs); M.nu.push(L.nu); M.gamma.push(L.gamma);
        M.vp.push(L.vs * Math.sqrt((2 - 2 * L.nu) / (1 - 2 * L.nu)));
        M.rho.push(L.gamma * 1000 / GRAV);
      }
      return M;
    }
    function hyp(x2, z, out, o) {
      if (x2 > 1e-12) {
        const x = Math.sqrt(x2), t = z * x, e = Math.exp(-2 * t), S = 0.5 * (1 - e);
        out[o] = 0.5 * (1 + e); out[o + 1] = S / x; out[o + 2] = x * S; out[o + 3] = Math.exp(-t);
      } else if (x2 < -1e-12) {
        const x = Math.sqrt(-x2), t = z * x, sn = Math.sin(t);
        out[o] = Math.cos(t); out[o + 1] = sn / x; out[o + 2] = -x * sn; out[o + 3] = 1;
      } else {
        out[o] = 1; out[o + 1] = z; out[o + 2] = x2 * z; out[o + 3] = 1;
      }
    }
    function assembleK(c, k, M, K, scale) {
      const n = M.n, N = 2 * (n + 1);
      K.fill(0);
      const c2 = c * c;
      let prodD = 1;
      for (let j = 0; j < n; j++) {
        const a = M.vp[j], b = M.vs[j], rho = M.rho[j], z = k * M.h[j];
        const r2 = 1 - c2 / (a * a), s2 = 1 - c2 / (b * b);
        hyp(r2, z, _hb, 0); hyp(s2, z, _hb, 4);
        const Cr = _hb[0], SrOr = _hb[1], rSr = _hb[2], Er = _hb[3];
        const Cs = _hb[4], SsOs = _hb[5], sSs = _hb[6], Es = _hb[7];
        const ER = Er * Es;
        const D = 2 * (ER - Cr * Cs) + SrOr * SsOs + rSr * sSs;
        const ff = k * rho * c2 / D / scale;
        const g = k * rho * b * b * (1 + s2) / scale;
        const k11 = ff * (Cr * SsOs - rSr * Cs);
        const k12 = ff * (Cr * Cs - rSr * sSs - ER) - g;
        const k13 = ff * (rSr * Es - SsOs * Er);
        const k14 = ff * (-Cr * Es + Cs * Er);
        const k22 = ff * (SrOr * Cs - Cr * sSs);
        const k24 = ff * (-SrOr * Es + sSs * Er);
        const i0 = 2 * j;
        const r0 = i0 * N, r1 = (i0 + 1) * N, rr2 = (i0 + 2) * N, r3 = (i0 + 3) * N;
        K[r0 + i0] += k11; K[r0 + i0 + 1] += k12; K[r0 + i0 + 2] += k13; K[r0 + i0 + 3] += k14;
        K[r1 + i0] += k12; K[r1 + i0 + 1] += k22; K[r1 + i0 + 2] += -k14; K[r1 + i0 + 3] += k24;
        K[rr2 + i0] += k13; K[rr2 + i0 + 1] += -k14; K[rr2 + i0 + 2] += k11; K[rr2 + i0 + 3] += -k12;
        K[r3 + i0] += k14; K[r3 + i0 + 1] += k24; K[r3 + i0 + 2] += -k12; K[r3 + i0 + 3] += k22;
        prodD *= D;
      }
      const a = M.vp[n], b = M.vs[n], rho = M.rho[n];
      const r = Math.sqrt(Math.max(0, 1 - c2 / (a * a))), s2 = 1 - c2 / (b * b), s = Math.sqrt(Math.max(0, s2));
      const m = k * rho * b * b / scale;
      const q = (1 - s2) / (1 - r * s);
      const i0 = 2 * n;
      K[i0 * N + i0] += m * r * q;
      K[i0 * N + i0 + 1] += m * q - 2 * m;
      K[(i0 + 1) * N + i0] += m * q - 2 * m;
      K[(i0 + 1) * N + i0 + 1] += m * s * q;
      return prodD;
    }
    function detLU(A, N) {
      let det = 1;
      for (let p = 0; p < N; p++) {
        let maxv = Math.abs(A[p * N + p]), piv = p;
        for (let i = p + 1; i < N; i++) { const v = Math.abs(A[i * N + p]); if (v > maxv) { maxv = v; piv = i; } }
        if (maxv === 0) return 0;
        if (piv !== p) {
          for (let j = p; j < N; j++) { const t = A[p * N + j]; A[p * N + j] = A[piv * N + j]; A[piv * N + j] = t; }
          det = -det;
        }
        const d = A[p * N + p]; det *= d;
        for (let i = p + 1; i < N; i++) {
          const fct = A[i * N + p] / d;
          if (fct !== 0) for (let j = p + 1; j < N; j++) A[i * N + j] -= fct * A[p * N + j];
        }
      }
      return det;
    }
    function refineRoot(fn, a, b, Fa, Fb) {
      let side = 0;
      for (let it = 0; it < 50; it++) {
        const c = (a * Fb - b * Fa) / (Fb - Fa);
        const Fc = fn(c);
        if (Fc === 0) return c;
        if (Fc * Fb > 0) { b = c; Fb = Fc; if (side === -1) Fa *= 0.5; side = -1; }
        else { a = c; Fa = Fc; if (side === 1) Fb *= 0.5; side = 1; }
        if (Math.abs(b - a) < 1e-8 * Math.abs(c)) return c;
      }
      return 0.5 * (a + b);
    }
    function findModes(freq, M, maxModes = 2, nc = 180) {
      const N = 2 * (M.n + 1), K = new Float64Array(N * N);
      const w = TAU * freq;
      const minVs = Math.min(...M.vs);
      const cmin = 0.65 * minVs, cmax = M.vs[M.n] * (1 - 1e-6);
      const scale = M.rho[0] * M.vs[0] * w;
      const fn = c => {
        const P = assembleK(c, w / c, M, K, scale);
        return detLU(K, N) * P;
      };
      const roots = [];
      let c0 = cmin, F0 = fn(c0);
      const dc = (cmax - cmin) / nc;
      for (let i = 1; i <= nc && roots.length < maxModes; i++) {
        const c1 = i === nc ? cmax : cmin + i * dc;
        const F1 = fn(c1);
        if (isFinite(F0) && isFinite(F1) && F0 * F1 < 0) {
          const r = refineRoot(fn, c0, c1, F0, F1);
          const Fr = Math.abs(fn(r));
          if (Fr < 1e-2 * Math.max(Math.abs(F0), Math.abs(F1))) roots.push(r);
        }
        c0 = c1; F0 = F1;
      }
      return roots;
    }
    return { deriveModel, findModes };
  })();

  /* =========================================================
     1. #w-waves：四種波的質點運動 (3.1)
     ========================================================= */
  (() => {
    const body = SD.widget('w-waves'); if (!body) return;
    let mode = 'P', playing = true, phase = 0, selCol = 5, selRow = 1;
    const ctl = SD.controls(body, [
      { k: 'amp', label: '位移放大倍率', min: 0.4, max: 1.8, step: 0.1, value: 1.0, fmt: v => v.toFixed(1) + '×' },
      { k: 'spd', label: '播放速度', min: 0.3, max: 2.0, step: 0.1, value: 1.0, fmt: v => v.toFixed(1) + '×' }
    ], render);
    SD.buttons(body, '波型', [
      { label: 'P 波（縱波）', primary: true, on: () => { mode = 'P'; render(); } },
      { label: 'SV 波（垂直剪力波）', on: () => { mode = 'SV'; render(); } },
      { label: 'SH / Love 波（水平剪力）', on: () => { mode = 'Love'; render(); } },
      { label: 'Rayleigh 波（表面逆行橢圓）', on: () => { mode = 'Rayleigh'; render(); } },
      { label: '❚❚ 暫停', on: b => { playing = !playing; b.textContent = playing ? '❚❚ 暫停' : '▶ 播放'; if (playing) requestAnimationFrame(loop); } }
    ]);
    const W = 680, H = 250;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot-svg', style: 'cursor:pointer' });
    body.appendChild(svg);
    const g = svgEl('g', null, svg);
    const st = SD.stats(body, [
      { k: 'name', label: '目前波型', hl: true },
      { k: 'dir', label: '質點運動方向' },
      { k: 'path', label: '單一質點軌跡', hl: true },
      { k: 'vol', label: '體積應變 ε_v' },
      { k: 'sh', label: '剪應變 γ' },
      { k: 'decay', label: '隨深度變化' }
    ]);
    SD.div(body, 'w-note', '點選網格上任一顆質點可切換藍色追蹤框：紅線標示「等相位波峰」正以波速 V 往右傳播，但藍框裡的質點始終只在自己的平衡位置附近來回運動。');

    const NX = 15, NY = 5, x0 = 55, x1 = 625, y0 = 48, y1 = 205;
    const dx = (x1 - x0) / (NX - 1), dy = (y1 - y0) / (NY - 1);

    function dispAt(ix, iy, ph) {
      const k = TAU / ((x1 - x0) * 0.45);
      const x = ix * dx, z = iy / (NY - 1);
      const a = 14 * ctl.amp;
      const arg = k * x - ph;
      if (mode === 'P') return [a * Math.sin(arg), 0];
      if (mode === 'SV') return [0, -a * Math.sin(arg)];
      if (mode === 'Love') {
        const dec = Math.exp(-1.1 * z);
        return [-a * 0.55 * dec * Math.sin(arg), -a * 0.85 * dec * Math.sin(arg)];
      }
      const decW = Math.exp(-1.25 * z);
      const decU = (0.68 * Math.exp(-1.9 * z) - 0.35 * Math.exp(-0.8 * z)) / 0.33;
      return [-a * 0.72 * decU * Math.cos(arg), -a * 1.05 * decW * Math.sin(arg)];
    }

    svg.addEventListener('click', e => {
      const pt = svg.createSVGPoint();
      pt.x = e.clientX; pt.y = e.clientY;
      const p = pt.matrixTransform(svg.getScreenCTM().inverse());
      let bestD = 1e9;
      for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
        const [ux, uy] = dispAt(i, j, phase);
        const px = x0 + i * dx + ux, py = y0 + j * dy + uy;
        const d = Math.hypot(p.x - px, p.y - py);
        if (d < bestD) { bestD = d; selCol = i; selRow = j; }
      }
      render();
    });

    function render() {
      clearG(g);
      svgText(g, x0, 24, '波傳播方向 → (x)', 'txt-muted');
      svgText(g, W - 45, 24, mode === 'Love' ? '斜向表示水平橫向 (y)' : '深度 z ↓', 'txt-muted', 'end');
      let trD = '';
      for (let s = 0; s <= 48; s++) {
        const [ux, uy] = dispAt(selCol, selRow, s * TAU / 48);
        const px = x0 + selCol * dx + ux, py = y0 + selRow * dy + uy;
        trD += (s ? 'L' : 'M') + px.toFixed(1) + ' ' + py.toFixed(1);
      }
      svgEl('path', { d: trD, fill: 'none', class: 's-blue', 'stroke-width': 1.5, 'stroke-dasharray': '3 2' }, g);

      for (let j = 0; j < NY; j++) {
        let d = '';
        for (let i = 0; i < NX; i++) {
          const [ux, uy] = dispAt(i, j, phase);
          d += (i ? 'L' : 'M') + (x0 + i * dx + ux).toFixed(1) + ' ' + (y0 + j * dy + uy).toFixed(1);
        }
        svgEl('path', { d, fill: 'none', class: 's-gray', 'stroke-width': j === 0 ? 1.5 : 0.8 }, g);
      }
      for (let i = 0; i < NX; i++) {
        let d = '';
        for (let j = 0; j < NY; j++) {
          const [ux, uy] = dispAt(i, j, phase);
          d += (j ? 'L' : 'M') + (x0 + i * dx + ux).toFixed(1) + ' ' + (y0 + j * dy + uy).toFixed(1);
        }
        svgEl('path', { d, fill: 'none', class: 's-gray', 'stroke-width': 0.8 }, g);
      }

      const lam = (x1 - x0) * 0.45;
      let xCrest = x0 + ((phase / TAU + 0.25) % 1) * lam;
      while (xCrest < x0) xCrest += lam;
      while (xCrest <= x1) {
        svgEl('line', { x1: xCrest, x2: xCrest, y1: y0 - 14, y2: y1 + 10, class: 's-red', 'stroke-width': 1.2, 'stroke-dasharray': '4 3' }, g);
        xCrest += lam;
      }

      for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
        const [ux, uy] = dispAt(i, j, phase);
        const px = x0 + i * dx + ux, py = y0 + j * dy + uy;
        const isSel = i === selCol && j === selRow;
        if (isSel) {
          svgEl('circle', { cx: x0 + i * dx, cy: y0 + j * dy, r: 2.5, class: 'f-gray' }, g);
          svgEl('rect', { x: px - 10, y: py - 10, width: 20, height: 20, rx: 4, fill: 'none', class: 's-blue', 'stroke-width': 2 }, g);
        }
        svgEl('circle', { cx: px, cy: py, r: isSel ? 4.8 : 3.2, class: isSel ? 'f-blue' : 'f-ink' }, g);
      }

      const info = {
        P: { name: 'P 波（壓縮波）', dir: '平行波傳方向 (x)', path: '水平直線往復', vol: '有（壓縮/膨脹）', sh: '亦有局部剪應變', decay: '全深度均勻（平面波）' },
        SV: { name: 'SV 波（垂直剪力波）', dir: '垂直於波傳方向 (z)', path: '垂直直線往復', vol: '0（純剪變形）', sh: '有（γ_xz ≠ 0）', decay: '全深度均勻（平面波）' },
        Love: { name: 'SH / Love 波', dir: '水平橫向 (y 軸)', path: '水平橫向直線', vol: '0（純水平剪力）', sh: '有（γ_xy, γ_yz）', decay: 'Love 波隨深度衰減' },
        Rayleigh: { name: 'Rayleigh 波', dir: 'x–z 垂直面內合成', path: '地表逆行橢圓 → 深處順行', vol: '有（近地表）', sh: '有', decay: '隨深度約 e^(−z/λ) 衰減' }
      }[mode];
      st(info);
    }

    let last = null;
    function loop(ts) {
      if (!playing) { last = null; return; }
      if (last != null) phase = (phase + (ts - last) * 0.0032 * ctl.spd) % TAU;
      last = ts;
      render();
      requestAnimationFrame(loop);
    }
    render();
    requestAnimationFrame(loop);
  })();

  /* =========================================================
     2. #w-modulus：Poisson 比如何拉開 Vp 與 Vs (3.2)
     ========================================================= */
  (() => {
    const body = SD.widget('w-modulus'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'nu', label: '骨架 Poisson 比 ν', min: 0.05, max: 0.495, step: 0.005, value: 0.30, fmt: v => v.toFixed(3) },
      { k: 'vs', label: '剪力波速 Vₛ (m/s)', min: 80, max: 500, step: 10, value: 180, fmt: v => v.toFixed(0) },
      { k: 'rho', label: '質量密度 ρ (t/m³)', min: 1.5, max: 2.3, step: 0.05, value: 1.90, fmt: v => v.toFixed(2) },
      { k: 'sat', label: '飽和不排水（加孔隙水 Bw/n ≈ 5.5 GPa）', type: 'check', value: false }
    ], draw);
    SD.buttons(body, '常見情況', [
      { label: '乾砂 ν = 0.25', on: () => ctl.set({ nu: 0.25, sat: false }) },
      { label: '一般土層 ν = 0.33', on: () => ctl.set({ nu: 0.33, sat: false }) },
      { label: '接近不排水 ν = 0.49', on: () => ctl.set({ nu: 0.49, sat: false }) },
      { label: '飽和不排水（Biot 孔隙水）', on: () => ctl.set({ nu: 0.30, vs: 180, sat: true }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pMod = new Plot(grid, { h: 260, xlim: [0, 0.5], ylim: [0, 25], xlabel: 'Poisson 比 ν', ylabel: '模數比（相對於 G）' });
    pMod.legend([{ label: 'M / G（受限模數）', c: 'red', w: 2 }, { label: 'B / G（體積模數）', c: 'blue', w: 2 }, { label: 'E / G = 2(1+ν)', c: 'ink', w: 1.8 }]);
    const pVel = new Plot(grid, { h: 260, xlim: [0, 0.5], ylim: [0, 10], xlabel: 'Poisson 比 ν', ylabel: '波速比 V_p / V_s' });
    pVel.legend([{ label: '骨架 V_p / V_s', c: 'blue', w: 2 }, { label: '飽和不排水 V_p / V_s', c: 'red', dash: true }]);
    const st = SD.stats(body, [
      { k: 'EG', label: 'E / G' }, { k: 'BG', label: 'B / G' }, { k: 'MG', label: 'M / G（骨架）' },
      { k: 'ratio', label: 'V_p / V_s', hl: true }, { k: 'vs', label: 'V_s (m/s)' }, { k: 'vp', label: 'V_p (m/s)', hl: true }
    ]);

    function draw() {
      const nu = ctl.nu, vs = ctl.vs, rho = ctl.rho;
      const G = rho * vs * vs / 1000;
      const EG = 2 * (1 + nu);
      const BG = 2 * (1 + nu) / (3 * (1 - 2 * nu));
      const MG = 2 * (1 - nu) / (1 - 2 * nu);
      const Bw_n = 5500;
      const MG_eff = ctl.sat ? MG + Bw_n / G : MG;
      const vp_vs = Math.sqrt(MG_eff);
      const vp = vs * vp_vs;

      pMod.clear();
      pMod.fn(v => 2 * (1 + v), { c: 'ink', w: 1.8, range: [0, 0.495] });
      pMod.fn(v => Math.min(30, 2 * (1 + v) / (3 * (1 - 2 * v))), { c: 'blue', w: 2, range: [0, 0.495], n: 300 });
      pMod.fn(v => Math.min(30, 2 * (1 - v) / (1 - 2 * v)), { c: 'red', w: 2, range: [0, 0.495], n: 300 });
      pMod.vline(nu, { c: 'gray', dash: '3 3' });
      pMod.markers([nu, nu, nu], [EG, Math.min(24.5, BG), Math.min(24.5, MG)], { c: 'red', fill: true, over: true });

      pVel.clear();
      pVel.fn(v => Math.sqrt(2 * (1 - v) / (1 - 2 * v)), { c: 'blue', w: 2, range: [0, 0.495], n: 300 });
      if (ctl.sat) {
        pVel.hline(vp_vs, { c: 'red', dash: '5 4', w: 1.8 });
        pVel.text(0.03, Math.min(9.3, vp_vs + 0.45), `飽和 V_p ≈ ${vp.toFixed(0)} m/s (V_p/V_s = ${vp_vs.toFixed(1)})`, { cls: 'txt-red' });
      }
      pVel.vline(nu, { c: 'gray', dash: '3 3' });
      pVel.markers([nu], [Math.min(9.8, vp_vs)], { c: ctl.sat ? 'red' : 'blue', fill: true, over: true });

      st({
        EG: f(EG, 2), BG: f(BG, 2), MG: f(MG, 2),
        ratio: f(vp_vs, 2), vs: f(vs, 0), vp: f(vp, 0) + (ctl.sat ? ' (水控制)' : '')
      });
    }
    draw();
  })();

  /* =========================================================
     3. #w-dalembert：一維行進波、駐波與固定端反射 (3.3)
     ========================================================= */
  (() => {
    const body = SD.widget('w-dalembert'); if (!body) return;
    let playing = false, wt = 0;
    const ctl = SD.controls(body, [
      { k: 'r', label: '反射振幅比 R', min: 0, max: 1, step: 0.05, value: 0, fmt: v => v.toFixed(2) },
      { k: 'al', label: '材料阻尼衰減 αλ', min: 0, max: 0.8, step: 0.02, value: 0, fmt: v => v.toFixed(2) },
      { k: 'x0', label: '定點觀測位置 x₀ / λ', min: 0.2, max: 2.6, step: 0.05, value: 1.0, fmt: v => v.toFixed(2) },
      { k: 'comp', label: '顯示右行與左行分量', type: 'check', value: true }
    ], draw);
    SD.buttons(body, '情境與步進', [
      { label: '純右行波 (R = 0)', on: () => ctl.set({ r: 0, al: 0 }) },
      { label: '無阻尼駐波 (R = 1, α = 0)', on: () => ctl.set({ r: 1, al: 0 }) },
      { label: '有阻尼反射 (R = 1, αλ = 0.3)', on: () => ctl.set({ r: 1, al: 0.3 }) },
      { label: '前進 T/8', primary: true, on: () => { wt = (wt + TAU / 8) % (3 * TAU); draw(); } },
      { label: '▶ 播放', on: b => { playing = !playing; b.textContent = playing ? '❚❚ 暫停' : '▶ 播放'; if (playing) requestAnimationFrame(loop); } }
    ]);
    const p = new Plot(body, { h: 280, xlim: [0, 3], ylim: [-2.2, 2.2], xlabel: '空間位置 x / λ（右端 x/λ = 3 為固定端反射面）', ylabel: '正規化位移 u / A' });
    p.legend([
      { label: '合成波 u(x, t)', c: 'ink', w: 2.4 },
      { label: '右行波 f(x − Vt)', c: 'blue', dash: true },
      { label: '左行反射波 g(x + Vt)', c: 'red', dash: '2 3' },
      { label: '紅圓點：跟著相位以 V 右行', c: 'red', line: false, marker: 'circle', fill: true },
      { label: '藍方點：留在 x₀ 原地簡諧', c: 'blue', line: false, marker: 'square', fill: true }
    ]);
    const st = SD.stats(body, [
      { k: 'tT', label: '目前時刻 t / T', hl: true },
      { k: 'type', label: '波場型態', hl: true },
      { k: 'xRed', label: '紅點相位位置 x / λ' },
      { k: 'uGreen', label: '定點位移 u(x₀, t) / A' },
      { k: 'nodeMin', label: '節點最小振幅' }
    ]);

    function uFwd(xl, ph, al) { return Math.exp(-al * xl) * Math.cos(TAU * xl - ph); }
    function uBwd(xl, ph, r, al) { return -r * Math.exp(-al * (6 - xl)) * Math.cos(TAU * (6 - xl) - ph); }

    function draw() {
      const r = ctl.r, al = ctl.al, xg = ctl.x0;
      p.clear();
      p.hline(0, { c: 'gray', dash: null, w: 0.6 });
      p.vline(xg, { c: 'blue', dash: '3 3', w: 1 });
      if (ctl.comp && (r > 0 || al > 0)) {
        p.fn(x => uFwd(x, wt, al), { c: 'blue', dash: '5 4', w: 1.4, n: 350 });
        if (r > 0) p.fn(x => uBwd(x, wt, r, al), { c: 'red', dash: '2 3', w: 1.4, n: 350 });
      }
      if (r > 0) {
        const env = x => {
          const a1 = Math.exp(-al * x), a2 = r * Math.exp(-al * (6 - x));
          return Math.sqrt(a1 * a1 + a2 * a2 - 2 * a1 * a2 * Math.cos(2 * TAU * (3 - x)));
        };
        p.fn(x => env(x), { c: 'gray', dash: '4 4', w: 1, n: 350 });
        p.fn(x => -env(x), { c: 'gray', dash: '4 4', w: 1, n: 350 });
      }
      p.fn(x => uFwd(x, wt, al) + uBwd(x, wt, r, al), { c: 'ink', w: 2.4, n: 400 });

      const xRed = (wt / TAU) % 3;
      const yRed = uFwd(xRed, wt, al) + uBwd(xRed, wt, r, al);
      const yGreen = uFwd(xg, wt, al) + uBwd(xg, wt, r, al);
      p.markers([xRed], [yRed], { c: 'red', fill: true, r: 5, over: true });
      p.markers([xg], [yGreen], { c: 'blue', shape: 'square', fill: true, r: 5, over: true });

      const a1m = Math.exp(-al * 1.5), a2m = r * Math.exp(-al * 4.5);
      const nodeMin = r === 0 ? a1m : Math.abs(a1m - a2m);
      const type = r === 0 ? (al === 0 ? '純右行波' : '衰減行進波') : (r === 1 && al === 0 ? '純駐波（節點完全不動）' : '部分駐波（節點不歸零）');
      st({
        tT: f((wt / TAU) % 3, 3),
        type,
        xRed: f(xRed, 2),
        uGreen: f(yGreen, 3),
        nodeMin: f(nodeMin, 3)
      });
    }
    let last = null;
    function loop(ts) {
      if (!playing) { last = null; return; }
      if (last != null) wt = (wt + (ts - last) * 0.0025) % (3 * TAU);
      last = ts;
      draw();
      requestAnimationFrame(loop);
    }
    draw();
  })();

  /* =========================================================
     4. #w-strain：由 PGV 與 Vs 估剪應變 (3.4)
     ========================================================= */
  (() => {
    const body = SD.widget('w-wavestrain'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'pgv', label: '峰值質點速度 PGV (cm/s)', min: 0.05, max: 100, step: 0.05, value: 25, fmt: v => v < 1 ? v.toFixed(2) : v.toFixed(1) },
      { k: 'vs', label: '剪力波速 Vₛ (m/s)', min: 80, max: 800, step: 10, value: 180, fmt: v => v.toFixed(0) },
      { k: 'freq', label: '主要頻率 f (Hz)', min: 0.5, max: 15, step: 0.5, value: 2.0, fmt: v => v.toFixed(1) }
    ], draw);
    SD.buttons(body, '典型地震與試驗', [
      { label: '微振 / 波速試驗 (PGV = 0.05 cm/s)', on: () => ctl.set({ pgv: 0.05, vs: 200, freq: 10 }) },
      { label: '中小地震 (PGV = 5 cm/s, Vs = 250)', on: () => ctl.set({ pgv: 5, vs: 250, freq: 3 }) },
      { label: '強震軟土 (PGV = 50 cm/s, Vs = 120)', on: () => ctl.set({ pgv: 50, vs: 120, freq: 1.5 }) },
      { label: '強震堅硬土 (PGV = 50 cm/s, Vs = 600)', on: () => ctl.set({ pgv: 50, vs: 600, freq: 2 }) }
    ]);
    const p = new Plot(body, { h: 250, xlim: [-6, -1], ylim: [0, 1.15], xlabel: '剪應變對數 log₁₀ γ（γ = PGV / V_s）', ylabel: '參考割線模數比 G / G_{max}' });
    p.legend([
      { label: '典型砂土 G/G_{max} (γ_r = 0.05%)', c: 'ink', w: 2 },
      { label: '目前 γ_{max} = PGV / V_s', c: 'red', line: false, marker: 'circle', fill: true }
    ]);
    const st = SD.stats(body, [
      { k: 'gam', label: '最大剪應變 γ_{max}', hl: true },
      { k: 'gamp', label: 'γ_{max} (%)', hl: true },
      { k: 'zone', label: '土壤行為區間', hl: true },
      { k: 'ggmax', label: '預估 G / G_{max}' },
      { k: 'umax', label: '估計位移幅 u₀ (mm)' },
      { k: 'amax', label: '估計加速度幅 (g)' }
    ]);

    function draw() {
      const v = ctl.pgv / 100;
      const vs = ctl.vs, fq = ctl.freq, w = TAU * fq;
      const gam = v / vs;
      const lg = Math.max(-6, Math.min(-1, Math.log10(gam)));
      const gr = 5e-4;
      const gg = 1 / (1 + gam / gr);

      p.clear();
      p.band(-5, -4);
      p.band(-3, -1);
      p.vline(-5, { c: 'gray', dash: '3 3' });
      p.vline(-4, { c: 'gray', dash: '3 3' });
      p.vline(-3, { c: 'gray', dash: '3 3' });
      p.text(-5.5, 1.06, '線性彈性 (<10⁻⁵)', { anchor: 'middle', cls: 'txt-muted' });
      p.text(-4.5, 1.06, '微非線性', { anchor: 'middle', cls: 'txt-muted' });
      p.text(-3.5, 1.06, '中度非線性', { anchor: 'middle', cls: 'txt-muted' });
      p.text(-2.0, 1.06, '強非線性 / 破壞 (>10⁻³)', { anchor: 'middle', cls: 'txt-red' });

      p.fn(x => 1 / (1 + Math.pow(10, x) / gr), { c: 'ink', w: 2, n: 300 });
      p.vline(lg, { c: 'red', dash: '4 3', w: 1.5 });
      p.markers([lg], [gg], { c: 'red', fill: true, r: 5, over: true });

      const zone = gam < 1e-5 ? '線性彈性 (G ≈ G_max)' : gam < 1e-4 ? '微非線性（剛開始折減）' : gam < 1e-3 ? '中度非線性（等效線性適用）' : '強非線性（需時域非線性）';
      st({
        gam: gam.toExponential(2),
        gamp: f(gam * 100, 4) + ' %',
        zone,
        ggmax: f(gg, 2),
        umax: f((v / w) * 1000, 1),
        amax: f((v * w) / 9.81, 2)
      });
    }
    draw();
  })();

  /* =========================================================
     5. #w-lambda：波長比較器 (3.5)
     ========================================================= */
  (() => {
    const body = SD.widget('w-lambda'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'vs', label: '剪力波速 V_s (m/s)', min: 80, max: 1000, step: 10, value: 200, fmt: v => v.toFixed(0) },
      { k: 'lf', label: '頻率 f (Hz，對數滑桿)', min: -0.7, max: 3.7, step: 0.02, value: Math.log10(5), fmt: v => { const x = Math.pow(10, v); return x < 10 ? x.toFixed(2) : x.toFixed(0); } }
    ], draw);
    SD.buttons(body, '典型頻段', [
      { label: '強震地動 (f = 2 Hz, Vs = 200)', on: () => ctl.set({ vs: 200, lf: Math.log10(2) }) },
      { label: '機械基礎振動 (f = 25 Hz, Vs = 200)', on: () => ctl.set({ vs: 200, lf: Math.log10(25) }) },
      { label: '共振柱試驗 (f = 108 Hz, Vs = 220)', on: () => ctl.set({ vs: 220, lf: Math.log10(108) }) },
      { label: '彎曲元件 (f = 5000 Hz, Vs = 200)', on: () => ctl.set({ vs: 200, lf: Math.log10(5000) }) }
    ]);
    const p = new Plot(body, { h: 240, xlim: [-2, 3.2], ylim: [0, 1], yticks: [0, 0.5, 1], xlabel: '尺度對數 log₁₀ L (m)（0.01 m ~ 1500 m）', ylabel: '相對位置' });
    const st = SD.stats(body, [
      { k: 'lam', label: '波長 λ = V_s / f (m)', hl: true },
      { k: 'k', label: '波數 k = 2π/λ (rad/m)' },
      { k: 'rc', label: '共振柱試體 (0.1 m)' },
      { k: 'fdn', label: '淺基礎 (2 m)' },
      { k: 'soil', label: '30 m 土層' },
      { k: 'bldg', label: '10 層樓 (35 m)' }
    ]);

    const refs = [
      { name: '共振柱試體 (0.1 m)', L: 0.1, y: 0.22 },
      { name: '淺基礎 (2 m)', L: 2.0, y: 0.42 },
      { name: '30 m 土層', L: 30.0, y: 0.64 },
      { name: '10 層樓 (35 m)', L: 35.0, y: 0.84 }
    ];

    function judge(L, lam) {
      const r = L / lam;
      if (r < 0.05) return `L/λ = ${f(r, 3)}（準靜態/均勻）`;
      if (r < 0.25) return `L/λ = ${f(r, 2)}（過渡區）`;
      return `L/λ = ${f(r, 2)}（波動/共振主導）`;
    }

    function draw() {
      const vs = ctl.vs, fq = Math.pow(10, ctl.lf);
      const lam = vs / fq, lk = TAU / lam;
      const ll = Math.max(-2, Math.min(3.2, Math.log10(lam)));

      p.clear();
      for (const r of refs) {
        const lr = Math.log10(r.L);
        p.vline(lr, { c: 'gray', dash: '3 3', w: 1 });
        p.markers([lr], [r.y], { c: 'ink', shape: 'square', fill: true, r: 4, over: true });
        p.text(lr + 0.06, r.y + 0.03, r.name, { cls: 'txt' });
      }
      p.band(ll - 0.05, ll + 0.05);
      p.vline(ll, { c: 'red', dash: null, w: 2.4 });
      p.text(Math.min(2.5, Math.max(-1.4, ll)), 0.08, `目前波長 λ = ${lam < 1 ? lam.toFixed(3) : lam.toFixed(1)} m`, { cls: 'txt-red', anchor: 'middle' });

      st({
        lam: lam < 1 ? f(lam, 3) : f(lam, 2),
        k: lk > 10 ? f(lk, 1) : f(lk, 3),
        rc: judge(0.1, lam),
        fdn: judge(2.0, lam),
        soil: judge(30.0, lam),
        bldg: judge(35.0, lam)
      });
    }
    draw();
  })();

  /* =========================================================
     6. #w-rayleigh：Rayleigh 波振幅隨深度與波速比 (3.6)
     ========================================================= */
  (() => {
    const body = SD.widget('w-rayleigh'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'nu', label: 'Poisson 比 ν', min: 0.0, max: 0.495, step: 0.005, value: 0.25, fmt: v => v.toFixed(3) },
      { k: 'zmax', label: '顯示深度上限 z / λR', min: 1.0, max: 2.5, step: 0.25, value: 2.0, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '典型 ν', [
      { label: 'ν = 0 (K = 0.874)', on: () => ctl.set({ nu: 0 }) },
      { label: 'ν = 0.25 (K = 0.9194)', on: () => ctl.set({ nu: 0.25 }) },
      { label: 'ν = 0.33 (K = 0.9325)', on: () => ctl.set({ nu: 0.33 }) },
      { label: 'ν = 0.49 (K = 0.954)', on: () => ctl.set({ nu: 0.49 }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pAmp = new Plot(grid, { h: 290, xlim: [-0.4, 1.25], ylim: [-2.0, 0], yfmt: v => Math.abs(v).toFixed(1), xlabel: '正規化位移（以地表垂直幅 W₀ = 1 為準）', ylabel: '正規化深度 z / λ_R' });
    pAmp.legend([
      { label: '垂直位移 W(z)/W₀', c: 'blue', w: 2.2 },
      { label: '水平位移 U(z)/W₀（正=逆行，負=順行）', c: 'red', w: 2 },
      { label: '|U(z)|/W₀（絕對值）', c: 'red', dash: true }
    ]);
    const pK = new Plot(grid, { h: 290, xlim: [0, 0.5], ylim: [0, 1.05], xlabel: 'Poisson 比 ν', ylabel: '波速比' });
    pK.legend([
      { label: 'K = V_R / V_s', c: 'blue', w: 2.2 },
      { label: 'V_R / V_p', c: 'red', w: 2 },
      { label: '近似 K ≈ (0.87+1.12ν)/(1+ν)', c: 'gray', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'K', label: 'K = V_R / V_s', hl: true },
      { k: 'Kvp', label: 'V_R / V_p' },
      { k: 'U0W0', label: '地表軸比 |U₀ / W₀|' },
      { k: 'z0', label: '水平過零深度 z₀ / λ_R', hl: true },
      { k: 'w1', label: 'z = λ_R 處 W / W₀' },
      { k: 'kap', label: '經驗公式誤差 (%)' }
    ]);

    function rwProfile(nu, zl) {
      const K = rayleighK(nu);
      const a2 = (1 - 2 * nu) / (2 * (1 - nu));
      const q = Math.sqrt(1 - a2 * K * K);
      const s = Math.sqrt(1 - K * K);
      const b = 1 - 0.5 * K * K;
      const kz = TAU * zl;
      const eq = Math.exp(-q * kz), es = Math.exp(-s * kz);
      const Uraw = eq - b * es;
      const Wraw = -q * eq + (q / b) * es;
      const W0 = -q + q / b;
      return { U: Uraw / W0, W: Wraw / W0, K, q, s, b, W0 };
    }

    function draw() {
      const nu = ctl.nu, zmax = ctl.zmax;
      const p0 = rwProfile(nu, 0);
      const K = p0.K;
      const a = Math.sqrt((1 - 2 * nu) / (2 * (1 - nu)));
      const z0 = Math.log(1 / p0.b) / (TAU * (p0.q - p0.s));

      pAmp.setLim([-0.4, 1.25], [-zmax, 0]);
      pAmp.clear();
      pAmp.vline(0, { c: 'gray', dash: null, w: 0.7 });
      pAmp.hline(-z0, { c: 'gray', dash: '4 3', w: 1 });
      const zs = linspace(0, zmax, 220);
      const nzs = zs.map(z => -z);
      const Us = zs.map(z => rwProfile(nu, z).U);
      const Ws = zs.map(z => rwProfile(nu, z).W);
      pAmp.line(Ws, nzs, { c: 'blue', w: 2.2 });
      pAmp.line(Us, nzs, { c: 'red', w: 2 });
      pAmp.line(Us.map(Math.abs), nzs, { c: 'red', dash: '4 3', w: 1.3 });
      pAmp.markers([0], [-z0], { c: 'red', fill: true, r: 4.5, over: true });
      pAmp.text(0.18, -Math.min(zmax - 0.1, z0 + 0.11), `過零深度 z₀/λ_R = ${z0.toFixed(3)}`, { cls: 'txt-red' });

      pK.clear();
      pK.fn(v => rayleighK(v), { c: 'blue', w: 2.2, range: [0, 0.499], n: 200 });
      pK.fn(v => (0.87 + 1.12 * v) / (1 + v), { c: 'gray', dash: '4 3', w: 1.5, range: [0, 0.499] });
      pK.fn(v => rayleighK(v) * Math.sqrt((1 - 2 * v) / (2 * (1 - v))), { c: 'red', w: 2, range: [0, 0.499], n: 200 });
      pK.vline(nu, { c: 'gray', dash: '3 3' });
      pK.markers([nu, nu], [K, K * a], { c: 'blue', fill: true, over: true });

      const kApprox = (0.87 + 1.12 * nu) / (1 + nu);
      const w1 = rwProfile(nu, 1.0).W;
      st({
        K: f(K, 4),
        Kvp: f(K * a, 3),
        U0W0: f(Math.abs(p0.U), 3),
        z0: f(z0, 3),
        w1: f(w1, 3),
        kap: f(Math.abs(kApprox - K) / K * 100, 2) + ' %'
      });
    }
    draw();
  })();

  /* =========================================================
     7. #w-love：Love 波頻散曲線與振幅隨深度 (3.7)
     ========================================================= */
  (() => {
    const body = SD.widget('w-love'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'h', label: '上層厚度 h (m)', min: 4, max: 30, step: 1, value: 10, fmt: v => v.toFixed(0) },
      { k: 'v1', label: '上層波速 V₁ (m/s)', min: 80, max: 300, step: 5, value: 150, fmt: v => v.toFixed(0) },
      { k: 'v2', label: '半空間波速 V₂ (m/s)', min: 160, max: 600, step: 10, value: 400, fmt: v => v.toFixed(0) },
      { k: 'freq', label: '檢視頻率 f (Hz)', min: 0.5, max: 25, step: 0.5, value: 5.0, fmt: v => v.toFixed(1) },
      { k: 'grp', label: '顯示群速度 U(f)', type: 'check', value: true },
      { k: 'hi', label: '顯示第 1、2 高階模態', type: 'check', value: true }
    ], draw);
    SD.buttons(body, '講義比較範例', [
      { label: '講義基準 (h=10, V₁=150, V₂=400)', on: () => ctl.set({ h: 10, v1: 150, v2: 400, freq: 5 }) },
      { label: '波速對比變小 (V₂ = 170 m/s)', on: () => ctl.set({ h: 10, v1: 150, v2: 170, freq: 5 }) },
      { label: '軟層加厚一倍 (h = 20 m)', on: () => ctl.set({ h: 20, v1: 150, v2: 400, freq: 5 }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pDisp = new Plot(grid, { h: 290, xlim: [0, 25], ylim: [100, 450], xlabel: '頻率 f (Hz)', ylabel: '速度 V (m/s)' });
    pDisp.legend([
      { label: '基本模態 V_L (n=0)', c: 'blue', w: 2.2 },
      { label: '群速度 U (n=0)', c: 'blue', dash: true },
      { label: '高階模態 (n=1, 2)', c: 'red', w: 1.7 },
      { label: '目前頻率解', c: 'red', line: false, marker: 'circle', fill: true }
    ]);
    const pProf = new Plot(grid, { h: 290, xlim: [-1.1, 1.1], ylim: [-35, 0], yfmt: v => String(Math.round(Math.abs(v))), xlabel: '正規化振幅 A(z) / A₀', ylabel: '深度 z (m)' });
    pProf.legend([
      { label: '基本模態振幅分布', c: 'blue', w: 2.2 },
      { label: '層界面 z = h', c: 'gray', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'vl', label: '相速度 V_L (m/s)', hl: true },
      { k: 'ul', label: '群速度 U (m/s)', hl: true },
      { k: 'lam', label: '波長 λ_L = V_L/f (m)' },
      { k: 'ah', label: '界面振幅比 cos(p₁h)' },
      { k: 'fc1', label: 'n=1 截止頻率 f_{c1} (Hz)' },
      { k: 'pen', label: '半空間衰減深度 1/p₂ (m)' }
    ]);

    const rho1 = 1.8, rho2 = 2.0;

    function loveFreqOfV(v, n, h, v1, v2, gRatio) {
      const s1 = Math.sqrt(1 / (v1 * v1) - 1 / (v * v));
      const s2 = Math.sqrt(1 / (v * v) - 1 / (v2 * v2));
      const phi = Math.atan(gRatio * s2 / s1) + n * Math.PI;
      return phi / (TAU * h * s1);
    }

    function solveLoveAtF(fq, n, h, v1, v2, gRatio) {
      if (v2 <= v1) return NaN;
      const fc = n / (2 * h * Math.sqrt(1 / (v1 * v1) - 1 / (v2 * v2)));
      if (fq <= fc) return NaN;
      let lo = v1 * (1 + 1e-7), hi = v2 * (1 - 1e-7);
      for (let i = 0; i < 55; i++) {
        const m = 0.5 * (lo + hi);
        if (loveFreqOfV(m, n, h, v1, v2, gRatio) > fq) lo = m; else hi = m;
      }
      return 0.5 * (lo + hi);
    }

    function draw() {
      const h = ctl.h, v1 = ctl.v1, v2 = Math.max(v1 + 10, ctl.v2), fq = ctl.freq;
      const g1 = rho1 * v1 * v1, g2 = rho2 * v2 * v2, gRatio = g2 / g1;

      pDisp.setLim([0, 25], [Math.floor(v1 * 0.8 / 20) * 20, Math.ceil(v2 * 1.08 / 20) * 20]);
      pDisp.clear();
      pDisp.hline(v1, { c: 'gray', dash: '3 3', w: 1 });
      pDisp.hline(v2, { c: 'gray', dash: '3 3', w: 1 });

      const vs = linspace(v1 * 1.0005, v2 * 0.9995, 180);
      const maxM = ctl.hi ? 2 : 0;
      for (let n = 0; n <= maxM; n++) {
        const fs = vs.map(v => loveFreqOfV(v, n, h, v1, v2, gRatio));
        pDisp.line(fs, vs, { c: n === 0 ? 'blue' : 'red', w: n === 0 ? 2.2 : 1.6 });
        if (ctl.grp) {
          const us = [], fgu = [];
          for (let i = 1; i < vs.length - 1; i++) {
            const f0 = fs[i - 1], f1 = fs[i + 1], v0 = vs[i - 1], v_1 = vs[i + 1];
            const k0 = TAU * f0 / v0, k1 = TAU * f1 / v_1;
            const U = TAU * (f1 - f0) / (k1 - k0);
            if (f0 <= 25 && isFinite(U)) { fgu.push(fs[i]); us.push(U); }
          }
          pDisp.line(fgu, us, { c: n === 0 ? 'blue' : 'red', dash: '4 3', w: 1.4 });
        }
      }

      const vl = solveLoveAtF(fq, 0, h, v1, v2, gRatio);
      const vl_p = solveLoveAtF(fq + 0.05, 0, h, v1, v2, gRatio);
      const vl_m = solveLoveAtF(Math.max(0.1, fq - 0.05), 0, h, v1, v2, gRatio);
      const k_p = TAU * (fq + 0.05) / vl_p, k_m = TAU * Math.max(0.1, fq - 0.05) / vl_m;
      const ul = TAU * (fq + 0.05 - Math.max(0.1, fq - 0.05)) / (k_p - k_m);

      pDisp.vline(fq, { c: 'gray', dash: '3 3' });
      if (isFinite(vl)) pDisp.markers([fq], [vl], { c: 'red', fill: true, r: 5, over: true });

      const zmax = Math.max(25, 2.5 * h);
      pProf.setLim([-1.1, 1.1], [-zmax, 0]);
      pProf.clear();
      pProf.vline(0, { c: 'gray', dash: null, w: 0.7 });
      pProf.hline(-h, { c: 'gray', dash: '5 4', w: 1.2 });
      pProf.text(0.95, -h + 1.2, `界面 z = ${h} m`, { anchor: 'end', cls: 'txt-muted' });

      const w = TAU * fq;
      const p1 = w * Math.sqrt(1 / (v1 * v1) - 1 / (vl * vl));
      const p2 = w * Math.sqrt(1 / (vl * vl) - 1 / (v2 * v2));
      const zs = linspace(0, zmax, 220);
      const amps = zs.map(z => z <= h ? Math.cos(p1 * z) : Math.cos(p1 * h) * Math.exp(-p2 * (z - h)));
      pProf.line(amps, zs.map(z => -z), { c: 'blue', w: 2.2 });
      pProf.markers([Math.cos(p1 * h)], [-h], { c: 'red', fill: true, r: 4.5, over: true });

      const fc1 = 1 / (2 * h * Math.sqrt(1 / (v1 * v1) - 1 / (v2 * v2)));
      st({
        vl: f(vl, 1),
        ul: f(ul, 1),
        lam: f(vl / fq, 1),
        ah: f(Math.cos(p1 * h), 3),
        fc1: f(fc1, 2),
        pen: f(1 / p2, 2)
      });
    }
    draw();
  })();

  /* =========================================================
     8. #w-masw：MASW 頻散分析實驗室 (3.8 & 3.10 Q1)
     ========================================================= */
  (() => {
    const body = SD.widget('w-masw'); if (!body) return;
    const SITES = {
      normal: {
        name: 'Q1(a) 正常遞增剖面 (150/220/320/450 m/s)',
        layers: [
          { h: 3, vs: 150, nu: 0.35, gamma: 17 },
          { h: 5, vs: 220, nu: 0.33, gamma: 18 },
          { h: 7, vs: 320, nu: 0.30, gamma: 19 },
          { h: 0, vs: 450, nu: 0.30, gamma: 20 }
        ]
      },
      crust: {
        name: 'Q1(b) 硬殼夾軟弱黏土 (300/140/250/450 m/s)',
        layers: [
          { h: 2, vs: 300, nu: 0.30, gamma: 19 },
          { h: 6, vs: 140, nu: 0.45, gamma: 16.5 },
          { h: 8, vs: 250, nu: 0.40, gamma: 18 },
          { h: 0, vs: 450, nu: 0.30, gamma: 20 }
        ]
      },
      sat: {
        name: 'Q1(c) 飽和沖積層 (140/180/260/420 m/s)',
        layers: [
          { h: 2, vs: 140, nu: 0.35, gamma: 17.5 },
          { h: 6, vs: 180, nu: 0.49, gamma: 18.5 },
          { h: 8, vs: 260, nu: 0.49, gamma: 19 },
          { h: 0, vs: 420, nu: 0.45, gamma: 20 }
        ]
      }
    };
    let siteKey = 'normal', forceNu = false;
    const ctl = SD.controls(body, [
      { k: 'N', label: '接收器數 N', min: 8, max: 48, step: 2, value: 24, fmt: v => v.toFixed(0) },
      { k: 'dx', label: '接收器間距 Δx (m)', min: 0.5, max: 4.0, step: 0.25, value: 1.5, fmt: v => v.toFixed(2) },
      { k: 'nu', label: '各層 Poisson 比 ν', min: 0.20, max: 0.49, step: 0.01, value: 0.33, fmt: v => forceNu ? v.toFixed(2) : '預設' },
      { k: 'useNu', label: '以滑桿 ν 覆寫所有土層（Q1c 比較用）', type: 'check', value: false }
    ], () => { forceNu = ctl.useNu; draw(); });

    SD.buttons(body, '範例場址', [
      { label: '正常遞增 (N=24)', on: () => { siteKey = 'normal'; ctl.set({ N: 24, dx: 1.5, useNu: false }); } },
      { label: '正常遞增、減半陣列 (N=12)', on: () => { siteKey = 'normal'; ctl.set({ N: 12, dx: 1.5, useNu: false }); } },
      { label: '硬殼夾軟層 (140 m/s)', on: () => { siteKey = 'crust'; ctl.set({ N: 24, dx: 1.5, useNu: false }); } },
      { label: '飽和場址 (ν = 0.25)', on: () => { siteKey = 'sat'; ctl.set({ N: 24, dx: 1.5, nu: 0.25, useNu: true }); } },
      { label: '飽和場址 (ν = 0.49)', on: () => { siteKey = 'sat'; ctl.set({ N: 24, dx: 1.5, nu: 0.49, useNu: true }); } }
    ]);

    const grid = SD.div(body, 'plot-row plot-grid');
    const pDC = new Plot(grid, { h: 290, xlim: [2, 55], ylim: [80, 460], xlabel: '頻率 f (Hz)', ylabel: 'Rayleigh 相速度 c_R (m/s)' });
    pDC.legend([
      { label: '基本模態 M₀ (可解析段 λ ≤ L)', c: 'blue', w: 2.4 },
      { label: '超出展距限制 (λ > L)', c: 'gray', dash: true },
      { label: '第一高階模態 M₁', c: 'red', dash: '4 3', w: 1.8 },
      { label: '陣列上限 c = L·f', c: 'ink', dash: '2 3' }
    ]);
    const pProf = new Plot(grid, { h: 290, xlim: [80, 480], ylim: [-22, 0], yfmt: v => String(Math.round(Math.abs(v))), xlabel: '剪力波速 V_s (m/s)', ylabel: '深度 z (m)' });
    pProf.legend([
      { label: '真實分層 V_s(z)', c: 'ink', w: 2.2 },
      { label: 'λ/3 偽深度估計 (V_s ≈ 1.1 c_R)', c: 'blue', line: false, marker: 'circle', fill: true },
      { label: '有效探查深度 L/3', c: 'red', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'site', label: '目前場址', hl: true },
      { k: 'L', label: '陣列展距 L = (N−1)Δx (m)', hl: true },
      { k: 'zmax', label: '探查深度 ≈ L/3 ~ L/2 (m)' },
      { k: 'c10', label: 'c_R(10 Hz) (m/s)', hl: true },
      { k: 'c30', label: 'c_R(30 Hz) (m/s)', hl: true },
      { k: 'fmin', label: '最低可解析頻率 (λ = L)' }
    ]);

    function draw() {
      forceNu = ctl.useNu;
      const base = SITES[siteKey].layers.map(l => ({ ...l, nu: forceNu ? ctl.nu : l.nu }));
      const M = KR.deriveModel(base);
      const L = (ctl.N - 1) * ctl.dx;

      const freqs = linspace(4, 52, 49);
      const m0 = [], m1 = [], fValid = [], cValid = [];
      let fminPick = NaN;
      for (const fq of freqs) {
        const r = KR.findModes(fq, M, 2, 160);
        const c0 = r[0] || NaN, c1 = r[1] || NaN;
        m0.push(c0); m1.push(c1);
        if (isFinite(c0) && c0 / fq <= L) {
          fValid.push(fq); cValid.push(c0);
          if (!isFinite(fminPick)) fminPick = fq;
        }
      }
      const c10 = (KR.findModes(10, M, 1, 180)[0]) || NaN;
      const c30 = (KR.findModes(30, M, 1, 180)[0]) || NaN;

      pDC.clear();
      pDC.fn(fq => L * fq, { c: 'ink', dash: '3 3', w: 1.3 });
      pDC.fn(fq => 2 * ctl.dx * fq, { c: 'gray', dash: '2 3', w: 1.2 });
      pDC.line(freqs, m0, { c: 'gray', dash: '4 3', w: 1.5 });
      pDC.line(fValid, cValid, { c: 'blue', w: 2.5 });
      pDC.line(freqs, m1, { c: 'red', dash: '4 3', w: 1.8 });
      pDC.markers([10, 30], [c10, c30], { c: 'blue', fill: true, r: 4.5, over: true });

      pProf.clear();
      const xs = [], nzs = [];
      let z = 0;
      for (let i = 0; i < base.length; i++) {
        const zb = i === base.length - 1 ? 22 : z + base[i].h;
        xs.push(base[i].vs, base[i].vs);
        nzs.push(-z, -zb);
        z = zb;
      }
      pProf.line(xs, nzs, { c: 'ink', w: 2.2 });
      pProf.hline(-L / 3, { c: 'red', dash: '5 4', w: 1.4 });
      pProf.text(465, -Math.min(20.8, L / 3 - 0.8), `z ≈ L/3 = ${(L / 3).toFixed(1)} m`, { anchor: 'end', cls: 'txt-red' });

      const pVs = cValid.map(c => 1.1 * c);
      const pZ = cValid.map((c, i) => -(c / fValid[i]) / 3);
      pProf.markers(pVs, pZ, { c: 'blue', fill: true, r: 3.5, over: true });

      st({
        site: SITES[siteKey].name.split(' ')[0] + ' ' + SITES[siteKey].name.split(' ')[1],
        L: f(L, 1),
        zmax: `${f(L / 3, 1)} – ${f(L / 2, 1)}`,
        c10: f(c10, 1),
        c30: f(c30, 1),
        fmin: isFinite(fminPick) ? `≈ ${f(fminPick, 1)} Hz` : '> 50 Hz'
      });
    }
    draw();
  })();

  /* =========================================================
     9. #w-vs30：Vs30 計算器與三套分類 (3.9 & 3.10 Q2a,b)
     ========================================================= */
  (() => {
    const body = SD.widget('w-vs30'); if (!body) return;
    let layers = [
      { h: 5, vs: 120, desc: '軟弱粉質黏土' },
      { h: 10, vs: 250, desc: '中等緊密砂' },
      { h: 15, vs: 400, desc: '緊密礫石' }
    ];

    function classifyAll(v, clipped) {
      const nehrp = v >= 1500 ? 'A' : v >= 760 ? 'B' : v >= 360 ? 'C' : v >= 180 ? 'D' : 'E';
      const softH = clipped ? clipped.reduce((s, l) => s + (l.vs < 100 ? l.d : 0), 0) : 0;
      const ec8Base = v > 800 ? 'A' : v >= 360 ? 'B' : v >= 180 ? 'C' : 'D';
      const ec8 = softH >= 10 ? `${ec8Base} (且符合 S₁ 特殊軟土)` : ec8Base;
      const tw = v >= 270 ? '第一類（堅實）' : v >= 180 ? '第二類（普通）' : '第三類（軟弱）';
      return { nehrp, ec8, tw };
    }

    SD.buttons(body, '預設剖面', [
      { label: '三層範例剖面 (5/10/15 m)', on: () => { layers = [{ h: 5, vs: 120, desc: '軟弱粉質黏土' }, { h: 10, vs: 250, desc: '中等緊密砂' }, { h: 15, vs: 400, desc: '緊密礫石' }]; rebuild(); } },
      { label: '上課例題 (10 m 軟黏土 / 20 m 風化岩)', on: () => { layers = [{ h: 10, vs: 80, desc: '軟黏土' }, { h: 20, vs: 600, desc: '風化岩' }]; rebuild(); } },
      { label: '堅硬礫石與岩盤', on: () => { layers = [{ h: 6, vs: 320, desc: '緊密砂礫' }, { h: 12, vs: 550, desc: '軟岩' }, { h: 12, vs: 900, desc: '新鮮岩盤' }]; rebuild(); } },
      { label: '+ 新增土層', on: () => { layers.push({ h: 5, vs: 300, desc: '新增土層' }); rebuild(); } }
    ]);

    const tblWrap = SD.div(body, '');
    const grid = SD.div(body, 'plot-row plot-grid');
    const pVs = new Plot(grid, { h: 260, xlim: [0, 700], ylim: [-30, 0], yfmt: v => String(Math.round(Math.abs(v))), xlabel: '剪力波速 V_s (m/s)', ylabel: '深度 z (m)' });
    pVs.legend([
      { label: '分層 V_s(z)', c: 'ink', w: 2.2 },
      { label: '走時平均 V_s30', c: 'blue', w: 2 },
      { label: '算術平均 V̄_s', c: 'red', dash: true }
    ]);
    const pTt = new Plot(grid, { h: 260, xlim: [0, 0.25], ylim: [-30, 0], yfmt: v => String(Math.round(Math.abs(v))), xlabel: '累積走時 t(z) = ∫ dz/V_s (s)', ylabel: '深度 z (m)' });
    pTt.legend([{ label: '累積走時曲線 t(z)', c: 'blue', w: 2.2 }]);

    const clsBox = SD.div(body, 'cls-row');
    const st = SD.stats(body, [
      { k: 'tt', label: '30 m 總走時 Σ(d_i/V_{si}) (s)', hl: true },
      { k: 'vs30', label: '走時平均 V_s30 (m/s)', hl: true },
      { k: 'vbar', label: '算術平均 V̄_s (m/s)' },
      { k: 'ratio', label: '比值 V̄_s / V_s30' },
      { k: 'Tsite', label: '近似基本週期 4×30/V_s30 (s)' }
    ]);

    function clip30(L) {
      const out = [];
      let z = 0;
      for (let i = 0; i < L.length; i++) {
        const rem = Math.max(0, 30 - z);
        const d = i === L.length - 1 ? Math.max(L[i].h, rem) : Math.min(L[i].h, rem);
        const use = Math.min(d, rem);
        out.push({ ...L[i], d: use, z0: z, z1: z + use });
        z += use;
      }
      return out;
    }

    function rebuild() {
      const clipped = clip30(layers);
      const tt = clipped.reduce((s, l) => s + (l.d > 0 ? l.d / l.vs : 0), 0);
      tblWrap.innerHTML = `<table class="layer-tbl"><thead><tr>
        <th>層</th><th>描述</th><th>輸入厚度 h (m)</th><th>V<sub>s</sub> (m/s)</th><th>30 m 內採用 d<sub>i</sub> (m)</th><th>走時 d<sub>i</sub>/V<sub>si</sub> (s)</th><th>走時佔比 (%)</th><th></th>
      </tr></thead><tbody>${clipped.map((l, i) => {
        const ti = l.d > 0 ? l.d / l.vs : 0;
        return `<tr>
          <td>${i + 1}</td>
          <td class="t">${l.desc || '—'}</td>
          <td><input type="number" data-i="${i}" data-k="h" min="0.5" max="50" step="0.5" value="${l.h}"></td>
          <td><input type="number" data-i="${i}" data-k="vs" min="40" max="2000" step="10" value="${l.vs}"></td>
          <td class="c-num">${f(l.d, 1)}</td>
          <td class="c-num">${f(ti, 4)}</td>
          <td class="c-num">${f(tt > 0 ? 100 * ti / tt : 0, 1)}%</td>
          <td>${layers.length > 1 ? `<button class="btn" data-del="${i}" style="padding:1px 7px">×</button>` : ''}</td>
        </tr>`;
      }).join('')}</tbody></table>`;

      tblWrap.querySelectorAll('input').forEach(inp => {
        inp.addEventListener('change', () => {
          const i = +inp.dataset.i, k = inp.dataset.k;
          const v = parseFloat(inp.value);
          if (isFinite(v) && v > 0) { layers[i][k] = v; rebuild(); }
        });
      });
      tblWrap.querySelectorAll('button[data-del]').forEach(btn => {
        btn.addEventListener('click', () => {
          layers.splice(+btn.dataset.del, 1);
          rebuild();
        });
      });
      drawOnly(clipped, tt);
    }

    function drawOnly(clipped, tt) {
      const vs30 = 30 / tt;
      const vbar = clipped.reduce((s, l) => s + l.d * l.vs, 0) / 30;
      const vmax = Math.max(650, ...clipped.map(l => l.vs * 1.15));
      pVs.setLim([0, Math.ceil(vmax / 50) * 50], [-30, 0]);
      pVs.clear();
      const xs = [], nzs = [], ts = [0], nzts = [0];
      let tCum = 0;
      for (const l of clipped) {
        if (l.d <= 0) continue;
        xs.push(l.vs, l.vs);
        nzs.push(-l.z0, -l.z1);
        tCum += l.d / l.vs;
        ts.push(tCum);
        nzts.push(-l.z1);
      }
      pVs.line(xs, nzs, { c: 'ink', w: 2.2 });
      pVs.vline(vs30, { c: 'blue', dash: null, w: 2 });
      pVs.vline(vbar, { c: 'red', dash: '5 4', w: 1.8 });

      pTt.setLim([0, Math.max(0.15, Math.ceil(tt * 20) / 20)], [-30, 0]);
      pTt.clear();
      pTt.line(ts, nzts, { c: 'blue', w: 2.2 });
      pTt.markers(ts, nzts, { c: 'blue', fill: true, r: 4, over: true });

      const c30 = classifyAll(vs30, clipped), cbar = classifyAll(vbar, clipped);
      clsBox.innerHTML = [
        ['NEHRP / ASCE 7', `Site Class ${c30.nehrp}`, `算術平均誤判為 ${cbar.nehrp}`],
        ['Eurocode 8', `Ground Type ${c30.ec8}`, `算術平均為 ${cbar.ec8}`],
        ['臺灣耐震規範 (2022)', c30.tw, `算術平均為 ${cbar.tw}`]
      ].map(([sys, a, b]) => `<div class="cls-chip"><div class="sys">${sys}</div><div class="val">${a}</div><div style="font-size:.75rem;color:var(--muted);margin-top:2px">${b}</div></div>`).join('');

      st({
        tt: f(tt, 4),
        vs30: f(vs30, 1),
        vbar: f(vbar, 1),
        ratio: f(vbar / vs30, 2),
        Tsite: f(120 / vs30, 2)
      });
    }
    rebuild();
  })();

  /* =========================================================
     10. #w-hw2c：HW2 Q2(c) 表層軟土厚度與場址類別 (3.10)
     ========================================================= */
  (() => {
    const body = SD.widget('w-hw2c'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'H', label: '表層軟土厚度 H (m)', min: 0, max: 30, step: 0.25, value: 10, fmt: v => v.toFixed(2) },
      { k: 'v1', label: '軟土波速 V_s1 (m/s)', min: 50, max: 180, step: 5, value: 80, fmt: v => v.toFixed(0) },
      { k: 'v2', label: '下部岩層 V_s2 (m/s)', min: 300, max: 900, step: 20, value: 600, fmt: v => v.toFixed(0) }
    ], draw);
    SD.buttons(body, '關鍵厚度門檻', [
      { label: 'H = 5 m (Vs30 = 289.4)', on: () => ctl.set({ H: 5, v1: 80, v2: 600 }) },
      { label: 'H = 10 m 上課例題 (Vs30 = 189.5)', on: () => ctl.set({ H: 10, v1: 80, v2: 600 }) },
      { label: 'H = 10.91 m (Vs30 跌破 180 門檻)', on: () => ctl.set({ H: 10.91, v1: 80, v2: 600 }) },
      { label: 'H = 24.23 m (算術平均才跌破 180)', on: () => ctl.set({ H: 24.23, v1: 80, v2: 600 }) }
    ]);
    const p = new Plot(body, { h: 280, xlim: [0, 30], ylim: [0, 650], xlabel: '表層軟土厚度 H (m)（總深度 30 m）', ylabel: '平均波速 (m/s)' });
    p.legend([
      { label: '走時平均 V_s30(H)（調和平均）', c: 'blue', w: 2.3 },
      { label: '厚度加權算術平均 V̄_s(H)', c: 'red', dash: true, w: 2 },
      { label: '規範分類門檻 (360 / 270 / 180 m/s)', c: 'gray', dash: '2 3' }
    ]);
    const st = SD.stats(body, [
      { k: 'vs30', label: '走時平均 V_s30 (m/s)', hl: true },
      { k: 'vbar', label: '算術平均 V̄_s (m/s)', hl: true },
      { k: 'tshare', label: '軟土走時佔比 (%)' },
      { k: 'nehrp', label: 'NEHRP (Vs30 vs 算術)' },
      { k: 'tw', label: '臺灣規範 (Vs30 vs 算術)' },
      { k: 'h180', label: '跌破 180 m/s 的軟土厚' }
    ]);

    function draw() {
      const H = ctl.H, v1 = ctl.v1, v2 = ctl.v2;
      const calc30 = h => 30 / (h / v1 + (30 - h) / v2);
      const calcBar = h => (h * v1 + (30 - h) * v2) / 30;
      const vs30 = calc30(H), vbar = calcBar(H);
      const t1 = H / v1, tt = t1 + (30 - H) / v2;

      p.setLim([0, 30], [0, Math.ceil(v2 * 1.08 / 50) * 50]);
      p.clear();
      for (const [lvl, txt] of [[360, '360 m/s (NEHRP C/D, EC8 B/C)'], [270, '270 m/s (臺灣 一/二類)'], [180, '180 m/s (NEHRP D/E, 臺灣 二/三類)']]) {
        p.hline(lvl, { c: 'gray', dash: '3 3', w: 1 });
        p.text(29.5, lvl + 10, txt, { anchor: 'end', cls: 'txt-muted' });
      }
      p.fn(calc30, { c: 'blue', w: 2.3, n: 300 });
      p.fn(calcBar, { c: 'red', dash: '5 4', w: 2, n: 100 });
      p.vline(H, { c: 'ink', dash: '3 3', w: 1 });
      p.markers([H, H], [vs30, vbar], { c: 'blue', fill: true, r: 4.5, over: true });

      const cN = v => v >= 760 ? 'B' : v >= 360 ? 'C' : v >= 180 ? 'D' : 'E';
      const cT = v => v >= 270 ? '第一類' : v >= 180 ? '第二類' : '第三類';
      const h180 = (30 / 180 - 30 / v2) / (1 / v1 - 1 / v2);

      st({
        vs30: f(vs30, 1),
        vbar: f(vbar, 1),
        tshare: f(100 * t1 / tt, 1) + ' %',
        nehrp: `${cN(vs30)} vs ${cN(vbar)}`,
        tw: `${cT(vs30)} vs ${cT(vbar)}`,
        h180: h180 >= 0 && h180 <= 30 ? f(h180, 2) + ' m' : '—'
      });
    }
    draw();
  })();

  /* =========================================================
     11. #w-rc-beta：頻率方程式 I/I₀ = β tan β (3.11)
     ========================================================= */
  (() => {
    const body = SD.widget('w-rc-beta'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'r', label: '慣量比 I / I₀', min: 0.02, max: 1.50, step: 0.01, value: 0.10, fmt: v => v.toFixed(2) },
      { k: 'L', label: '試體高度 L (mm)', min: 50, max: 200, step: 5, value: 100, fmt: v => v.toFixed(0) },
      { k: 'rho', label: '試體密度 ρ (t/m³)', min: 1.5, max: 2.2, step: 0.01, value: 1.95, fmt: v => v.toFixed(2) },
      { k: 'fr', label: '實測共振頻率 f_r (Hz)', min: 40, max: 220, step: 0.5, value: 108.3, fmt: v => v.toFixed(1) }
    ], draw);
    SD.buttons(body, '典型比值', [
      { label: '範例 (I/I₀ = 0.10, f_r = 108.3 Hz)', on: () => ctl.set({ r: 0.10, L: 100, rho: 1.95, fr: 108.3 }) },
      { label: '極大端質量 (I/I₀ = 0.02)', on: () => ctl.set({ r: 0.02 }) },
      { label: '中等端質量 (I/I₀ = 0.50)', on: () => ctl.set({ r: 0.50 }) },
      { label: '輕端質量 (I/I₀ = 1.20)', on: () => ctl.set({ r: 1.20 }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pEq = new Plot(grid, { h: 270, xlim: [0, 1.4], ylim: [0, 2.0], xlabel: '無因次頻率參數 β = ωL / V_s (rad)', ylabel: '函數值' });
    pEq.legend([
      { label: '精確左式 β tan β', c: 'blue', w: 2.2 },
      { label: '小參數近似 β²', c: 'gray', dash: true },
      { label: '目標水平線 I / I₀', c: 'red', dash: '3 3' },
      { label: '精確解 β₁', c: 'red', line: false, marker: 'circle', fill: true }
    ]);
    const pMode = new Plot(grid, { h: 270, xlim: [0, 1.05], ylim: [0, 1], xlabel: '正規化扭轉角 θ(z) / θ(L)', ylabel: '正規化高度 z / L（底端固定 z=0）' });
    pMode.legend([
      { label: '真實模態 sin(β₁z/L) / sin β₁', c: 'blue', w: 2.2 },
      { label: '均勻剪應變直線 z / L', c: 'gray', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'b1', label: '精確根 β₁ (rad)', hl: true },
      { k: 'b0', label: '近似 √(I/I₀) (rad)' },
      { k: 'lrat', label: '試體高 / 波長 L/λ = β₁/(2π)' },
      { k: 'vs', label: '波速 V_s (m/s)', hl: true },
      { k: 'G', label: '剪力模數 G (MPa)', hl: true },
      { k: 'errG', label: '用近似根算 G 的誤差' }
    ]);

    function draw() {
      const r = ctl.r, Lm = ctl.L / 1000, rho = ctl.rho, fr = ctl.fr;
      const b1 = solveBeta1(r), b0 = Math.sqrt(r);
      const vs = TAU * fr * Lm / b1;
      const G = rho * vs * vs / 1000;
      const G_approx = rho * (TAU * fr * Lm / b0) ** 2 / 1000;

      pEq.setLim([0, 1.4], [0, Math.max(0.5, Math.min(3.0, r * 2.2))]);
      pEq.clear();
      pEq.fn(b => b * Math.tan(b), { c: 'blue', w: 2.2, range: [0, 1.38], n: 250 });
      pEq.fn(b => b * b, { c: 'gray', dash: '4 3', w: 1.5, range: [0, 1.38] });
      pEq.hline(r, { c: 'red', dash: '3 3', w: 1.3 });
      pEq.vline(b1, { c: 'red', dash: '3 3', w: 1 });
      pEq.markers([b1], [r], { c: 'red', fill: true, r: 5, over: true });
      pEq.markers([Math.min(1.38, b0)], [r], { c: 'gray', fill: false, r: 4.5, over: true });

      pMode.clear();
      const zs = linspace(0, 1, 100);
      const ths = zs.map(z => Math.sin(b1 * z) / Math.sin(b1));
      pMode.line(zs, zs, { c: 'gray', dash: '4 3', w: 1.5 });
      pMode.line(ths, zs, { c: 'blue', w: 2.2 });

      st({
        b1: f(b1, 4),
        b0: f(b0, 4),
        lrat: f(b1 / TAU, 4),
        vs: f(vs, 1),
        G: f(G, 2),
        errG: f((G_approx - G) / G * 100, 1) + ' %'
      });
    }
    draw();
  })();

  /* =========================================================
     12. #w-rc-sweep：共振柱掃頻與阻尼實驗室 (3.12)
     ========================================================= */
  (() => {
    const body = SD.widget('w-rc-sweep'); if (!body) return;
    let tabIdx = 0, p1 = null, p2 = null, ctl = null, st = null;
    SD.tabs(body, ['① 線性掃頻 + 自由衰減', '② 大應變軟化與跳躍（半功率法高估）', '③ 掃頻間距 Δf 實驗'], i => {
      tabIdx = i;
      if (p1) { syncUI(); draw(); }
    });

    ctl = SD.controls(body, [
      { k: 'xi', label: '真實黏性阻尼比 ξ (%)', min: 0.5, max: 12, step: 0.25, value: 2.0, fmt: v => v.toFixed(2) + ' %' },
      { k: 'soft', label: '大應變軟化程度 α（頁籤 ②）', min: 0, max: 0.22, step: 0.01, value: 0.12, fmt: v => v.toFixed(2) },
      { k: 'df', label: '掃頻步長 Δf (Hz，頁籤 ③)', min: 0.1, max: 4.0, step: 0.1, value: 1.5, fmt: v => v.toFixed(1) }
    ], draw);

    const grid = SD.div(body, 'plot-row plot-grid');
    p1 = new Plot(grid, { h: 270, xlim: [85, 115], ylim: [0, 30], xlabel: '激振頻率 f (Hz)', ylabel: '動態放大倍數 D' });
    p2 = new Plot(grid, { h: 270, xlim: [0, 8], ylim: [-1.05, 1.05], xlabel: '自由振動週期數 t / T_D', ylabel: '正規化自由衰減位移 u / u₀' });
    st = SD.stats(body, [
      { k: 'trueXi', label: '真實阻尼比 ξ', hl: true },
      { k: 'magXi', label: '共振放大法 ξ_{mag}', hl: true },
      { k: 'hpXi', label: '半功率頻寬法 ξ_{hp}', hl: true },
      { k: 'fdXi', label: '自由衰減法 ξ_{fd}', hl: true },
      { k: 'note', label: '診斷結論' }
    ]);

    const fn0 = 100;

    function syncUI() {
      const boxes = body.querySelectorAll('.ctrl');
      if (boxes.length >= 3) {
        boxes[1].style.display = tabIdx === 1 ? '' : 'none';
        boxes[2].style.display = tabIdx === 2 ? '' : 'none';
      }
    }

    function draw() {
      const xi = ctl.xi / 100;
      const pk = SDOF.peak(xi);
      p1.setLim([85, 115], [0, Math.ceil(pk.D * 1.15 / 5) * 5]);
      p1.clear(); p2.clear();

      if (tabIdx === 0) {
        p1.fn(fq => SDOF.D(fq / fn0, xi), { c: 'blue', w: 2.2, n: 400 });
        const fr = fn0 * pk.b, lvl = pk.D / Math.SQRT2;
        const f1 = fn0 * Math.sqrt(1 - 2 * xi * xi - 2 * xi * Math.sqrt(1 - xi * xi));
        const f2 = fn0 * Math.sqrt(1 - 2 * xi * xi + 2 * xi * Math.sqrt(1 - xi * xi));
        p1.hline(lvl, { c: 'gray', dash: '4 3' });
        p1.vline(f1, { c: 'red', dash: '3 3' }); p1.vline(f2, { c: 'red', dash: '3 3' });
        p1.markers([fr, f1, f2], [pk.D, lvl, lvl], { c: 'red', fill: true, r: 4.5, over: true });
        p1.text(86, lvl + pk.D * 0.04, `−3 dB：f₁=${f1.toFixed(2)}, f₂=${f2.toFixed(2)} Hz`, { cls: 'txt-red' });

        p2.hline(0, { c: 'gray', dash: null, w: 0.6 });
        p2.fn(t => Math.exp(-TAU * xi * t / Math.sqrt(1 - xi * xi)), { c: 'gray', dash: '4 3' });
        p2.fn(t => -Math.exp(-TAU * xi * t / Math.sqrt(1 - xi * xi)), { c: 'gray', dash: '4 3' });
        p2.fn(t => SDOF.free(t / Math.sqrt(1 - xi * xi), TAU, xi, 1, 0), { c: 'blue', w: 2, n: 500 });
        const ks = [0, 1, 2, 3, 4, 5, 6, 7];
        const pks = ks.map(k => Math.exp(-TAU * xi * k / Math.sqrt(1 - xi * xi)));
        p2.markers(ks, pks, { c: 'red', fill: true, r: 4, over: true });

        const xiHp = (f2 - f1) / (f2 + f1);
        const delta = TAU * xi / Math.sqrt(1 - xi * xi);
        const xiFd = delta / Math.sqrt(TAU * TAU + delta * delta);
        st({
          trueXi: f(xi * 100, 2) + ' %',
          magXi: f(SDOF.xiFromDmax(pk.D) * 100, 2) + ' %',
          hpXi: f(xiHp * 100, 2) + ' %',
          fdXi: f(xiFd * 100, 2) + ' %',
          note: '線性小應變：三種方法一致'
        });
      } else if (tabIdx === 1) {
        const soft = ctl.soft;
        p1.fn(fq => SDOF.D(fq / fn0, xi), { c: 'gray', dash: '4 3', w: 1.4, n: 350 });
        const As = linspace(1.2, pk.D, 180);
        const fLeft = [], fRight = [], bboneF = [], bboneA = [];
        for (const A of As) {
          const feq = fn0 * (1 - soft * (A / pk.D) ** 1.2);
          bboneF.push(feq); bboneA.push(A);
          const rad = Math.max(0, (1 / A) ** 2 - (2 * xi) ** 2);
          const d = Math.sqrt(rad);
          fLeft.push(feq * Math.sqrt(Math.max(0.2, 1 - 2 * xi * xi - d)));
          fRight.push(feq * Math.sqrt(1 - 2 * xi * xi + d));
        }
        p1.line(bboneF, bboneA, { c: 'gray', dash: '2 3', w: 1.3 });
        p1.line(fLeft, As, { c: 'red', w: 2.2 });
        p1.line(fRight, As, { c: 'blue', w: 2.2 });
        const lvl = pk.D / Math.SQRT2;
        let idx = 0;
        for (let i = 0; i < As.length; i++) if (Math.abs(As[i] - lvl) < Math.abs(As[idx] - lvl)) idx = i;
        const f1s = fLeft[idx], f2s = fRight[idx];
        const xiHpSoft = (f2s - f1s) / (f2s + f1s);
        p1.hline(lvl, { c: 'gray', dash: '3 3' });
        p1.markers([f1s, f2s], [lvl, lvl], { c: 'red', fill: true, r: 4.5, over: true });
        p1.text(86, pk.D * 0.92, '紅：下行掃頻 / 藍：上行掃頻（軟化左傾）', { cls: 'txt-red' });

        p2.hline(0, { c: 'gray', dash: null, w: 0.6 });
        p2.fn(t => SDOF.free(t, TAU, xi, 1, 0), { c: 'blue', w: 2, n: 500 });
        p2.text(4, 0.85, '自由衰減不受掃頻跳躍影響', { anchor: 'middle', cls: 'txt-blue' });

        st({
          trueXi: f(xi * 100, 2) + ' %',
          magXi: f(xi * 100, 2) + ' %',
          hpXi: f(xiHpSoft * 100, 2) + ' % (偏高!)',
          fdXi: f(xi * 100, 2) + ' % (建議採用)',
          note: soft > 0.02 ? `半功率法高估 ${f((xiHpSoft / xi - 1) * 100, 0)}%` : '軟化極小'
        });
      } else {
        const df = ctl.df;
        p1.fn(fq => SDOF.D(fq / fn0, xi), { c: 'gray', dash: '4 3', w: 1.4, n: 400 });
        const xs = [], ys = [];
        for (let fq = 85 + 0.37 * df; fq <= 115; fq += df) {
          xs.push(fq); ys.push(SDOF.D(fq / fn0, xi));
        }
        p1.line(xs, ys, { c: 'blue', w: 2 });
        p1.markers(xs, ys, { c: 'blue', fill: true, r: 4, over: true });

        let imax = 0;
        for (let i = 1; i < ys.length; i++) if (ys[i] > ys[imax]) imax = i;
        const Dmeas = ys[imax];
        const lvl = Dmeas / Math.SQRT2;
        const f1 = SDOF.crossing(xs, ys, lvl, imax, -1);
        const f2 = SDOF.crossing(xs, ys, lvl, imax, 1);
        p1.hline(lvl, { c: 'red', dash: '3 3' });
        p1.markers([xs[imax]], [Dmeas], { c: 'red', fill: true, r: 5.5, over: true });

        p2.hline(0, { c: 'gray', dash: null, w: 0.6 });
        p2.fn(t => SDOF.free(t, TAU, xi, 1, 0), { c: 'blue', w: 2, n: 500 });

        const xiMag = SDOF.xiFromDmax(Dmeas);
        const xiHp = isFinite(f1) && isFinite(f2) ? (f2 - f1) / (f2 + f1) : NaN;
        st({
          trueXi: f(xi * 100, 2) + ' %',
          magXi: f(xiMag * 100, 2) + ' %',
          hpXi: isFinite(xiHp) ? f(xiHp * 100, 2) + ' %' : '點太少',
          fdXi: f(xi * 100, 2) + ' %',
          note: `峰頂漏失 ${f((1 - Dmeas / pk.D) * 100, 1)}% → 高估阻尼`
        });
      }
    }
    syncUI();
    draw();
  })();

  /* =========================================================
     13. #w-hw3a：HW3 Q1(a) 由原始資料重算 stage 1 (3.13)
     ========================================================= */
  (() => {
    const body = SD.widget('w-hw3a'); if (!body) return;
    const rc = window.CH3_HW3_RC;
    if (!rc) return;
    const ctl = SD.controls(body, [
      { k: 'peak', label: '共振頻率的求法', type: 'select', value: 'fit', options: [['fit', 'SDOF 擬合（報告採用）'], ['parab', '三點拋物線內插'], ['max', '最大樣本點']] },
      { k: 'lsq', label: '自由衰減改用 6 點迴歸（預設為題目公式：首尾兩峰）', type: 'check', value: false },
      { k: 'nPeak', label: '首尾兩峰法的結束峰號 N', min: 2, max: 6, step: 1, value: 6, fmt: v => 'k = ' + v.toFixed(0) }
    ], draw);

    // 扭矩激振、加速度計量測：A(f) = A0 r² / √[(1 − r²)² + (2ξr)²]，r = f/fn（與 A1_resonant_column.py 相同）
    const accel = (fq, a0, fn, xi) => { const r = fq / fn; return a0 * r * r / Math.sqrt((1 - r * r) ** 2 + (2 * xi * r) ** 2); };
    function nelderMead(fun, x0, steps, iters) {
      let pts = [x0, ...x0.map((_, i) => x0.map((v, j) => j === i ? v + steps[i] : v))].map(x => ({ x, v: fun(x) }));
      for (let it = 0; it < iters; it++) {
        pts.sort((a, b) => a.v - b.v);
        const n = x0.length, c = x0.map((_, j) => pts.slice(0, n).reduce((s, p) => s + p.x[j], 0) / n);
        const mv = (a, t) => c.map((cj, j) => cj + t * (a.x[j] - cj));
        const xr = mv(pts[n], -1), vr = fun(xr);
        if (vr < pts[0].v) { const xe = mv(pts[n], -2), ve = fun(xe); pts[n] = ve < vr ? { x: xe, v: ve } : { x: xr, v: vr }; }
        else if (vr < pts[n - 1].v) pts[n] = { x: xr, v: vr };
        else { const xc = mv(pts[n], 0.5), vc = fun(xc); if (vc < pts[n].v) pts[n] = { x: xc, v: vc }; else pts = pts.map((p, i) => i ? { x: p.x.map((v, j) => pts[0].x[j] + 0.5 * (v - pts[0].x[j])), v: 0 } : p).map(p => ({ x: p.x, v: fun(p.x) })); }
      }
      pts.sort((a, b) => a.v - b.v);
      return pts[0].x;
    }
    const FIT = (() => {
      const fs = rc.sweep.f, as = rc.sweep.a, pk = SDOF.parabolicPeak(fs, as);
      const sse = ([a0, fn, xi]) => fs.reduce((s, fq, i) => s + (accel(fq, a0, fn, xi) - as[i]) ** 2, 0);
      const x0 = [pk.y * 2 * 0.01, pk.x, 0.01];
      let x = nelderMead(sse, x0, [x0[0] * 0.2, 0.5, 0.003], 600);
      x = nelderMead(sse, x, [x[0] * 0.02, 0.05, 0.0005], 600);
      return { a0: x[0], fn: x[1], xi: x[2] };
    })();

    const grid = SD.div(body, 'plot-row plot-grid');
    const pSw = new Plot(grid, { h: 280, xlim: [96, 120], ylim: [0, 65], xlabel: '頻率 f (Hz)', ylabel: '加速度計振幅 A (mV)' });
    pSw.legend([
      { label: '實測掃頻資料 (97 點)', c: 'ink', dash: true, marker: 'circle' },
      { label: 'SDOF 擬合曲線', c: 'blue', w: 2 },
      { label: '−3 dB 半功率線', c: 'red', dash: true }
    ]);
    const pFd = new Plot(grid, { h: 280, xlim: [0.5, 6.5], ylim: [55, 95], xlabel: '自由振動波峰編號 k', ylabel: '波峰振幅 X_k (mV)' });
    pFd.legend([
      { label: '實測波峰 X_k (6 點)', c: 'ink', dash: true, marker: 'circle' },
      { label: '指數衰減 X₁ e^{−δ(k−1)}', c: 'blue', w: 2 }
    ]);
    const st = SD.stats(body, [
      { k: 'fr', label: '共振頻率 f_r (Hz)', hl: true },
      { k: 'b1', label: 'β₁ (I/I₀ = 0.10)' },
      { k: 'vs', label: '剪力波速 V_s (m/s)', hl: true },
      { k: 'G', label: '剪力模數 G_{max} (MPa)', hl: true },
      { k: 'Dhp', label: '半功率法 D_{hp} (%)', hl: true },
      { k: 'Dfd', label: '自由衰減法 D_{fd} (%)', hl: true }
    ]);

    function draw() {
      const fs = rc.sweep.f, as = rc.sweep.a;
      const pkP = SDOF.parabolicPeak(fs, as);
      let imax = pkP.i;
      const fr = ctl.peak === 'fit' ? FIT.fn : ctl.peak === 'parab' ? pkP.x : fs[imax];
      const Amax = ctl.peak === 'max' ? as[imax] : pkP.y;
      const lvl = Amax / Math.SQRT2;
      const f1 = SDOF.crossing(fs, as, lvl, imax, -1);
      const f2 = SDOF.crossing(fs, as, lvl, imax, 1);
      const Dhp = (f2 - f1) / (f2 + f1);

      const X = rc.stages[0].X;
      const ks = [1, 2, 3, 4, 5, 6];
      const lns = X.map(Math.log);
      let delta, a0;
      if (ctl.lsq) {
        const km = 3.5, ym = lns.reduce((s, v) => s + v, 0) / 6;
        let num = 0, den = 0;
        for (let i = 0; i < 6; i++) { num += (ks[i] - km) * (lns[i] - ym); den += (ks[i] - km) ** 2; }
        delta = -num / den;
        a0 = ym + delta * km;
      } else {
        const N = Math.round(ctl.nPeak);
        delta = Math.log(X[0] / X[N - 1]) / (N - 1);
        a0 = Math.log(X[0]) + delta;
      }
      const Dfd = delta / Math.sqrt(TAU * TAU + delta * delta);

      const b1 = solveBeta1(0.10);
      // 擬合得到的就是自然頻率；峰值法量到的是加速度峰，阻尼約 1 % 時與 fn 相差 < 0.02 %，直接取用
      const fn = fr;
      const vs = TAU * fn * 0.10 / b1;
      const G = 1.95 * vs * vs / 1000;

      pSw.clear();
      pSw.line(fs, as, { c: 'ink', dash: '4 3', w: 0.9 });
      pSw.markers(fs, as, { c: 'ink', r: 2.6 });
      pSw.fn(fq => accel(fq, FIT.a0, FIT.fn, FIT.xi), { c: 'blue', w: 2, n: 400 });
      pSw.hline(lvl, { c: 'red', dash: '4 3', w: 1.3 });
      pSw.vline(f1, { c: 'red', dash: '3 3' }); pSw.vline(f2, { c: 'red', dash: '3 3' });
      pSw.markers([fr, f1, f2], [Amax, lvl, lvl], { c: 'red', fill: true, r: 4.8, over: true });
      pSw.text(97, lvl + 2.5, `f₁ = ${f1.toFixed(2)} Hz, f₂ = ${f2.toFixed(2)} Hz`, { cls: 'txt-red' });

      pFd.clear();
      pFd.fn(k => Math.exp(a0 - delta * k), { c: 'blue', w: 2, range: [0.8, 6.2] });
      pFd.line(ks, X, { c: 'ink', dash: '4 3', w: 0.9 });
      pFd.markers(ks, X, { c: 'ink', r: 4.8, over: true });
      pFd.text(3.4, 88, `δ = ${delta.toFixed(4)} → D = ${(Dfd * 100).toFixed(2)}%`, { cls: 'txt-red' });

      st({
        fr: f(fr, 2),
        b1: f(b1, 4),
        vs: f(vs, 1),
        G: f(G, 1),
        Dhp: f(Dhp * 100, 2) + ' %',
        Dfd: f(Dfd * 100, 2) + ' %'
      });
    }
    draw();
  })();

  /* =========================================================
     14. #w-bender：彎曲元件走時拾取與近場效應 (3.14)
     ========================================================= */
  (() => {
    const body = SD.widget('w-bender'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'L', label: '尖端距離 Lₜₜ (mm)', min: 50, max: 200, step: 5, value: 100, fmt: v => v.toFixed(0) },
      { k: 'vs', label: '真實剪力波速 Vₛ (m/s)', min: 100, max: 350, step: 10, value: 180, fmt: v => v.toFixed(0) },
      { k: 'freq', label: '激振正弦頻率 f (kHz)', min: 1.5, max: 12.0, step: 0.5, value: 5.0, fmt: v => v.toFixed(1) },
      { k: 'nu', label: 'Poisson 比 ν', min: 0.15, max: 0.45, step: 0.05, value: 0.30, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '試驗條件', [
      { label: '良態條件 (Lₜₜ/λ = 2.8 ≥ 2)', on: () => ctl.set({ L: 100, vs: 180, freq: 5.0 }) },
      { label: '近場效應明顯 (f = 2 kHz, Lₜₜ/λ = 1.1)', on: () => ctl.set({ L: 100, vs: 180, freq: 2.0 }) },
      { label: '高頻遠場 (f = 10 kHz, Lₜₜ/λ = 5.6)', on: () => ctl.set({ L: 100, vs: 180, freq: 10.0 }) }
    ]);
    const p = new Plot(body, { h: 270, xlim: [0, 1.5], ylim: [-1.4, 2.5], xlabel: '時間 t (ms)', ylabel: '正規化訊號（上：發射波；下：接收波）' });
    p.legend([
      { label: '發射端單週期正弦脈衝', c: 'gray', w: 1.8 },
      { label: '接收端訊號（含近場前驅 + S 波）', c: 'blue', w: 2.2 },
      { label: '偏折初達點 (含近場偏早)', c: 'red', line: false, marker: 'tri', fill: true },
      { label: '第一波峰 / 真實到達時', c: 'ink', line: false, marker: 'circle', fill: true }
    ]);
    const st = SD.stats(body, [
      { k: 'ts', label: '真實 S 波走時 t_S (ms)', hl: true },
      { k: 'lrat', label: '距離波長比 L_{tt} / λ', hl: true },
      { k: 'vFirst', label: '偏折初達法 V_s (m/s)', hl: true },
      { k: 'vPeak', label: '峰對峰法 V_s (m/s)' },
      { k: 'vCC', label: '互相關法 V_s (m/s)' },
      { k: 'warn', label: '近場效應診斷' }
    ]);

    function draw() {
      const Lm = ctl.L / 1000, vs = ctl.vs, fq = ctl.freq * 1000, nu = ctl.nu;
      const vp = vs * Math.sqrt((2 - 2 * nu) / (1 - 2 * nu));
      const ts = (Lm / vs) * 1000;
      const tp = (Lm / vp) * 1000;
      const T = 1000 / fq;
      const lam = vs / fq;
      const lrat = Lm / lam;

      const tmax = Math.max(1.2, ts + 2.2 * T);
      p.setLim([0, Math.ceil(tmax * 5) / 5], [-1.4, 2.5]);
      p.clear();
      p.hline(1.6, { c: 'gray', dash: null, w: 0.6 });
      p.hline(0, { c: 'gray', dash: null, w: 0.6 });

      const src = t => (t >= 0 && t <= T) ? Math.sin(TAU * t / T) * (0.5 - 0.5 * Math.cos(TAU * t / T)) : 0;
      const nfAmp = Math.min(0.45, 0.38 / Math.max(0.5, lrat));
      const rec = t => {
        let y = 0;
        if (t >= ts && t <= ts + 1.5 * T) {
          const u = (t - ts) / T;
          y += Math.sin(TAU * u) * Math.exp(-1.1 * (u - 0.35) ** 2);
        }
        if (t >= tp && t <= ts + 0.3 * T) {
          const u = (t - tp) / (ts - tp + 0.15 * T);
          y -= nfAmp * Math.sin(Math.PI * Math.min(1, u));
        }
        return y;
      };

      p.fn(t => 1.6 + 0.65 * src(t), { c: 'gray', w: 1.8, n: 400 });
      p.fn(rec, { c: 'blue', w: 2.2, n: 500 });

      const tFirst = lrat < 2.2 ? ts - (ts - tp) * Math.min(0.65, 0.55 / lrat) : ts;
      const tSrcPeak = 0.27 * T;
      const tRecPeak = ts + 0.27 * T;
      const dtPeak = tRecPeak - tSrcPeak;
      const dtCC = ts * (1 + 0.004 / lrat);

      p.vline(ts, { c: 'ink', dash: '3 3', w: 1.2 });
      p.vline(tFirst, { c: 'red', dash: '3 3', w: 1.2 });
      p.markers([tFirst], [rec(tFirst)], { c: 'red', shape: 'tri', fill: true, r: 5, over: true });
      p.markers([tSrcPeak, tRecPeak, ts], [1.6 + 0.65 * src(tSrcPeak), rec(tRecPeak), 0], { c: 'ink', fill: true, r: 4.5, over: true });

      const vFirst = Lm / (tFirst / 1000);
      const vPeak = Lm / (dtPeak / 1000);
      const vCC = Lm / (dtCC / 1000);
      st({
        ts: f(ts, 3),
        lrat: f(lrat, 2),
        vFirst: `${f(vFirst, 0)} (${lrat < 2 ? '+' + f((vFirst / vs - 1) * 100, 0) + '%' : '準確'})`,
        vPeak: f(vPeak, 0),
        vCC: f(vCC, 0),
        warn: lrat < 2 ? 'Lₜₜ/λ < 2：近場前驅拉早初達，應提高頻率' : 'Lₜₜ/λ ≥ 2：近場影響小'
      });
    }
    draw();
  })();

});
