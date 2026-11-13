/* =========================================================
   Ch6 地盤反應分析：互動元件
   ========================================================= */
document.addEventListener('sd:ready', () => {
  const TAU = Math.PI * 2, f = SD.f;
  const deg = r => r * 180 / Math.PI;
  const rad = d => d * Math.PI / 180;
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

  function clearG(g) { while (g.firstChild) g.removeChild(g.firstChild); }
  function svgText(g, x, y, s, cls, anchor) {
    const t = svgEl('text', { x, y, class: cls || 'txt', 'text-anchor': anchor || 'start' }, g);
    svgRich(t, s);
    return t;
  }

  /* ---- 複數與波傳輔助函式 ---- */
  function cCos(xr, xi) {
    return [Math.cos(xr) * Math.cosh(xi), -Math.sin(xr) * Math.sinh(xi)];
  }
  function cSin(xr, xi) {
    return [Math.sin(xr) * Math.cosh(xi), Math.cos(xr) * Math.sinh(xi)];
  }
  function cDiv(ar, ai, br, bi) {
    const d = br * br + bi * bi || 1e-30;
    return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d];
  }
  function cMul(ar, ai, br, bi) {
    return [ar * br - ai * bi, ar * bi + ai * br];
  }

  function kStarH(kH, xi) {
    const m2 = Math.hypot(1, 2 * xi);
    const sr = Math.sqrt(0.5 * (m2 + 1)), si = Math.sqrt(0.5 * (m2 - 1));
    return cDiv(kH, 0, sr, si);
  }
  function tfSingle(kH, xi, az = 0) {
    const [kr, ki] = kStarH(kH, xi);
    const [cr, ci] = cCos(kr, ki);
    const [sr, si] = cSin(kr, ki);
    const [f2r, f2i] = cDiv(1, 0, cr, ci);
    const m2 = Math.hypot(1, 2 * xi);
    const ar = az * Math.sqrt(0.5 * (m2 + 1)), ai = az * Math.sqrt(0.5 * (m2 - 1));
    const [asr, asi] = cMul(ar, ai, sr, si);
    const denR = cr - asi, denI = ci + asr;
    const [f3r, f3i] = cDiv(1, 0, denR, denI);
    return {
      f2r, f2i, f2: Math.hypot(f2r, f2i),
      f3r, f3i, f3: Math.hypot(f3r, f3i),
      ub: Math.hypot(cr, ci) * Math.hypot(f3r, f3i)
    };
  }

  function tfMulti(freq, layers) {
    const w = TAU * freq;
    const n = layers.length;
    const vsC = layers.map(l => {
      const m2 = Math.hypot(1, 2 * (l.xi || 0));
      return [l.vs * Math.sqrt(0.5 * (m2 + 1)), l.vs * Math.sqrt(0.5 * (m2 - 1))];
    });
    let ar = 1, ai = 0, br = 1, bi = 0;
    for (let m = 0; m < n - 1; m++) {
      const h = layers[m].h;
      const [vr, vi] = vsC[m];
      const [kr, ki] = cDiv(w * h, 0, vr, vi);
      const ep = Math.exp(-ki), em = Math.exp(ki);
      const epr = ep * Math.cos(kr), epi = ep * Math.sin(kr);
      const emr = em * Math.cos(kr), emi = -em * Math.sin(kr);
      const rRatio = layers[m].rho / layers[m + 1].rho;
      const [alr, ali] = cDiv(rRatio * vr, rRatio * vi, vsC[m + 1][0], vsC[m + 1][1]);
      const p1r = 0.5 * (1 + alr), p1i = 0.5 * ali;
      const m1r = 0.5 * (1 - alr), m1i = -0.5 * ali;
      const [aeR, aeI] = cMul(ar, ai, epr, epi);
      const [beR, beI] = cMul(br, bi, emr, emi);
      const [t1r, t1i] = cMul(p1r, p1i, aeR, aeI);
      const [t2r, t2i] = cMul(m1r, m1i, beR, beI);
      const [t3r, t3i] = cMul(m1r, m1i, aeR, aeI);
      const [t4r, t4i] = cMul(p1r, p1i, beR, beI);
      ar = t1r + t2r; ai = t1i + t2i;
      br = t3r + t4r; bi = t3i + t4i;
    }
    const [f3r, f3i] = cDiv(1, 0, ar, ai);
    const [f2r, f2i] = cDiv(2, 0, ar + br, ai + bi);
    return { f3r, f3i, f3: Math.hypot(f3r, f3i), f2r, f2i, f2: Math.hypot(f2r, f2i) };
  }

  function fft(re, im, sign) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = sign * TAU / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let j = 0; j < len / 2; j++) {
          const a = i + j, b = a + len / 2;
          const xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
          re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
          const ncr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = ncr;
        }
      }
    }
  }

  /* =========================================================
     1. #w-snell：Snell 定律射線折射 (6.1)
     ========================================================= */
  (() => {
    const body = SD.widget('w-snell'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'th0', label: '基底岩盤入射角 θ₀ (°)', min: 5, max: 75, step: 1, value: 40, fmt: v => v.toFixed(0) + '°' },
      { k: 'v0', label: '岩盤波速 V₀ (m/s)', min: 800, max: 2500, step: 50, value: 1500, fmt: v => v.toFixed(0) },
      { k: 'v1', label: '深部土層 V₁ (m/s)', min: 300, max: 900, step: 20, value: 500, fmt: v => v.toFixed(0) },
      { k: 'v2', label: '中部土層 V₂ (m/s)', min: 180, max: 500, step: 10, value: 300, fmt: v => v.toFixed(0) },
      { k: 'v3', label: '表層軟土 V₃ (m/s)', min: 60, max: 250, step: 10, value: 150, fmt: v => v.toFixed(0) }
    ], render);
    SD.buttons(body, '範例', [
      { label: '預設 (40° 入射, 1500 → 150 m/s → 3.7°)', on: () => ctl.set({ th0: 40, v0: 1500, v1: 500, v2: 300, v3: 150 }) },
      { label: '斜入射 (70° 入射)', on: () => ctl.set({ th0: 70, v0: 1500, v1: 500, v2: 300, v3: 150 }) },
      { label: '極軟表土 (V₃ = 80 m/s)', on: () => ctl.set({ th0: 40, v0: 1500, v1: 500, v2: 300, v3: 80 }) }
    ]);
    const W = 680, H = 260;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot-svg' });
    body.appendChild(svg);
    const g = svgEl('g', null, svg);
    const st = SD.stats(body, [
      { k: 'th0', label: '岩盤入射角 θ₀' },
      { k: 'th1', label: '深層折射角 θ₁' },
      { k: 'th2', label: '中層折射角 θ₂' },
      { k: 'th3', label: '表層折射角 θ₃', hl: true },
      { k: 'cx', label: '水平視速度 c_x (m/s)' },
      { k: 'cos3', label: '鉛直分量比 cos θ₃', hl: true }
    ]);

    function render() {
      clearG(g);
      const vs = [ctl.v3, ctl.v2, ctl.v1, ctl.v0];
      const names = ['表層軟土 V₃', '中部土層 V₂', '深部土層 V₁', '基底岩盤 V₀'];
      const p = Math.sin(rad(ctl.th0)) / ctl.v0;
      const ths = vs.map(v => Math.asin(Math.min(0.999, v * p)));

      const x0 = 50, x1 = 630, yTop = 28, dy = 52;
      for (let i = 0; i < 4; i++) {
        const yA = yTop + i * dy;
        const alpha = 0.05 + i * 0.04;
        svgEl('rect', { x: x0, y: yA, width: x1 - x0, height: dy, fill: `rgba(31,95,168,${alpha})` }, g);
        svgEl('line', { x1: x0, x2: x1, y1: yA, y2: yA, class: i === 0 ? 's-ink' : 's-gray', 'stroke-width': i === 0 ? 1.8 : 1, 'stroke-dasharray': i === 0 ? null : '5 4' }, g);
        svgText(g, x0 + 10, yA + 22, `${names[i]} = ${vs[i]} m/s`, 'txt');
        svgText(g, x1 - 10, yA + 22, `θ = ${deg(ths[i]).toFixed(1)}°`, i === 0 ? 'txt-red' : 'txt-muted', 'end');
      }

      let rx = 210, ry = yTop + 4 * dy;
      const pts = [[rx, ry]];
      for (let i = 3; i >= 0; i--) {
        const yNext = yTop + i * dy;
        const xNext = rx + dy * Math.tan(ths[i]);
        svgEl('line', { x1: rx, x2: rx, y1: ry - dy, y2: ry, class: 's-gray', 'stroke-width': 1, 'stroke-dasharray': '2 3' }, g);
        const mx = 0.5 * (rx + xNext), my = 0.5 * (ry + yNext);
        const wx = 18 * Math.cos(ths[i]), wy = 18 * Math.sin(ths[i]);
        svgEl('line', { x1: mx - wx, x2: mx + wx, y1: my - wy, y2: my + wy, class: 's-red', 'stroke-width': 1.6 }, g);
        rx = xNext; ry = yNext;
        pts.push([rx, ry]);
      }
      const d = pts.map((pt, idx) => (idx ? 'L' : 'M') + pt[0].toFixed(1) + ' ' + pt[1].toFixed(1)).join('');
      svgEl('path', { d, fill: 'none', class: 's-blue', 'stroke-width': 2.6 }, g);
      pts.forEach(([px, py]) => svgEl('circle', { cx: px, cy: py, r: 3.8, class: 'f-blue' }, g));

      st({
        th0: f(ctl.th0, 1) + '°',
        th1: f(deg(ths[2]), 1) + '°',
        th2: f(deg(ths[1]), 1) + '°',
        th3: f(deg(ths[0]), 1) + '°',
        cx: f(1 / p, 0),
        cos3: f(Math.cos(ths[0]), 4)
      });
    }
    render();
  })();

  /* =========================================================
     2. #w-outcrop：地表、岩盤頂、露頭各動多少 (6.2)
     ========================================================= */
  (() => {
    const body = SD.widget('w-outcrop'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'az', label: '阻抗比 αz = ρsVs / (ρrVr)', min: 0.0, max: 0.8, step: 0.02, value: 0.25, fmt: v => v.toFixed(2) },
      { k: 'xi', label: '土層阻尼比 ξ (%)', min: 0, max: 15, step: 0.5, value: 2.0, fmt: v => v.toFixed(1) + ' %' },
      { k: 'khp', label: '檢視無因次頻率 kH / π', min: 0.05, max: 2.5, step: 0.05, value: 0.50, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '特徵頻率', [
      { label: '一階共振 (kH = π/2)', on: () => ctl.set({ khp: 0.50, xi: 2.0, az: 0.25 }) },
      { label: '反共振鎖死岩盤頂 (kH = π, ξ = 0)', on: () => ctl.set({ khp: 1.00, xi: 0.0, az: 0.25 }) },
      { label: '二階共振 (kH = 3π/2)', on: () => ctl.set({ khp: 1.50, xi: 2.0, az: 0.25 }) },
      { label: '剛性岩盤極限 (αz = 0)', on: () => ctl.set({ az: 0.0, xi: 5.0, khp: 0.50 }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pFreq = new Plot(grid, { h: 280, xlim: [0, 2.5], ylim: [0, 10], xlabel: '無因次頻率 kH / π（0.5, 1.5 為土層共振；1.0, 2.0 為反共振）', ylabel: '相對入射波 A_r 的振幅比' });
    pFreq.legend([
      { label: '土層地表 |u_s| / A_r = 2|F₃|', c: 'blue', w: 2.2 },
      { label: '岩盤頂 (土層底) |u_b| / A_r', c: 'red', dash: true, w: 2 },
      { label: '岩盤露頭 |u_{outcrop}| / A_r ≡ 2', c: 'gray', dash: '2 3' }
    ]);
    const pProf = new Plot(grid, { h: 280, xlim: [0, 10], ylim: [-1, 0], yfmt: v => Math.abs(v).toFixed(1), xlabel: '振幅包絡 |u(z)| / A_r', ylabel: '正規化深度 z / H（0=地表，1=岩盤頂）' });
    pProf.legend([
      { label: '土柱內振幅包絡 |u(z)|/A_r', c: 'blue', w: 2.2 },
      { label: '岩盤露頭基準 2A_r', c: 'gray', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'out', label: '岩盤露頭 |u_{out}| / A_r' },
      { k: 'us', label: '土層地表 |u_s| / A_r', hl: true },
      { k: 'ub', label: '岩盤頂 |u_b| / A_r', hl: true },
      { k: 'f2', label: '|F₂| = |u_s / u_b|', hl: true },
      { k: 'f3', label: '|F₃| = |u_s / u_{out}|', hl: true }
    ]);

    function draw() {
      const az = ctl.az, xi = ctl.xi / 100, kh = ctl.khp * Math.PI;
      const cur = tfSingle(kh, xi, az);
      const ymax = Math.min(16, Math.max(6, Math.ceil(2 * tfSingle(Math.PI / 2, Math.max(0.01, xi), az).f3 * 1.15)));

      pFreq.setLim([0, 2.5], [0, ymax]);
      pFreq.clear();
      pFreq.hline(2, { c: 'gray', dash: '3 3', w: 1.4 });
      pFreq.fn(x => Math.min(ymax * 1.2, 2 * tfSingle(x * Math.PI, xi, az).f3), { c: 'blue', w: 2.2, n: 350 });
      pFreq.fn(x => 2 * tfSingle(x * Math.PI, xi, az).ub, { c: 'red', dash: '5 4', w: 2, n: 350 });
      pFreq.vline(ctl.khp, { c: 'gray', dash: '3 3' });
      pFreq.markers([ctl.khp, ctl.khp], [Math.min(ymax, 2 * cur.f3), 2 * cur.ub], { c: 'blue', fill: true, r: 4.5, over: true });

      pProf.setLim([0, ymax], [-1, 0]);
      pProf.clear();
      pProf.vline(2, { c: 'gray', dash: '3 3', w: 1.4 });
      const zs = linspace(0, 1, 120);
      const [kr, ki] = kStarH(kh, xi);
      const us = zs.map(z => {
        const [cr, ci] = cCos(kr * z, ki * z);
        return Math.min(ymax, 2 * cur.f3 * Math.hypot(cr, ci));
      });
      pProf.line(us, zs.map(z => -z), { c: 'blue', w: 2.2 });
      pProf.markers([Math.min(ymax, 2 * cur.f3), 2 * cur.ub], [0, -1], { c: 'red', fill: true, r: 5, over: true });

      st({
        out: '2.00',
        us: f(2 * cur.f3, 2),
        ub: f(2 * cur.ub, 2),
        f2: cur.f2 > 99 ? '∞ (節點)' : f(cur.f2, 2),
        f3: f(cur.f3, 2)
      });
    }
    draw();
  })();

  /* =========================================================
     3. #w-tf1：均勻土層、剛性岩盤的轉換函數 F₁ 與 F₂ (6.3)
     ========================================================= */
  (() => {
    const body = SD.widget('w-tf1'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'H', label: '土層厚度 H (m)', min: 5, max: 100, step: 5, value: 30, fmt: v => v.toFixed(0) },
      { k: 'vs', label: '剪力波速 V_s (m/s)', min: 80, max: 600, step: 10, value: 240, fmt: v => v.toFixed(0) },
      { k: 'xi', label: '阻尼比 ξ (%)', min: 1, max: 20, step: 0.5, value: 5.0, fmt: v => v.toFixed(1) + ' %' },
      { k: 'showF1', label: '疊加無阻尼 |F₁(f)|', type: 'check', value: true }
    ], draw);
    SD.buttons(body, '預設剖面', [
      { label: '基準剖面 (H=30 m, Vs=240 → f₁=2 Hz)', on: () => ctl.set({ H: 30, vs: 240, xi: 5 }) },
      { label: '深厚軟土 (H=60 m, Vs=150 → f₁=0.625 Hz)', on: () => ctl.set({ H: 60, vs: 150, xi: 5 }) },
      { label: '大阻尼 (ξ = 15%)', on: () => ctl.set({ xi: 15 }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pMag = new Plot(grid, { h: 270, xlim: [0, 12], ylim: [0, 16], xlabel: '頻率 f (Hz)', ylabel: '振幅轉換函數 |F(f)|' });
    pMag.legend([
      { label: '有阻尼 |F₂(f)|', c: 'blue', w: 2.2 },
      { label: '無阻尼 |F₁(f)| = 1/|cos(kH)|', c: 'gray', dash: true },
      { label: '共振峰漸近線 2/(πξ·(2n−1))', c: 'red', dash: '2 3' }
    ]);
    const pPh = new Plot(grid, { h: 270, xlim: [0, 12], ylim: [-540, 20], xlabel: '頻率 f (Hz)', ylabel: '展開相位角 ∠F₂(f) (°)' });
    pPh.legend([
      { label: '相位角 ∠F₂(f)', c: 'blue', w: 2.2 },
      { label: '共振相位 −90°, −270°, −450°', c: 'red', line: false, marker: 'circle', fill: true }
    ]);
    const st = SD.stats(body, [
      { k: 'f1', label: '基本頻率 f₁ = V_s/(4H) (Hz)', hl: true },
      { k: 'T1', label: '基本週期 T₁ = 4H/V_s (s)', hl: true },
      { k: 'p1', label: '第 1 峰 |F₂(f₁)| ≈ 2/(πξ)', hl: true },
      { k: 'p2', label: '第 2 峰 |F₂(3f₁)|' },
      { k: 'p3', label: '第 3 峰 |F₂(5f₁)|' }
    ]);

    function draw() {
      const H = ctl.H, vs = ctl.vs, xi = ctl.xi / 100;
      const f1 = vs / (4 * H), fmax = 6 * f1;
      const ymax = Math.min(25, Math.max(8, Math.ceil((2 / (Math.PI * xi)) * 1.15)));

      pMag.setLim([0, fmax], [0, ymax]);
      pMag.clear();
      if (ctl.showF1) {
        pMag.fn(fq => {
          const c = Math.abs(Math.cos(TAU * fq * H / vs));
          return c < 1 / (ymax * 1.2) ? NaN : 1 / c;
        }, { c: 'gray', dash: '4 3', w: 1.3, n: 600 });
      }
      pMag.fn(fq => {
        if (fq < 0.2 * f1) return NaN;
        return Math.min(ymax * 1.1, 1 / (xi * TAU * fq * H / vs));
      }, { c: 'red', dash: '2 3', w: 1.3, n: 300 });
      pMag.fn(fq => tfSingle(TAU * fq * H / vs, xi, 0).f2, { c: 'blue', w: 2.2, n: 500 });

      const fModes = [f1, 3 * f1, 5 * f1];
      const pModes = fModes.map(fq => tfSingle(TAU * fq * H / vs, xi, 0).f2);
      pMag.markers(fModes, pModes, { c: 'red', fill: true, r: 4.5, over: true });

      pPh.setLim([0, fmax], [-540, 20]);
      pPh.clear();
      for (const y of [-90, -270, -450]) pPh.hline(y, { c: 'gray', dash: '3 3' });
      const fs = linspace(0, fmax, 400);
      const phs = [];
      let prev = 0, wrap = 0;
      for (const fq of fs) {
        const t = tfSingle(TAU * fq * H / vs, xi, 0);
        let a = Math.atan2(t.f2i, t.f2r);
        if (a - prev > Math.PI) wrap -= TAU;
        if (a - prev < -Math.PI) wrap += TAU;
        prev = a;
        phs.push(deg(a + wrap));
      }
      pPh.line(fs, phs, { c: 'blue', w: 2.2 });
      pPh.markers(fModes, [-90, -270, -450], { c: 'red', fill: true, r: 4.5, over: true });

      st({
        f1: f(f1, 3),
        T1: f(1 / f1, 2),
        p1: f(pModes[0], 2),
        p2: f(pModes[1], 2),
        p3: f(pModes[2], 2)
      });
    }
    draw();
  })();

  /* =========================================================
     4. #w-fft：FFT → 乘上轉換函數 → 逆 FFT (6.4)
     ========================================================= */
  (() => {
    const body = SD.widget('w-fft'); if (!body) return;
    const mot = window.CH6_MOTION;
    if (!mot) return;

    const N = mot.a.length, dt = mot.dt, df = 1 / (N * dt);
    const ain = Float64Array.from(mot.a, v => v * mot.unit);
    const baseRe = Float64Array.from(ain), baseIm = new Float64Array(N);
    fft(baseRe, baseIm, -1);

    const ctl = SD.controls(body, [
      { k: 'H', label: '土層厚度 H (m)', min: 3, max: 200, step: 1, value: 165, fmt: v => v.toFixed(0) },
      { k: 'vs', label: '剪力波速 Vₛ (m/s)', min: 100, max: 600, step: 10, value: 457, fmt: v => v.toFixed(0) },
      { k: 'xi', label: '阻尼比 ξ (%)', min: 2, max: 20, step: 1, value: 5, fmt: v => v.toFixed(0) + ' %' },
      { k: 'az', label: '阻抗比 αz (0=剛性岩盤)', min: 0, max: 0.6, step: 0.05, value: 0, fmt: v => v.toFixed(2) },
      { k: 'noPhase', label: '只乘 |F|、丟掉相位 ∠F（示範失去相位後的波形失真）', type: 'check', value: false }
    ], draw);

    SD.buttons(body, '講義 3.2 節案例', [
      { label: '10 ft 薄層 (H = 3 m, Vs = 457 → f₁ ≈ 38 Hz)', on: () => ctl.set({ H: 3, vs: 457, xi: 5, az: 0, noPhase: false }) },
      { label: '540 ft 厚層・剛性岩盤 (H = 165 m, f₁ ≈ 0.69 Hz)', on: () => ctl.set({ H: 165, vs: 457, xi: 5, az: 0, noPhase: false }) },
      { label: '540 ft 厚層・彈性岩盤 (αz = 0.25)', on: () => ctl.set({ H: 165, vs: 457, xi: 5, az: 0.25, noPhase: false }) },
      { label: '丟掉相位 ∠F（錯誤示範）', on: () => ctl.set({ H: 165, vs: 457, xi: 5, az: 0, noPhase: true }) }
    ]);

    const grid = SD.div(body, 'plot-row plot-grid');
    const pSpec = new Plot(grid, { h: 270, xlim: [0, 10], ylim: [0, 1], xlabel: '頻率 f (Hz)', ylabel: '正規化傅立葉振幅譜 & 轉換函數' });
    pSpec.legend([
      { label: '輸入振幅譜 |A_{in}(f)|', c: 'gray', w: 1.4 },
      { label: '地表振幅譜 |A_{out}(f)|', c: 'blue', w: 1.8 },
      { label: '轉換函數 |F(f)| (縮放顯示)', c: 'red', dash: true, w: 1.8 }
    ]);
    const pTime = new Plot(grid, { h: 270, xlim: [0, 41], ylim: [-0.6, 0.6], xlabel: '時間 t (s)', ylabel: '加速度 a(t) (g)' });
    pTime.legend([
      { label: '輸入運動 a_{in}(t)', c: 'gray', w: 1.2 },
      { label: '地表反應 a_{out}(t)', c: 'blue', w: 1.5 }
    ]);
    const st = SD.stats(body, [
      { k: 'f1', label: '土層基本頻率 f₁ (Hz)', hl: true },
      { k: 'pgaIn', label: '輸入 PGA_{in} (g)' },
      { k: 'pgaOut', label: '地表 PGA_{out} (g)', hl: true },
      { k: 'amp', label: 'PGA 放大比', hl: true },
      { k: 'fmax', label: '轉換函數峰值 |F|_{max}' },
      { k: 'phMode', label: '相位處理' }
    ]);

    function draw() {
      const H = ctl.H, vs = ctl.vs, xi = ctl.xi / 100, az = ctl.az;
      const re = new Float64Array(N), im = new Float64Array(N);
      re[0] = baseRe[0]; im[0] = 0;
      let fMaxVal = 1;
      for (let k = 1; k <= N / 2; k++) {
        const fq = k * df;
        const tf = tfSingle(TAU * fq * H / vs, xi, az);
        const fr = az === 0 ? tf.f2r : tf.f3r;
        const fi = az === 0 ? tf.f2i : tf.f3i;
        const fm = Math.hypot(fr, fi);
        if (fm > fMaxVal) fMaxVal = fm;
        if (ctl.noPhase) {
          re[k] = baseRe[k] * fm;
          im[k] = baseIm[k] * fm;
        } else {
          re[k] = baseRe[k] * fr - baseIm[k] * fi;
          im[k] = baseRe[k] * fi + baseIm[k] * fr;
        }
        if (k < N / 2) {
          re[N - k] = re[k];
          im[N - k] = -im[k];
        }
      }
      const outRe = Float64Array.from(re), outIm = Float64Array.from(im);
      fft(outRe, outIm, 1);
      const aout = new Float64Array(N);
      let pgaOut = 0;
      for (let i = 0; i < N; i++) {
        aout[i] = outRe[i] / N;
        if (Math.abs(aout[i]) > pgaOut) pgaOut = Math.abs(aout[i]);
      }

      const kMax = Math.floor(10 / df);
      const fs = [], sIn = [], sOut = [];
      let sMax = 1e-9;
      for (let k = 1; k <= kMax; k += 2) {
        const mIn = Math.hypot(baseRe[k], baseIm[k]) * dt;
        const mOut = Math.hypot(re[k], im[k]) * dt;
        if (mOut > sMax) sMax = mOut;
        if (mIn > sMax) sMax = mIn;
        fs.push(k * df); sIn.push(mIn); sOut.push(mOut);
      }
      pSpec.clear();
      pSpec.line(fs, sIn.map(v => v / sMax), { c: 'gray', w: 1.3 });
      pSpec.line(fs, sOut.map(v => v / sMax), { c: 'blue', w: 1.7 });
      pSpec.fn(fq => {
        const tf = tfSingle(TAU * fq * H / vs, xi, az);
        return (az === 0 ? tf.f2 : tf.f3) / fMaxVal;
      }, { c: 'red', dash: '4 3', w: 1.8, n: 300 });

      const ym = Math.max(0.3, Math.ceil(Math.max(mot.pga0, pgaOut) * 10) / 10);
      pTime.setLim([0, 41], [-ym, ym]);
      pTime.clear();
      const ts = [], yIn = [], yOut = [];
      for (let i = 0; i < N; i += 4) {
        ts.push(i * dt); yIn.push(ain[i]); yOut.push(aout[i]);
      }
      pTime.line(ts, yIn, { c: 'gray', w: 1.1 });
      pTime.line(ts, yOut, { c: ctl.noPhase ? 'red' : 'blue', w: 1.4 });

      st({
        f1: f(vs / (4 * H), 2),
        pgaIn: f(mot.pga0, 3),
        pgaOut: f(pgaOut, 3),
        amp: f(pgaOut / mot.pga0, 2) + '×',
        fmax: f(fMaxVal, 2),
        phMode: ctl.noPhase ? '丟棄相位（失真）' : '完整複數 F*(ω)'
      });
    }
    draw();
  })();

  /* =========================================================
     5. #w-f2：|F₂| 的無因次解圖（橫軸 kH） (6.5)
     ========================================================= */
  (() => {
    const body = SD.widget('w-f2'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'H', label: '土層厚度 H (m)', min: 5, max: 120, step: 5, value: 25, fmt: v => v.toFixed(0) },
      { k: 'vs', label: '剪力波速 V_s (m/s)', min: 80, max: 600, step: 10, value: 200, fmt: v => v.toFixed(0) },
      { k: 'freq', label: '查詢頻率 f (Hz)', min: 0.2, max: 12, step: 0.1, value: 2.0, fmt: v => v.toFixed(1) },
      { k: 'xi', label: '阻尼比 ξ (%)', min: 1, max: 25, step: 1, value: 5, fmt: v => v.toFixed(0) + ' %' }
    ], draw);
    SD.buttons(body, '快速定位', [
      { label: '一階共振 (kH = π/2)', on: () => ctl.set({ freq: ctl.vs / (4 * ctl.H) }) },
      { label: '反共振波谷 (kH = π)', on: () => ctl.set({ freq: 2 * ctl.vs / (4 * ctl.H) }) },
      { label: '二階共振 (kH = 3π/2)', on: () => ctl.set({ freq: 3 * ctl.vs / (4 * ctl.H) }) }
    ]);
    const p = new Plot(body, { h: 270, xlim: [0, 10], ylim: [0, 16], xlabel: '無因次波數 kH = ωH / V_s (rad)', ylabel: '放大倍數 |F₂(kH, ξ)|' });
    p.legend([
      { label: '目前阻尼比 ξ 曲線', c: 'blue', w: 2.4 },
      { label: '背景參考曲線 (ξ = 2%, 5%, 10%, 20%)', c: 'gray', dash: true },
      { label: '目前查詢點 (kH, |F₂|)', c: 'red', line: false, marker: 'circle', fill: true }
    ]);
    const st = SD.stats(body, [
      { k: 'kh', label: '無因次波數 kH (rad)', hl: true },
      { k: 'khp', label: 'kH / π' },
      { k: 'f2', label: '精確解 |F₂|', hl: true },
      { k: 'f2ap', label: '簡化式 1/√(cos²kH + (ξkH)²)' },
      { k: 'err', label: '簡化式誤差 (%)' },
      { k: 'f1', label: '對應 f₁ = V_s/(4H) (Hz)' }
    ]);

    function draw() {
      const H = ctl.H, vs = ctl.vs, fq = ctl.freq, xi = ctl.xi / 100;
      const kh = Math.min(10, TAU * fq * H / vs);
      const exact = tfSingle(kh, xi, 0).f2;
      const approx = 1 / Math.hypot(Math.cos(kh), xi * kh);

      p.clear();
      for (const xr of [0.02, 0.05, 0.10, 0.20]) {
        p.fn(x => Math.min(15.8, tfSingle(x, xr, 0).f2), { c: 'gray', dash: '3 3', w: 1.1, n: 350 });
      }
      p.fn(x => Math.min(15.8, tfSingle(x, xi, 0).f2), { c: 'blue', w: 2.4, n: 400 });
      p.vline(kh, { c: 'red', dash: '3 3', w: 1.2 });
      p.markers([kh], [Math.min(15.8, exact)], { c: 'red', fill: true, r: 5, over: true });

      st({
        kh: f(kh, 3),
        khp: f(kh / Math.PI, 2),
        f2: f(exact, 3),
        f2ap: f(approx, 3),
        err: f(Math.abs(approx - exact) / exact * 100, 2) + ' %',
        f1: f(vs / (4 * H), 2)
      });
    }
    draw();
  })();

  /* =========================================================
     6. #w-modes：土層振態與任一頻率下的變形形狀 (6.6)
     ========================================================= */
  (() => {
    const body = SD.widget('w-modes'); if (!body) return;
    let playing = true, ph = 0;
    const ctl = SD.controls(body, [
      { k: 'r', label: '激振頻率比 f / f₁（= 2kH / π）', min: 0.1, max: 7.5, step: 0.05, value: 1.0, fmt: v => v.toFixed(2) },
      { k: 'xi', label: '阻尼比 ξ (%)', min: 2, max: 20, step: 1, value: 5, fmt: v => v.toFixed(0) + ' %' }
    ], render);
    SD.buttons(body, '自然振態與反共振', [
      { label: 'n = 1 (f/f₁ = 1)', on: () => ctl.set({ r: 1.0 }) },
      { label: '反共振 (f/f₁ = 2)', on: () => ctl.set({ r: 2.0 }) },
      { label: 'n = 2 (f/f₁ = 3)', on: () => ctl.set({ r: 3.0 }) },
      { label: 'n = 3 (f/f₁ = 5)', on: () => ctl.set({ r: 5.0 }) },
      { label: 'n = 4 (f/f₁ = 7)', on: () => ctl.set({ r: 7.0 }) },
      { label: '❚❚ 暫停', primary: true, on: b => { playing = !playing; b.textContent = playing ? '❚❚ 暫停' : '▶ 播放'; if (playing) requestAnimationFrame(loop); } }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pM = new Plot(grid, { h: 270, xlim: [-1.15, 1.15], ylim: [-1, 0], yfmt: v => Math.abs(v).toFixed(1), xlabel: '正規化振態 φ_n(z) = cos((2n−1)πz / 2H)', ylabel: '正規化深度 z / H' });
    pM.legend([
      { label: '第 1 振態 (n=1)', c: 'blue', w: 2.2 },
      { label: '第 2 振態 (n=2)', c: 'red', w: 1.8 },
      { label: '第 3 振態 (n=3)', c: 'ink', w: 1.6 },
      { label: '第 4 振態 (n=4)', c: 'gray', dash: true }
    ]);
    const pD = new Plot(grid, { h: 270, xlim: [-14, 14], ylim: [-1, 0], yfmt: v => Math.abs(v).toFixed(1), xlabel: '穩態位移比 u(z, t) / u_{rock}', ylabel: '正規化深度 z / H' });
    pD.legend([
      { label: '瞬時變形剖面 u(z, t)/u_{rock}', c: 'blue', w: 2.4 },
      { label: '振幅包絡線 ±|F₂(z)|', c: 'gray', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'r', label: '頻率比 f / f₁', hl: true },
      { k: 'kh', label: 'kH / π' },
      { k: 'surf', label: '地表放大 |u(0)/u_{rock}|', hl: true },
      { k: 'maxIn', label: '土柱內最大放大' },
      { k: 'nodes', label: '內部節點數 (擬零點)' }
    ]);

    const zs = linspace(0, 1, 140), nzs = zs.map(z => -z);
    pM.vline(0, { c: 'gray', dash: null, w: 0.7 });
    const cols = ['blue', 'red', 'ink', 'gray'];
    [1, 2, 3, 4].forEach((n, idx) => {
      const xs = zs.map(z => Math.cos((2 * n - 1) * Math.PI * 0.5 * z));
      pM.line(xs, nzs, { c: cols[idx], w: idx === 0 ? 2.2 : 1.6, dash: idx === 3 ? '4 3' : null });
    });

    function render() {
      const r = ctl.r, xi = ctl.xi / 100, kh = r * Math.PI / 2;
      const [kr, ki] = kStarH(kh, xi);
      const [cbR, cbI] = cCos(kr, ki);
      const env = [], inst = [];
      let maxIn = 0;
      for (const z of zs) {
        const [czR, czI] = cCos(kr * z, ki * z);
        const [fr, fi] = cDiv(czR, czI, cbR, cbI);
        const m = Math.hypot(fr, fi);
        if (m > maxIn) maxIn = m;
        env.push(m);
        inst.push(fr * Math.cos(ph) - fi * Math.sin(ph));
      }
      const xm = Math.max(3, Math.ceil(maxIn * 1.15));
      pD.setLim([-xm, xm], [-1, 0]);
      pD.clear();
      pD.vline(0, { c: 'gray', dash: null, w: 0.7 });
      pD.line(env, nzs, { c: 'gray', dash: '4 3', w: 1.4 });
      pD.line(env.map(v => -v), nzs, { c: 'gray', dash: '4 3', w: 1.4 });
      pD.line(inst, nzs, { c: 'blue', w: 2.4 });
      pD.markers([inst[0], inst[inst.length - 1]], [0, -1], { c: 'red', fill: true, r: 4.5, over: true });

      const nNodes = Math.floor((r + 1) / 2) - 1;
      st({
        r: f(r, 2),
        kh: f(r / 2, 2),
        surf: f(env[0], 2),
        maxIn: f(maxIn, 2),
        nodes: String(Math.max(0, nNodes))
      });
    }
    let last = null;
    function loop(ts) {
      if (!playing) { last = null; return; }
      if (last != null) ph = (ph + (ts - last) * 0.004) % TAU;
      last = ts;
      render();
      requestAnimationFrame(loop);
    }
    render();
    requestAnimationFrame(loop);
  })();

  /* =========================================================
     7. #w-f3：阻抗比與輻射阻尼 F₃ (6.7)
     ========================================================= */
  (() => {
    const body = SD.widget('w-f3'); if (!body) return;
    let showErrata = false;
    const ctl = SD.controls(body, [
      { k: 'az', label: '阻抗比 α_z = ρ_s V_s / (ρ_r V_r)', min: 0.0, max: 1.0, step: 0.01, value: 0.25, fmt: v => v.toFixed(2) },
      { k: 'xi', label: '土層材料阻尼比 ξ_s (%)', min: 0.0, max: 15.0, step: 0.5, value: 5.0, fmt: v => v.toFixed(1) + ' %' }
    ], () => { draw(); });

    SD.buttons(body, '情境與講義對照', [
      { label: '典型軟土 / 岩盤 (α_z = 0.25, ξ_s = 5%)', on: () => { showErrata = false; ctl.set({ az: 0.25, xi: 5.0 }); } },
      { label: '純輻射阻尼 (ξ_s = 0, α_z = 0.25 → 峰高 1/α_z = 4)', on: () => { showErrata = false; ctl.set({ az: 0.25, xi: 0.0 }); } },
      { label: '講義 2.7 圖（ξ_s = 0，對照 α_z 與 α_z² 標籤勘誤）', primary: true, on: () => { showErrata = true; ctl.set({ az: 0.316, xi: 0.0 }); } },
      { label: '均勻半空間無界面 (α_z = 1.0 → |F₃| ≡ 1)', on: () => { showErrata = false; ctl.set({ az: 1.0, xi: 0.0 }); } }
    ]);

    const p = new Plot(body, { h: 280, xlim: [0, 10], ylim: [0, 10], xlabel: '無因次波數 kH = ωH / V_s (rad)', ylabel: '轉換函數振幅 |F₃(kH)|' });
    p.legend([
      { label: '彈性岩盤 |F₃(kH)| (相對露頭)', c: 'blue', w: 2.4 },
      { label: '剛性岩盤 |F₂(kH)| (α_z = 0 對照)', c: 'gray', dash: true },
      { label: '講義 2.7 圖四條參考線 (ξ_s = 0)', c: 'red', dash: '2 3' }
    ]);
    const st = SD.stats(body, [
      { k: 'az', label: '真實阻抗比 α_z', hl: true },
      { k: 'az2', label: '講義圖誤標值 α_z²', hl: true },
      { k: 'xirad', label: '等效輻射阻尼 ≈ 2α_z/π' },
      { k: 'p1', label: '第 1 峰 |F₃(π/2)|', hl: true },
      { k: 'p2', label: '第 2 峰 |F₃(3π/2)|' },
      { k: 'p0', label: 'ξ_s=0 理論峰高 1/α_z' }
    ]);

    function draw() {
      const az = ctl.az, xi = ctl.xi / 100;
      const p1 = tfSingle(Math.PI / 2, xi, az).f3;
      const p2 = tfSingle(3 * Math.PI / 2, xi, az).f3;
      p.clear();

      if (showErrata) {
        for (const [a2, pk] of [[0.5, 1.41], [0.1, 3.16], [0.05, 4.47], [0.02, 7.07]]) {
          const aTrue = Math.sqrt(a2);
          p.fn(x => tfSingle(x, 0, aTrue).f3, { c: 'red', dash: '2 3', w: 1.3, n: 350 });
          p.text(Math.PI / 2 + 0.15, Math.min(9.5, pk), `峰 ${pk.toFixed(2)} (真 α_z=${aTrue.toFixed(2)}, 圖標 α_z²=${a2})`, { cls: 'txt-red' });
        }
      }
      if (xi > 0) {
        p.fn(x => Math.min(9.9, tfSingle(x, xi, 0).f2), { c: 'gray', dash: '5 4', w: 1.5, n: 350 });
      }
      p.fn(x => Math.min(9.9, tfSingle(x, xi, az).f3), { c: 'blue', w: 2.4, n: 400 });
      p.markers([Math.PI / 2, 3 * Math.PI / 2], [Math.min(9.8, p1), Math.min(9.8, p2)], { c: 'blue', fill: true, r: 4.5, over: true });

      st({
        az: f(az, 3),
        az2: f(az * az, 3),
        xirad: f((2 * az / Math.PI) * 100, 1) + ' %',
        p1: p1 > 99 ? '∞' : f(p1, 2),
        p2: p2 > 99 ? '∞' : f(p2, 2),
        p0: az > 0 ? f(1 / az, 2) : '∞'
      });
    }
    draw();
  })();

  /* =========================================================
     8. #w-layers：層狀剖面的轉換函數（可編輯） (6.8)
     ========================================================= */
  (() => {
    const body = SD.widget('w-layers'); if (!body) return;
    let layers = [
      { h: 6, vs: 140, xi: 0.05, rho: 1.75 },
      { h: 12, vs: 240, xi: 0.05, rho: 1.90 },
      { h: 18, vs: 360, xi: 0.04, rho: 2.05 },
      { h: 0, vs: 1000, xi: 0.01, rho: 2.40 }
    ];
    let showEq = true;

    SD.buttons(body, '預設層狀剖面', [
      { label: '正常遞增三層剖面 (36 m)', on: () => { layers = [{ h: 6, vs: 140, xi: 0.05, rho: 1.75 }, { h: 12, vs: 240, xi: 0.05, rho: 1.90 }, { h: 18, vs: 360, xi: 0.04, rho: 2.05 }, { h: 0, vs: 1000, xi: 0.01, rho: 2.40 }]; rebuild(); } },
      { label: '講義 3.2 節三層剖面 (540 ft ≈ 165 m)', on: () => { layers = [{ h: 30, vs: 240, xi: 0.05, rho: 1.90 }, { h: 60, vs: 420, xi: 0.05, rho: 2.00 }, { h: 75, vs: 600, xi: 0.04, rho: 2.15 }, { h: 0, vs: 1500, xi: 0.01, rho: 2.50 }]; rebuild(); } },
      { label: '含軟弱夾層剖面 (速度逆轉)', on: () => { layers = [{ h: 6, vs: 280, xi: 0.04, rho: 1.95 }, { h: 10, vs: 120, xi: 0.06, rho: 1.70 }, { h: 14, vs: 340, xi: 0.04, rho: 2.05 }, { h: 0, vs: 900, xi: 0.01, rho: 2.35 }]; rebuild(); } },
      { label: '+ 新增土層', on: () => { layers.splice(layers.length - 1, 0, { h: 8, vs: 250, xi: 0.05, rho: 1.90 }); rebuild(); } },
      { label: '切換「同走時等效單一層」比較', on: () => { showEq = !showEq; rebuild(); } }
    ]);

    const tblWrap = SD.div(body, '');
    const grid = SD.div(body, 'plot-row plot-grid');
    const pProf = new Plot(grid, { h: 270, xlim: [0, 800], ylim: [-50, 0], yfmt: v => String(Math.round(Math.abs(v))), xlabel: '剪力波速 V_s (m/s)', ylabel: '深度 z (m)' });
    pProf.legend([
      { label: '分層波速 V_s(z)', c: 'ink', w: 2.2 },
      { label: '同走時等效單一層 V_{s,eq}', c: 'gray', dash: true }
    ]);
    const pTF = new Plot(grid, { h: 270, xlim: [0, 10], ylim: [0, 10], xlabel: '頻率 f (Hz)', ylabel: '轉換函數振幅 |F(f)|' });
    pTF.legend([
      { label: '多層遞迴 |F₃(f)| (相對露頭)', c: 'blue', w: 2.3 },
      { label: '多層遞迴 |F₂(f)| (相對岩盤頂)', c: 'red', dash: '4 3', w: 1.6 },
      { label: '同走時均勻單層 |F_{3,eq}(f)|', c: 'gray', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'H', label: '土柱總厚度 H (m)' },
      { k: 'ts', label: '單程走時 t_s = Σ(h_m/V_{sm}) (s)' },
      { k: 'veq', label: '走時平均波速 V_{s,eq} (m/s)', hl: true },
      { k: 'f1eq', label: '均勻近似 1/(4t_s) (Hz)' },
      { k: 'fpk', label: '實際第一峰頻率 f_{peak} (Hz)', hl: true },
      { k: 'amp', label: '第一峰放大倍率 |F₃|_{max}', hl: true }
    ]);

    function rebuild() {
      tblWrap.innerHTML = `<table class="edit"><thead><tr>
        <th>層號</th><th>厚度 h_m (m)</th><th>波速 V_sm (m/s)</th><th>阻尼比 ξ_m (%)</th><th>密度 ρ_m (t/m³)</th><th></th>
      </tr></thead><tbody>${layers.map((l, i) => {
        const hs = i === layers.length - 1;
        return `<tr>
          <td>${hs ? '岩盤半空間' : '土層 ' + (i + 1)}</td>
          <td>${hs ? '∞' : `<input type="number" data-i="${i}" data-k="h" min="1" max="150" step="1" value="${l.h}">`}</td>
          <td><input type="number" data-i="${i}" data-k="vs" min="50" max="3000" step="10" value="${l.vs}"></td>
          <td><input type="number" data-i="${i}" data-k="xi" min="0" max="25" step="0.5" value="${+(l.xi * 100).toFixed(1)}"></td>
          <td><input type="number" data-i="${i}" data-k="rho" min="1.3" max="2.8" step="0.05" value="${l.rho}"></td>
          <td>${(!hs && layers.length > 2) ? `<button class="btn" data-del="${i}" style="padding:1px 7px">×</button>` : ''}</td>
        </tr>`;
      }).join('')}</tbody></table>`;

      tblWrap.querySelectorAll('input').forEach(inp => {
        inp.addEventListener('change', () => {
          const i = +inp.dataset.i, k = inp.dataset.k, v = parseFloat(inp.value);
          if (isFinite(v) && v >= 0) {
            layers[i][k] = k === 'xi' ? v / 100 : v;
            drawOnly();
          }
        });
      });
      tblWrap.querySelectorAll('button[data-del]').forEach(btn => {
        btn.addEventListener('click', () => { layers.splice(+btn.dataset.del, 1); rebuild(); });
      });
      drawOnly();
    }

    function drawOnly() {
      const soils = layers.slice(0, -1), rock = layers[layers.length - 1];
      const H = soils.reduce((s, l) => s + l.h, 0);
      const ts = soils.reduce((s, l) => s + l.h / l.vs, 0);
      const veq = H / ts;
      const rhoEq = soils.reduce((s, l) => s + l.h * l.rho, 0) / H;
      const xiEq = soils.reduce((s, l) => s + (l.h / l.vs) * l.xi, 0) / ts;
      const eqLayers = [{ h: H, vs: veq, xi: xiEq, rho: rhoEq }, rock];

      const zmax = Math.ceil(H * 1.2 / 5) * 5;
      const vmax = Math.ceil(Math.max(500, ...soils.map(l => l.vs * 1.25)) / 50) * 50;
      pProf.setLim([0, vmax], [-zmax, 0]);
      pProf.clear();
      const xs = [], nzs = [];
      let z = 0;
      for (const l of soils) { xs.push(l.vs, l.vs); nzs.push(-z, -(z + l.h)); z += l.h; }
      xs.push(Math.min(vmax, rock.vs), Math.min(vmax, rock.vs)); nzs.push(-z, -zmax);
      pProf.line(xs, nzs, { c: 'ink', w: 2.2 });
      if (showEq) pProf.vline(veq, { c: 'gray', dash: '5 4', w: 1.6 });

      const fmax = Math.max(4, Math.min(15, Math.ceil(4.5 / ts)));
      const fs = linspace(0.02, fmax, 360);
      const f3s = [], f2s = [], f3eqs = [];
      let fpk = 0, amax = 0;
      for (const fq of fs) {
        const r = tfMulti(fq, layers);
        f3s.push(r.f3); f2s.push(r.f2);
        if (r.f3 > amax && fq < 2.5 / ts) { amax = r.f3; fpk = fq; }
        if (showEq) f3eqs.push(tfMulti(fq, eqLayers).f3);
      }
      const ym = Math.min(18, Math.max(6, Math.ceil(Math.max(...f3s) * 1.15)));
      pTF.setLim([0, fmax], [0, ym]);
      pTF.clear();
      if (showEq) pTF.line(fs, f3eqs, { c: 'gray', dash: '4 3', w: 1.6 });
      pTF.line(fs, f2s.map(v => Math.min(ym * 1.05, v)), { c: 'red', dash: '4 3', w: 1.5 });
      pTF.line(fs, f3s, { c: 'blue', w: 2.3 });
      pTF.markers([fpk], [amax], { c: 'blue', fill: true, r: 4.5, over: true });

      st({
        H: f(H, 1),
        ts: f(ts, 3),
        veq: f(veq, 1),
        f1eq: f(1 / (4 * ts), 2),
        fpk: f(fpk, 2),
        amp: f(amax, 2)
      });
    }
    rebuild();
  })();

  /* =========================================================
     9. #w-sites：軟弱盆地 vs 堅硬場址（課堂活動 1） (6.9)
     ========================================================= */
  (() => {
    const body = SD.widget('w-sites'); if (!body) return;
    const mot = window.CH6_MOTION;
    if (!mot) return;

    const N = mot.a.length, dt = mot.dt, df = 1 / (N * dt);
    const ain = Float64Array.from(mot.a, v => v * mot.unit);
    const baseRe = Float64Array.from(ain), baseIm = new Float64Array(N);
    fft(baseRe, baseIm, -1);

    const ctl = SD.controls(body, [
      { k: 'hA', label: '軟土盆地 A 厚度 H_A (m)', min: 15, max: 80, step: 5, value: 40, fmt: v => v.toFixed(0) },
      { k: 'vsA', label: '軟土盆地 A 波速 V_sA (m/s)', min: 70, max: 220, step: 10, value: 120, fmt: v => v.toFixed(0) },
      { k: 'hB', label: '堅硬場址 B 厚度 H_B (m)', min: 5, max: 40, step: 5, value: 15, fmt: v => v.toFixed(0) },
      { k: 'vsB', label: '堅硬場址 B 波速 V_sB (m/s)', min: 250, max: 600, step: 10, value: 380, fmt: v => v.toFixed(0) },
      { k: 'nStory', label: '目標建物樓層數 N (T_b≈0.1N)', min: 1, max: 24, step: 1, value: 12, fmt: v => `${v} 層` }
    ], draw);

    SD.buttons(body, '課堂活動情境', [
      { label: '12 層大樓 (T_b = 1.2 s，軟土盆地強烈共振)', on: () => ctl.set({ hA: 40, vsA: 120, hB: 15, vsB: 380, nStory: 12 }) },
      { label: '2 層低矮建物 (T_b = 0.2 s，堅硬場址較不利)', on: () => ctl.set({ hA: 40, vsA: 120, hB: 15, vsB: 380, nStory: 2 }) },
      { label: '墨西哥城湖相軟泥 (H_A = 45 m, VsA = 75 m/s → T₁ = 2.4 s)', on: () => ctl.set({ hA: 45, vsA: 75, hB: 15, vsB: 400, nStory: 20 }) }
    ]);

    const grid = SD.div(body, 'plot-row plot-grid');
    const pTime = new Plot(grid, { h: 270, xlim: [5, 35], ylim: [-0.6, 0.6], xlabel: '時間 t (s)（擷取強震段）', ylabel: '地表加速度 a(t) (g)' });
    pTime.legend([
      { label: '軟弱盆地 A 地表歷時', c: 'red', w: 1.5 },
      { label: '堅硬場址 B 地表歷時', c: 'blue', w: 1.3 }
    ]);
    const pSa = new Plot(grid, { h: 270, xlim: [0.05, 2.5], ylim: [0, 1.8], xlabel: '結構自然週期 T (s)', ylabel: '5% 阻尼擬加速度反應譜 S_a(T) (g)' });
    pSa.legend([
      { label: '軟弱盆地 A 反應譜', c: 'red', w: 2.2 },
      { label: '堅硬場址 B 反應譜', c: 'blue', w: 2.2 },
      { label: '岩盤露頭輸入反應譜', c: 'gray', dash: true },
      { label: '目標建物週期 T_b = 0.1 N', c: 'ink', dash: '2 3' }
    ]);
    const st = SD.stats(body, [
      { k: 'TA', label: '場址 A 基本週期 T₁A (s)', hl: true },
      { k: 'TB', label: '場址 B 基本週期 T₁B (s)' },
      { k: 'pgaA', label: '場址 A 地表 PGA (g)' },
      { k: 'pgaB', label: '場址 B 地表 PGA (g)' },
      { k: 'saA', label: '建物在 A 的 S_a(T_b) (g)', hl: true },
      { k: 'saB', label: '建物在 B 的 S_a(T_b) (g)', hl: true }
    ]);

    function synthMotion(H, vs, xi, az) {
      const re = new Float64Array(N), im = new Float64Array(N);
      re[0] = baseRe[0];
      for (let k = 1; k <= N / 2; k++) {
        const fq = k * df;
        const tf = tfSingle(TAU * fq * H / vs, xi, az);
        re[k] = baseRe[k] * tf.f3r - baseIm[k] * tf.f3i;
        im[k] = baseRe[k] * tf.f3i + baseIm[k] * tf.f3r;
        if (k < N / 2) { re[N - k] = re[k]; im[N - k] = -im[k]; }
      }
      fft(re, im, 1);
      const out = new Float64Array(N);
      let pga = 0;
      for (let i = 0; i < N; i++) {
        out[i] = re[i] / N;
        if (Math.abs(out[i]) > pga) pga = Math.abs(out[i]);
      }
      return { out, pga };
    }

    function respSpec(acc, Ts) {
      const step = 2, dtSub = dt * step;
      const ag = new Float64Array(N / step);
      for (let i = 0; i < ag.length; i++) ag[i] = acc[i * step];
      return Ts.map(T => {
        const w = TAU / T;
        return w * w * SDOF.peakDisp(ag, dtSub, T, 0.05);
      });
    }

    const Ts = [0.06, 0.1, 0.14, 0.18, 0.22, 0.28, 0.35, 0.42, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.45, 1.6, 1.8, 2.0, 2.25, 2.5];
    const saRock = respSpec(ain, Ts);

    function draw() {
      const TA = 4 * ctl.hA / ctl.vsA, TB = 4 * ctl.hB / ctl.vsB;
      const azA = (1.75 * ctl.vsA) / (2.4 * 1000), azB = (2.05 * ctl.vsB) / (2.4 * 1000);
      const mA = synthMotion(ctl.hA, ctl.vsA, 0.05, azA);
      const mB = synthMotion(ctl.hB, ctl.vsB, 0.04, azB);

      const ym = Math.max(0.4, Math.ceil(Math.max(mA.pga, mB.pga) * 10) / 10);
      pTime.setLim([5, 35], [-ym, ym]);
      pTime.clear();
      const ts = [], yA = [], yB = [];
      for (let i = 500; i <= 3500; i += 3) {
        ts.push(i * dt); yA.push(mA.out[i]); yB.push(mB.out[i]);
      }
      pTime.line(ts, yB, { c: 'blue', w: 1.2 });
      pTime.line(ts, yA, { c: 'red', w: 1.4 });

      const saA = respSpec(mA.out, Ts);
      const saB = respSpec(mB.out, Ts);
      const Tb = 0.1 * ctl.nStory;
      const saAtA = (TAU / Tb) ** 2 * SDOF.peakDisp(mA.out, dt, Tb, 0.05);
      const saAtB = (TAU / Tb) ** 2 * SDOF.peakDisp(mB.out, dt, Tb, 0.05);

      const saMax = Math.max(1.2, Math.ceil(Math.max(...saA, ...saB, saAtA, saAtB) * 1.15 * 5) / 5);
      pSa.setLim([0.05, 2.5], [0, saMax]);
      pSa.clear();
      pSa.line(Ts, saRock, { c: 'gray', dash: '4 3', w: 1.4 });
      pSa.line(Ts, saB, { c: 'blue', w: 2.2 });
      pSa.line(Ts, saA, { c: 'red', w: 2.2 });
      pSa.vline(Tb, { c: 'ink', dash: '3 3', w: 1.5 });
      pSa.markers([Tb, Tb], [saAtA, saAtB], { c: 'red', fill: true, r: 5, over: true });

      st({
        TA: f(TA, 2),
        TB: f(TB, 2),
        pgaA: f(mA.pga, 3),
        pgaB: f(mB.pga, 3),
        saA: f(saAtA, 2),
        saB: f(saAtB, 2)
      });
    }
    draw();
  })();

  /* =========================================================
     10. #w-eql：多層等值線性分析（講義 2.21 剖面） (6.10)
     ========================================================= */
  (() => {
    const body = SD.widget('w-eql'); if (!body) return;
    const mot = window.CH6_MOTION;
    if (!mot) return;

    const N = mot.a.length, dt = mot.dt, df = 1 / (N * dt);
    const baseRe = Float64Array.from(mot.a, v => v * mot.unit / mot.pga0);
    const baseIm = new Float64Array(N);
    fft(baseRe, baseIm, -1);

    const ctl = SD.controls(body, [
      { k: 'pga', label: '輸入露頭峰值加速度 PGA_in (g)', min: 0.02, max: 0.60, step: 0.02, value: 0.25, fmt: v => v.toFixed(2) + ' g' },
      { k: 'M', label: '地震規模 M（決定 R_γ = (M−1)/10）', min: 5.0, max: 8.0, step: 0.2, value: 6.5, fmt: v => v.toFixed(1) }
    ], draw);

    SD.buttons(body, '強度比較', [
      { label: '微弱地震 (PGA = 0.03 g，接近線性)', on: () => ctl.set({ pga: 0.03, M: 5.5 }) },
      { label: '中強地震 (PGA = 0.25 g, M = 6.5)', on: () => ctl.set({ pga: 0.25, M: 6.5 }) },
      { label: '強烈地震 (PGA = 0.50 g, M = 7.5，顯著軟化)', on: () => ctl.set({ pga: 0.50, M: 7.5 }) }
    ]);

    const grid = SD.div(body, 'plot-row plot-grid');
    const pTF = new Plot(grid, { h: 270, xlim: [0, 8], ylim: [0, 10], xlabel: '頻率 f (Hz)', ylabel: '轉換函數 |F₃(f)|' });
    pTF.legend([
      { label: '初始線性 (第 1 輪 G_{max})', c: 'gray', dash: true, w: 1.6 },
      { label: '第 2 輪迭代', c: 'red', dash: '2 3', w: 1.5 },
      { label: '最終收斂解 (G_{eq}, ξ_{eq})', c: 'blue', w: 2.3 }
    ]);
    const pCurve = new Plot(grid, { h: 270, xlim: [-4, 0], ylim: [0, 1.1], xlabel: '有效剪應變 log₁₀(γ_{eff} [%])（10⁻⁴% ~ 1%）', ylabel: 'G/G_{max} 與 阻尼比 ξ/0.30' });
    pCurve.legend([
      { label: '模數折減曲線 G/G_{max}(γ)', c: 'blue', w: 2 },
      { label: '阻尼曲線 ξ(γ) / 30%', c: 'red', dash: true, w: 2 },
      { label: '三層收斂工作點 (層 1, 2, 3)', c: 'blue', line: false, marker: 'circle', fill: true }
    ]);

    const logTbl = SD.div(body, '');
    const st = SD.stats(body, [
      { k: 'rg', label: '折減係數 R_γ = (M−1)/10' },
      { k: 'nit', label: '收斂迭代輪數 (<1%)', hl: true },
      { k: 'f0', label: '初始基本頻率 f₁^{(0)} (Hz)' },
      { k: 'fconv', label: '收斂基本頻率 f₁^{(conv)} (Hz)', hl: true },
      { k: 'pgaOut', label: '收斂地表 PGA_{out} (g)', hl: true },
      { k: 'gmax3', label: '底層有效剪應變 γ_{eff} (%)' }
    ]);

    const BASE_LAYERS = [
      { h: 6, vs0: 160, xi0: 0.02, rho: 1.80, gr: 0.04 },
      { h: 8, vs0: 240, xi0: 0.02, rho: 1.92, gr: 0.06 },
      { h: 10, vs0: 340, xi0: 0.015, rho: 2.05, gr: 0.08 }
    ];
    const ROCK = { h: 0, vs: 1000, xi: 0.01, rho: 2.40 };

    function ggCurve(gamPct, gr) { return 1 / (1 + gamPct / gr); }
    function xiCurve(gamPct, gr, xi0) { return xi0 + 0.22 * (1 - 1 / (1 + gamPct / gr)); }

    function draw() {
      const pgaIn = ctl.pga, Rg = (ctl.M - 1) / 10;
      const cur = BASE_LAYERS.map(l => ({ ...l, vs: l.vs0, xi: l.xi0, gg: 1, gamMax: 0, gamEff: 0 }));
      const hist = [];

      for (let it = 1; it <= 8; it++) {
        const full = [...cur, ROCK];
        const ts = cur.reduce((s, l) => s + l.h / l.vs, 0);
        const f1 = 1 / (4 * ts);
        const re = new Float64Array(N), im = new Float64Array(N);
        for (let k = 1; k <= N / 2; k += 2) {
          const fq = k * df;
          if (fq > 25) break;
          const tf = tfMulti(fq, full);
          re[k] = pgaIn * (baseRe[k] * tf.f3r - baseIm[k] * tf.f3i);
          im[k] = pgaIn * (baseRe[k] * tf.f3i + baseIm[k] * tf.f3r);
          re[N - k] = re[k]; im[N - k] = -im[k];
        }
        fft(re, im, 1);
        let pgaOut = 0;
        for (let i = 0; i < N; i++) {
          const a = Math.abs(re[i] / N) * 2;
          if (a > pgaOut) pgaOut = a;
        }

        let zTop = 0, sigV = 0, maxErr = 0;
        const snap = [];
        for (let m = 0; m < cur.length; m++) {
          const l = cur[m];
          const zMid = zTop + 0.5 * l.h;
          const sigMid = sigV + 0.5 * l.h * l.rho * 9.81;
          const rd = 1 - 0.012 * zMid;
          const tauMax = 0.65 * sigMid * (0.55 * pgaOut + 0.45 * pgaIn) * rd;
          const G = l.rho * l.vs * l.vs;
          const gamMaxPct = (tauMax / G) * 100;
          const gamEffPct = Rg * gamMaxPct;
          const ggNew = ggCurve(gamEffPct, l.gr);
          const xiNew = xiCurve(gamEffPct, l.gr, l.xi0);
          const err = Math.max(Math.abs(ggNew - l.gg) / l.gg, Math.abs(xiNew - l.xi) / l.xi) * 100;
          if (err > maxErr) maxErr = err;
          l.gamMax = gamMaxPct; l.gamEff = gamEffPct;
          snap.push({ m: m + 1, gamMax: gamMaxPct, gamEff: gamEffPct, gg: l.gg, xi: l.xi * 100, vs: l.vs });
          l.gg = ggNew; l.vs = l.vs0 * Math.sqrt(ggNew); l.xi = xiNew;
          zTop += l.h; sigV += l.h * l.rho * 9.81;
        }
        hist.push({ it, f1, pgaOut, maxErr, snap, layers: JSON.parse(JSON.stringify(full)) });
        if (it >= 2 && maxErr < 1.0) break;
      }

      const first = hist[0], last = hist[hist.length - 1];
      const fs = linspace(0.05, 8, 240);
      pTF.clear();
      pTF.line(fs, fs.map(fq => tfMulti(fq, first.layers).f3), { c: 'gray', dash: '4 3', w: 1.6 });
      if (hist.length > 1) pTF.line(fs, fs.map(fq => tfMulti(fq, hist[1].layers).f3), { c: 'red', dash: '2 3', w: 1.5 });
      pTF.line(fs, fs.map(fq => tfMulti(fq, last.layers).f3), { c: 'blue', w: 2.3 });

      pCurve.clear();
      pCurve.fn(lg => ggCurve(Math.pow(10, lg), 0.06), { c: 'blue', w: 2, n: 200 });
      pCurve.fn(lg => xiCurve(Math.pow(10, lg), 0.06, 0.02) / 0.30, { c: 'red', dash: '4 3', w: 2, n: 200 });
      const lgs = last.snap.map(s => Math.max(-4, Math.min(0, Math.log10(s.gamEff))));
      pCurve.markers(lgs, last.snap.map(s => s.gg), { c: 'blue', fill: true, r: 5, over: true });
      pCurve.markers(lgs, last.snap.map(s => (s.xi / 100) / 0.30), { c: 'red', fill: true, r: 4.5, over: true });

      logTbl.innerHTML = `<table class="data"><thead><tr>
        <th>迭代輪次</th><th>基本頻率 f₁ (Hz)</th><th class="num">層 1 γ<sub>eff</sub> (%)</th><th class="num">層 2 γ<sub>eff</sub> (%)</th><th class="num">層 3 γ<sub>eff</sub> (%)</th><th class="num">層 3 G/G<sub>max</sub></th><th class="num">層 3 ξ (%)</th><th class="num">最大差值 (%)</th>
      </tr></thead><tbody>${hist.map(h => `<tr>
        <td>Iteration ${h.it}</td>
        <td class="num">${f(h.f1, 2)}</td>
        <td class="num">${f(h.snap[0].gamEff, 4)}</td>
        <td class="num">${f(h.snap[1].gamEff, 4)}</td>
        <td class="num">${f(h.snap[2].gamEff, 4)}</td>
        <td class="num">${f(h.snap[2].gg, 3)}</td>
        <td class="num">${f(h.snap[2].xi, 2)}</td>
        <td class="num">${f(h.maxErr, 2)}</td>
      </tr>`).join('')}</tbody></table>`;

      st({
        rg: f(Rg, 2),
        nit: `${hist.length} 輪`,
        f0: f(first.f1, 2),
        fconv: f(last.f1, 2),
        pgaOut: f(last.pgaOut, 3),
        gmax3: f(last.snap[2].gamEff, 3) + ' %'
      });
    }
    draw();
  })();

  /* =========================================================
     11. #w-check：子層厚度與隱含強度檢核 (6.11)
     ========================================================= */
  (() => {
    const body = SD.widget('w-check'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'h', label: '子層厚度 h (m)', min: 0.5, max: 12.0, step: 0.5, value: 6.0, fmt: v => v.toFixed(1) },
      { k: 'vs', label: '折減後波速 V_s (m/s)', min: 60, max: 300, step: 5, value: 100, fmt: v => v.toFixed(0) },
      { k: 'fmax', label: '目標最高頻率 f_max (Hz)', min: 10, max: 30, step: 1, value: 20, fmt: v => v.toFixed(0) },
      { k: 'sigv', label: '有效垂直應力 σ\'_v (kPa)', min: 30, max: 300, step: 10, value: 100, fmt: v => v.toFixed(0) },
      { k: 'ocr', label: '過壓密比 OCR', min: 1.0, max: 4.0, step: 0.2, value: 1.0, fmt: v => v.toFixed(1) },
      { k: 'gam', label: '等值線性峰值剪應變 γ_max (%)', min: 0.05, max: 1.80, step: 0.05, value: 0.80, fmt: v => v.toFixed(2) + ' %' }
    ], draw);

    SD.buttons(body, '檢核案例', [
      { label: 'NC 黏土例 (h=6m, σ\'v=100, OCR=1, γ=0.8% → 雙重違規)', on: () => ctl.set({ h: 6.0, vs: 100, fmax: 20, sigv: 100, ocr: 1.0, gam: 0.80 }) },
      { label: '切細子層並提高強度 (h=1.2m, OCR=2.4, γ=0.25% → 通過)', on: () => ctl.set({ h: 1.2, vs: 120, fmax: 20, sigv: 100, ocr: 2.4, gam: 0.25 }) }
    ]);

    const p = new Plot(body, { h: 270, xlim: [0, 2.0], ylim: [0, 80], xlabel: '剪應變 γ (%)', ylabel: '剪應力 τ (kPa)' });
    p.legend([
      { label: '真實非線性骨架曲線 τ(γ)', c: 'blue', w: 2.2 },
      { label: '等值線性割線 τ = G_{eq} γ', c: 'red', w: 2 },
      { label: '靜態不排水強度上限 s_u = 0.22 σ\'_v OCR^{0.8}', c: 'ink', dash: true }
    ]);
    const st = SD.stats(body, [
      { k: 'fcut', label: '子層截止頻率 V_s/(4h) (Hz)', hl: true },
      { k: 'hmax', label: '允許最大層厚 V_s/(4f_{max}) (m)' },
      { k: 'su', label: '靜態不排水強度 s_u (kPa)', hl: true },
      { k: 'tauEql', label: '等值線性峰值應力 τ_{EQL} (kPa)', hl: true },
      { k: 'ratio', label: '應力強度比 τ_{EQL} / s_u' },
      { k: 'status', label: '檢核結論', hl: true }
    ]);

    function draw() {
      const h = ctl.h, vs = ctl.vs, fmax = ctl.fmax, sigv = ctl.sigv, ocr = ctl.ocr, gam = ctl.gam;
      const rho = 1.8;
      const Geq = rho * vs * vs;
      const Gmax = Geq * 2.4;
      const su = 0.22 * sigv * Math.pow(ocr, 0.8);
      const tauEql = Geq * (gam / 100);
      const fcut = vs / (4 * h);
      const hAllow = vs / (4 * fmax);

      const ym = Math.max(40, Math.ceil(Math.max(su, tauEql) * 1.25 / 10) * 10);
      p.setLim([0, 2.0], [0, ym]);
      p.clear();
      p.hline(su, { c: 'ink', dash: '5 4', w: 1.6 });
      p.text(1.95, su + ym * 0.03, `靜態強度上限 s_u = ${su.toFixed(1)} kPa`, { anchor: 'end', cls: 'txt' });
      p.fn(gPct => {
        const g = gPct / 100;
        return (Gmax * g) / (1 + (Gmax * g) / su);
      }, { c: 'blue', w: 2.2, n: 250 });
      p.line([0, gam], [0, tauEql], { c: 'red', w: 2.2 });
      p.markers([gam], [tauEql], { c: 'red', fill: true, r: 5, over: true });

      const okFreq = fcut >= fmax;
      const okStr = tauEql <= su;
      const msg = (okFreq && okStr) ? '全數通過' : (!okFreq && !okStr) ? '子層太厚 + 應力超限!' : !okFreq ? `子層應切為 ≥${Math.ceil(h / hAllow)} 層` : 'τ_EQL > s_u（建議改跑非線性）';

      st({
        fcut: `${f(fcut, 1)} (${okFreq ? '≥' : '<'} ${fmax})`,
        hmax: f(hAllow, 2),
        su: f(su, 1),
        tauEql: f(tauEql, 1),
        ratio: f(tauEql / su, 2),
        status: msg
      });
    }
    draw();
  })();

  /* =========================================================
     12. #w-nl：同一個土柱：等值線性 vs 非線性（顯式差分） (6.12)
     ========================================================= */
  (() => {
    const body = SD.widget('w-nl'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'pga', label: '輸入基底脈衝振幅 PGA_in (g)', min: 0.05, max: 0.60, step: 0.05, value: 0.35, fmt: v => v.toFixed(2) + ' g' },
      { k: 'su', label: '土層不排水剪力強度 s_u (kPa)', min: 25, max: 100, step: 5, value: 45, fmt: v => v.toFixed(0) },
      { k: 'xi0', label: '小應變黏性阻尼 ξ₀ (%)', min: 1, max: 8, step: 0.5, value: 2.0, fmt: v => v.toFixed(1) + ' %' }
    ], draw);

    SD.buttons(body, '強度對照', [
      { label: '弱震輸入 (PGA = 0.05 g，兩者幾乎重合)', on: () => ctl.set({ pga: 0.05, su: 50 }) },
      { label: '強震輸入 (PGA = 0.35 g，非線性削頂、等值線性超限)', on: () => ctl.set({ pga: 0.35, su: 45 }) },
      { label: '軟土極強震 (PGA = 0.55 g, s_u = 30 kPa，基底滑移隔震)', on: () => ctl.set({ pga: 0.55, su: 30 }) }
    ]);

    const grid = SD.div(body, 'plot-row plot-grid');
    const pAcc = new Plot(grid, { h: 270, xlim: [0, 3.0], ylim: [-0.8, 0.8], xlabel: '時間 t (s)', ylabel: '地表加速度 a_{surface}(t) (g)' });
    pAcc.legend([
      { label: '時域非線性 (Masing + 顯式差分)', c: 'red', w: 2.2 },
      { label: '等值線性 (固定割線 G_{eq}, ξ_{eq})', c: 'blue', dash: true, w: 1.8 },
      { label: '基底輸入波 a_{in}(t)', c: 'gray', dash: '2 3' }
    ]);
    const pLoop = new Plot(grid, { h: 270, xlim: [-1.2, 1.2], ylim: [-70, 70], xlabel: '中層剪應變 γ (%)', ylabel: '中層剪應力 τ (kPa)' });
    pLoop.legend([
      { label: '時域非線性真實遲滯路徑', c: 'red', w: 2 },
      { label: '等值線性割線與橢圓', c: 'blue', dash: true, w: 1.8 },
      { label: '強度上限 ±s_u', c: 'gray', dash: '2 3' }
    ]);
    const st = SD.stats(body, [
      { k: 'pgaIn', label: '輸入 PGA_{in} (g)' },
      { k: 'pgaEql', label: '等值線性地表 PGA (g)', hl: true },
      { k: 'pgaNl', label: '非線性地表 PGA (g)', hl: true },
      { k: 'tauEql', label: '等值線性峰值剪應力 (kPa)', hl: true },
      { k: 'tauNl', label: '非線性峰值剪應力 (kPa)', hl: true },
      { k: 'gamNl', label: '非線性中層峰值應變 (%)' }
    ]);

    function draw() {
      const H = 16, NL = 8, dz = H / NL, rho = 1.8, vs0 = 160;
      const Gmax = rho * vs0 * vs0;
      const su = ctl.su, gr = su / Gmax;
      const pga = ctl.pga, xi0 = ctl.xi0 / 100;
      const w1 = TAU * vs0 / (4 * H);
      const eta = 2 * xi0 * Gmax / w1;

      const dt = 0.001, nt = 3000;
      const baseA = new Float64Array(nt);
      const baseV = new Float64Array(nt);
      for (let n = 0; n < nt; n++) {
        const t = n * dt;
        const env = t < 2.2 ? Math.sin(Math.PI * t / 2.2) ** 2 : 0;
        baseA[n] = pga * 9.81 * env * Math.sin(TAU * 2.1 * t);
        if (n > 0) baseV[n] = baseV[n - 1] + 0.5 * (baseA[n - 1] + baseA[n]) * dt;
      }

      const v = new Float64Array(NL + 1);
      const gam = new Float64Array(NL);
      const tau = new Float64Array(NL);
      const gRev = new Float64Array(NL);
      const tRev = new Float64Array(NL);
      const dir = new Int8Array(NL);
      const gMaxEnv = new Float64Array(NL);

      const ts = [], aSurfNl = [], aInPlot = [], loopG = [], loopT = [];
      let pgaNl = 0, tauMaxNl = 0, gamMaxMid = 0;
      const mid = Math.floor(NL / 2);

      for (let n = 0; n < nt; n++) {
        v[NL] = baseV[n];
        for (let j = 0; j < NL; j++) {
          const dg = (v[j] - v[j + 1]) / dz * dt;
          const gNew = gam[j] + dg;
          const d = dg >= 0 ? 1 : -1;
          if (dir[j] !== 0 && d !== dir[j]) {
            gRev[j] = gam[j];
            tRev[j] = tau[j];
          }
          dir[j] = d;
          let tHyst;
          if (Math.abs(gNew) >= gMaxEnv[j]) {
            gMaxEnv[j] = Math.abs(gNew);
            tHyst = (Gmax * gNew) / (1 + Math.abs(gNew) / gr);
          } else {
            const dG = gNew - gRev[j];
            tHyst = tRev[j] + (Gmax * dG) / (1 + Math.abs(dG) / (2 * gr));
            const tEnv = (Gmax * gNew) / (1 + Math.abs(gNew) / gr);
            if (d > 0 && tHyst > tEnv) tHyst = tEnv;
            if (d < 0 && tHyst < tEnv) tHyst = tEnv;
          }
          gam[j] = gNew;
          tau[j] = tHyst;
        }
        let a0 = 0;
        for (let j = 0; j < NL; j++) {
          const tBot = tau[j] + eta * (v[j] - v[j + 1]) / dz;
          const tTop = j === 0 ? 0 : tau[j - 1] + eta * (v[j - 1] - v[j]) / dz;
          const acc = (tTop - tBot) / (rho * (j === 0 ? 0.5 * dz : dz));
          if (j === 0) a0 = acc / 9.81;
          v[j] += acc * dt;
        }
        if (Math.abs(a0) > pgaNl) pgaNl = Math.abs(a0);
        if (Math.abs(tau[mid]) > tauMaxNl) tauMaxNl = Math.abs(tau[mid]);
        if (Math.abs(gam[mid]) > gamMaxMid) gamMaxMid = Math.abs(gam[mid]);

        if (n % 5 === 0) {
          ts.push(n * dt);
          aSurfNl.push(a0);
          aInPlot.push(baseA[n] / 9.81);
          loopG.push(gam[mid] * 100);
          loopT.push(tau[mid]);
        }
      }

      const gEff = 0.65 * Math.max(1e-5, gamMaxMid);
      const ggEq = 1 / (1 + gEff / gr);
      const xiEq = xi0 + (2 / Math.PI) * (1 - ggEq) * 0.6;
      const Geq = Gmax * ggEq;
      const etaEq = 2 * xiEq * Geq / (w1 * Math.sqrt(ggEq));

      const vE = new Float64Array(NL + 1), gE = new Float64Array(NL);
      const aSurfEql = [];
      let pgaEql = 0, gMaxEql = 0, tauMaxEql = 0;
      for (let n = 0; n < nt; n++) {
        vE[NL] = baseV[n];
        for (let j = 0; j < NL; j++) gE[j] += (vE[j] - vE[j + 1]) / dz * dt;
        let a0 = 0;
        for (let j = 0; j < NL; j++) {
          const tBot = Geq * gE[j] + etaEq * (vE[j] - vE[j + 1]) / dz;
          const tTop = j === 0 ? 0 : Geq * gE[j - 1] + etaEq * (vE[j - 1] - vE[j]) / dz;
          const acc = (tTop - tBot) / (rho * (j === 0 ? 0.5 * dz : dz));
          if (j === 0) a0 = acc / 9.81;
          vE[j] += acc * dt;
        }
        if (Math.abs(a0) > pgaEql) pgaEql = Math.abs(a0);
        if (Math.abs(gE[mid]) > gMaxEql) gMaxEql = Math.abs(gE[mid]);
        const tCur = Math.abs(Geq * gE[mid]);
        if (tCur > tauMaxEql) tauMaxEql = tCur;
        if (n % 5 === 0) aSurfEql.push(a0);
      }

      const ymA = Math.max(0.2, Math.ceil(Math.max(pgaNl, pgaEql, pga) * 1.15 * 10) / 10);
      pAcc.setLim([0, 3.0], [-ymA, ymA]);
      pAcc.clear();
      pAcc.line(ts, aInPlot, { c: 'gray', dash: '2 3', w: 1.2 });
      pAcc.line(ts, aSurfEql, { c: 'blue', dash: '5 4', w: 1.7 });
      pAcc.line(ts, aSurfNl, { c: 'red', w: 2.1 });

      const xmG = Math.max(0.15, Math.ceil(Math.max(gamMaxMid, gMaxEql) * 100 * 1.2 * 10) / 10);
      const ymT = Math.max(40, Math.ceil(Math.max(su, tauMaxEql) * 1.2 / 10) * 10);
      pLoop.setLim([-xmG, xmG], [-ymT, ymT]);
      pLoop.clear();
      pLoop.hline(su, { c: 'gray', dash: '3 3' }); pLoop.hline(-su, { c: 'gray', dash: '3 3' });
      const gEpct = gMaxEql * 100;
      const ths = linspace(0, TAU, 90);
      pLoop.line(ths.map(t => gEpct * Math.cos(t)), ths.map(t => Geq * gMaxEql * Math.cos(t) + 2 * xiEq * Geq * gMaxEql * Math.sin(t)), { c: 'blue', dash: '5 4', w: 1.8 });
      pLoop.line([-gEpct, gEpct], [-Geq * gMaxEql, Geq * gMaxEql], { c: 'blue', dash: '2 3', w: 1.3 });
      pLoop.line(loopG, loopT, { c: 'red', w: 2 });

      st({
        pgaIn: f(pga, 2),
        pgaEql: f(pgaEql, 2),
        pgaNl: f(pgaNl, 2),
        tauEql: `${f(tauMaxEql, 1)}${tauMaxEql > su ? ' (>s_u!)' : ''}`,
        tauNl: `${f(tauMaxNl, 1)} (≤s_u)`,
        gamNl: f(gamMaxMid * 100, 2) + ' %'
      });
    }
    draw();
  })();

  /* =========================================================
     13. #w-dam：土壩剪力梁振態與 Rayleigh 阻尼 (6.13)
     ========================================================= */
  (() => {
    const body = SD.widget('w-dam'); if (!body) return;
    let tabIdx = 0, p1 = null, p2 = null, ctl = null, st = null;
    SD.tabs(body, ['① 土壩剪力梁振態 (Gazetas 1982)', '② 有限元素 Rayleigh 阻尼 ξ(f)'], i => {
      tabIdx = i;
      if (p1) { syncUI(); draw(); }
    });

    ctl = SD.controls(body, [
      { k: 'modeN', label: '檢視振態階數 n（頁籤 ①）', min: 1, max: 3, step: 1, value: 1, fmt: v => `第 ${v} 振態` },
      { k: 'xiT', label: '目標阻尼比 ξ_target (%)（頁籤 ②）', min: 1, max: 15, step: 0.5, value: 5.0, fmt: v => v.toFixed(1) + ' %' },
      { k: 'fi', label: '匹配頻率 f_i (Hz)（頁籤 ②）', min: 0.3, max: 5.0, step: 0.1, value: 1.0, fmt: v => v.toFixed(1) },
      { k: 'fj', label: '第二匹配頻率 f_j (Hz)（頁籤 ②）', min: 1.0, max: 20.0, step: 0.5, value: 6.0, fmt: v => v.toFixed(1) }
    ], draw);

    const grid = SD.div(body, 'plot-row plot-grid');
    p1 = new Plot(grid, { h: 270, xlim: [-1.15, 1.15], ylim: [-1, 0], yfmt: v => Math.abs(v).toFixed(1), xlabel: '正規化振態位移 φ_n(z)', ylabel: '正規化深度 z / H（0=壩頂/地表，1=基底）' });
    p2 = new Plot(grid, { h: 270, xlim: [0, 15], ylim: [0, 15], xlabel: '頻率 f (Hz)', ylabel: '阻尼比 ξ(f) (%)' });
    st = SD.stats(body, [
      { k: 's1', label: '特徵參數 1', hl: true },
      { k: 's2', label: '特徵參數 2', hl: true },
      { k: 's3', label: '頻率指標 1' },
      { k: 's4', label: '頻率指標 2' },
      { k: 's5', label: '工程重點', hl: true }
    ]);

    function besselJ0(x) {
      const ax = Math.abs(x);
      if (ax < 8) {
        const y = x * x;
        const a1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
        const a2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
        return a1 / a2;
      }
      const z = 8 / ax, xx = ax - 0.785398164;
      return Math.sqrt(0.636619772 / ax) * Math.cos(xx);
    }

    function syncUI() {
      const boxes = body.querySelectorAll('.ctrl');
      if (boxes.length >= 4) {
        boxes[0].style.display = tabIdx === 0 ? '' : 'none';
        boxes[1].style.display = tabIdx === 1 ? '' : 'none';
        boxes[2].style.display = tabIdx === 1 ? '' : 'none';
        boxes[3].style.display = tabIdx === 1 ? '' : 'none';
      }
    }

    function draw() {
      p1.clear(); p2.clear();
      const zs = linspace(0, 1, 150), nzs = zs.map(z => -z);
      const bRoots = [2.4048, 5.5201, 8.6537];

      if (tabIdx === 0) {
        const n = Math.round(ctl.modeN);
        const bn = bRoots[n - 1];
        const kn = (2 * n - 1) * Math.PI / 2;
        p1.legend([{ label: '土壩（楔形剪力梁）J₀(β_n z/H)', c: 'red', w: 2.3 }, { label: '均勻水平土柱 cos((2n−1)πz/2H)', c: 'blue', dash: true }]);
        p1.vline(0, { c: 'gray', dash: null, w: 0.7 });
        p1.line(zs.map(z => besselJ0(bn * z)), nzs, { c: 'red', w: 2.3 });
        p1.line(zs.map(z => Math.cos(kn * z)), nzs, { c: 'blue', dash: '5 4', w: 2 });

        p2.o.xlabel = '振態階數 n'; p2.o.ylabel = '自然頻率比 f_n / f_1'; p2.o.xticks = [1, 2, 3, 4];
        p2.setLim([0.5, 4.5], [0, 8]);
        p2.clear();
        p2.legend([{ label: '土壩 1 : 2.30 : 3.60 : 4.90（較密集）', c: 'red', w: 2.2, marker: 'circle', fill: true }, { label: '均勻土柱 1 : 3 : 5 : 7', c: 'blue', dash: true, marker: 'circle', fill: true }]);
        const ns = [1, 2, 3, 4];
        const rSoil = [1, 3, 5, 7];
        const rDam = [1, 5.5201 / 2.4048, 8.6537 / 2.4048, 11.7915 / 2.4048];
        p2.line(ns, rSoil, { c: 'blue', dash: '4 3', w: 2 });
        p2.line(ns, rDam, { c: 'red', w: 2.2 });
        p2.markers(ns, rSoil, { c: 'blue', fill: true, r: 4.5, over: true });
        p2.markers(ns, rDam, { c: 'red', fill: true, r: 5, over: true });

        const sub = ['', '₁', '₂', '₃'][n] || n;
        st({
          s1: `土壩 β${sub} = ${f(bn, 3)}`,
          s2: `土柱 (2n−1)π/2 = ${f(kn, 3)}`,
          s3: `土壩 f${sub}/f₁ = ${f(bn / bRoots[0], 2)}`,
          s4: `土柱 f${sub}/f₁ = ${2 * n - 1}`,
          s5: '壩頂寬度趨零 → 鞭梢放大顯著'
        });
      } else {
        const xiT = ctl.xiT, fi = ctl.fi, fj = Math.max(fi + 0.2, ctl.fj);
        const wi = TAU * fi, wj = TAU * fj;
        const a0 = xiT * (2 * wi * wj) / (wi + wj);
        const a1 = xiT * 2 / (wi + wj);

        p1.legend([{ label: '土壩第 1 振態', c: 'red', w: 2.2 }, { label: '均勻土柱第 1 振態', c: 'blue', dash: true }]);
        p1.vline(0, { c: 'gray', dash: null, w: 0.7 });
        p1.line(zs.map(z => besselJ0(bRoots[0] * z)), nzs, { c: 'red', w: 2.2 });
        p1.line(zs.map(z => Math.cos(Math.PI * 0.5 * z)), nzs, { c: 'blue', dash: '5 4', w: 2 });

        const fmax = Math.max(12, Math.ceil(fj * 1.8));
        const ym = Math.max(12, Math.ceil(xiT * 2.6));
        p2.o.xlabel = '頻率 f (Hz)'; p2.o.ylabel = '阻尼比 ξ(f) (%)'; p2.o.xticks = null;
        p2.setLim([0, fmax], [0, ym]);
        p2.clear();
        p2.legend([{ label: 'Rayleigh 總阻尼', c: 'ink', w: 2.5 }, { label: '質量比例 a₀/(2ω)', c: 'blue', dash: true }, { label: '勁度比例 a₁ω/2', c: 'red', dash: true }, { label: '目標 ξ', c: 'gray', dash: '3 3' }]);
        p2.hline(xiT, { c: 'gray', dash: '3 3', w: 1.2 });
        p2.fn(fq => fq < 0.15 ? NaN : a0 / (2 * TAU * fq), { c: 'blue', dash: '4 3', w: 1.5, n: 300 });
        p2.fn(fq => 0.5 * a1 * TAU * fq, { c: 'red', dash: '4 3', w: 1.5, n: 200 });
        p2.fn(fq => fq < 0.15 ? NaN : 0.5 * xiT * (fi / fq + fq / fi), { c: 'gray', dash: '2 3', w: 1.4, n: 300 });
        p2.fn(fq => fq < 0.15 ? NaN : a0 / (2 * TAU * fq) + 0.5 * a1 * TAU * fq, { c: 'ink', w: 2.5, n: 350 });
        p2.markers([fi, fj], [xiT, xiT], { c: 'red', fill: true, r: 5, over: true });

        const fmin = Math.sqrt(fi * fj);
        const xiMin = a0 / (2 * TAU * fmin) + 0.5 * a1 * TAU * fmin;
        st({
          s1: `a₀ = ${f(a0 / 100, 3)} s⁻¹`,
          s2: `a₁ = ${f(a1 / 100 * 1000, 2)} ms`,
          s3: `最低點頻率 √(f_i f_j) = ${f(fmin, 2)} Hz`,
          s4: `鞍部最小阻尼 = ${f(xiMin, 2)} %`,
          s5: `單頻在 f_j 處會高估至 ${f(0.5 * xiT * (fi / fj + fj / fi), 1)} %`
        });
      }
    }
    syncUI();
    draw();
  })();

});
