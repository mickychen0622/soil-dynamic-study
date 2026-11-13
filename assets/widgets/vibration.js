/* =========================================================
   Ch4 振動理論：互動元件
   ========================================================= */
document.addEventListener('sd:ready', () => {
  const S = SDOF, TAU = Math.PI * 2, f = SD.f;
  const deg = r => r * 180 / Math.PI;
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

  /* SVG 箭頭 */
  function arrow(g, x1, y1, x2, y2, c, w) {
    const L = Math.hypot(x2 - x1, y2 - y1);
    if (L < 0.5) return;
    const ux = (x2 - x1) / L, uy = (y2 - y1) / L, h = Math.min(9, L * 0.45), hw = h * 0.45;
    svgEl('line', { x1, y1, x2: x2 - ux * h * 0.6, y2: y2 - uy * h * 0.6, class: 's-' + c, 'stroke-width': w || 2, 'stroke-linecap': 'round' }, g);
    svgEl('path', { d: `M${x2} ${y2}L${x2 - ux * h - uy * hw} ${y2 - uy * h + ux * hw}L${x2 - ux * h + uy * hw} ${y2 - uy * h - ux * hw}Z`, class: 'f-' + c }, g);
  }
  function svgText(g, x, y, s, cls, anchor) {
    const t = svgEl('text', { x, y, class: cls || 'txt', 'text-anchor': anchor || 'start' }, g);
    svgRich(t, s);
    return t;
  }
  function clearG(g) { while (g.firstChild) g.removeChild(g.firstChild); }

  /* HW3 Q4：由實測資料求隔振墊參數（與 D1_sdof_isolation.py 相同的演算法） */
  const HW3 = !window.HW3_Q4 ? null : (() => {
    const d = window.HW3_Q4, M = 2000, P0 = 500;
    const pk = S.parabolicPeak(d.f, d.u);
    const log = [];
    let fn = pk.x, xi = 0;
    for (let it = 1; it < 30; it++) {
      const k = M * (TAU * fn) ** 2, rs = 1000 * P0 / k, Dm = pk.y / rs;
      const xn = S.xiFromDmax(Dm);
      log.push({ it, fn, rs, Dm, xi: xn });
      if (Math.abs(xn - xi) < 1e-9) break;
      xi = xn; fn = pk.x / Math.sqrt(1 - 2 * xi * xi);
    }
    const lvl = pk.y / Math.SQRT2;
    const f1 = S.crossing(d.f, d.u, lvl, pk.i, -1), f2 = S.crossing(d.f, d.u, lvl, pk.i, 1);
    const xi2 = (f2 - f1) / (f2 + f1), fn2 = pk.x / Math.sqrt(1 - 2 * xi2 * xi2);
    const fa = (fn + fn2) / 2, xa = (xi + xi2) / 2, ka = M * (TAU * fa) ** 2;
    return { d, M, P0, pk, log, m1: { fn, xi, k: M * (TAU * fn) ** 2 }, m2: { fn: fn2, xi: xi2, f1, f2, lvl, k: M * (TAU * fn2) ** 2 },
      fn: fa, xi: xa, k: ka, c: 2 * xa * Math.sqrt(ka * M), rs: 1000 * P0 / ka };
  })();

  /* ---------------- W1 無阻尼自由振動 ---------------- */
  (() => {
    const body = SD.widget('w-free0'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'm', label: '質量 m (kg)', min: 500, max: 5000, step: 50, value: 2000, fmt: v => v.toFixed(0) },
      { k: 'k', label: '勁度 k (kN/m)', min: 100, max: 2000, step: 10, value: 800, fmt: v => v.toFixed(0) },
      { k: 'u0', label: '初始位移 u(0) (mm)', min: -20, max: 20, step: 0.5, value: 10, fmt: v => v.toFixed(1) },
      { k: 'v0', label: '初速 u̇(0) (m/s)', min: -0.5, max: 0.5, step: 0.01, value: 0.3, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '範例', [
      { label: '講義例：2000 kg', on: () => ctl.set({ m: 2000, k: 800, u0: 10, v0: 0.3 }) },
      { label: '只拉開、不給初速', on: () => ctl.set({ u0: 10, v0: 0 }) },
      { label: '只給初速', on: () => ctl.set({ u0: 0, v0: 0.3 }) }
    ]);
    const p = new Plot(body, { h: 300, xlim: [0, 1], ylim: [-40, 40], xlabel: '時間 t (s)', ylabel: '位移 u (mm)' });
    p.legend([{ label: '位移 u(t)', c: 'blue', w: 2 }, { label: '振幅 ±ρ', c: 'gray', dash: true }, { label: '第一個峰', c: 'red', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [
      { k: 'w', label: '角頻率 ω (rad/s)' }, { k: 'f', label: '頻率 f (Hz)' }, { k: 'T', label: '週期 T (s)' },
      { k: 'rho', label: '振幅 ρ (mm)', hl: true }, { k: 'phi', label: '相位角 φ (°)' }, { k: 'vmax', label: '最大速度 (m/s)' }
    ]);
    SD.div(body, 'w-note', '試試看：只改 u(0) 或 u̇(0)，週期完全不變；只改 m 或 k，振幅會因 ω 改變而變（因為 A = u̇(0)/ω）。');
    function draw() {
      const w = Math.sqrt(ctl.k * 1000 / ctl.m), u0 = ctl.u0 / 1000, v0 = ctl.v0;
      const rho = Math.hypot(v0 / w, u0), ph = Math.atan2(v0 / w, u0);
      p.clear();
      p.hline(0, { c: 'gray', dash: null, w: 0.6 });
      p.hline(rho * 1000, { dash: '5 4' }); p.hline(-rho * 1000, { dash: '5 4' });
      p.fn(t => 1000 * S.free(t, w, 0, u0, v0), { c: 'blue', w: 2, n: 600 });
      let tp = ph / w; if (tp < 0) tp += TAU / w;
      p.markers([tp], [rho * 1000], { c: 'red', fill: true, over: true });
      st({ w: f(w, 2), f: f(w / TAU, 2), T: f(TAU / w, 3), rho: f(rho * 1000, 1), phi: f(deg(ph), 1), vmax: f(w * rho, 3) });
    }
    draw();
  })();

  /* ---------------- W2 有阻尼自由振動 ---------------- */
  (() => {
    const body = SD.widget('w-free1'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'xi', label: '阻尼比 ξ', min: 0, max: 2, step: 0.01, value: 0.1, fmt: v => v.toFixed(2) },
      { k: 'v0', label: '初速 u̇(0)/(ωu(0))', min: -1, max: 1, step: 0.05, value: 0, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, 'ξ =', [0, 0.05, 0.2, 0.5, 1, 2].map(x => ({ label: String(x), on: () => ctl.set({ xi: x }) })));
    const p = new Plot(body, { h: 300, xlim: [0, 5], ylim: [-1.1, 1.3], xlabel: '時間 t/T（T = 2π/ω，無阻尼週期）', ylabel: 'u / u(0)' });
    p.legend([{ label: '目前 ξ', c: 'blue', w: 2 }, { label: '包絡線 ±e^{−ξωt}', c: 'gray', dash: true }, { label: '臨界阻尼 ξ = 1（參考）', c: 'ink', dash: '2 3' }]);
    const st = SD.stats(body, [
      { k: 'reg', label: '情形', hl: true }, { k: 'wd', label: 'ω_D / ω' }, { k: 'td', label: 'T_D / T' },
      { k: 'ratio', label: '相鄰峰比 u_n/u_{n+1}' }, { k: 'ncyc', label: '衰減到 1/e 的圈數' }, { k: 't5', label: '|u| < 5 % 的時間 (t/T)' }
    ]);
    function draw() {
      const x = ctl.xi, w = TAU, v0 = ctl.v0 * w;
      p.clear();
      p.hline(0, { c: 'gray', dash: null, w: 0.6 });
      if (x > 0 && x < 1) {
        const amp = Math.hypot((v0 + x * w) / (w * Math.sqrt(1 - x * x)), 1);
        p.fn(t => amp * Math.exp(-x * w * t), { c: 'gray', dash: '5 4', w: 1 });
        p.fn(t => -amp * Math.exp(-x * w * t), { c: 'gray', dash: '5 4', w: 1 });
      }
      if (Math.abs(x - 1) > 0.005) p.fn(t => S.free(t, w, 1, 1, v0), { c: 'ink', dash: '2 3', w: 1.2 });
      p.fn(t => S.free(t, w, x, 1, v0), { c: 'blue', w: 2, n: 800 });
      let t5 = NaN;
      for (let t = 5; t >= 0; t -= 0.002) if (Math.abs(S.free(t, w, x, 1, v0)) > 0.05) { t5 = t; break; }
      const reg = x === 0 ? '無阻尼' : x < 1 ? '次阻尼' : Math.abs(x - 1) < 1e-9 ? '臨界阻尼' : '超阻尼';
      const wd = x < 1 ? Math.sqrt(1 - x * x) : NaN;
      st({
        reg, wd: x < 1 ? f(wd, 3) : '不振動', td: x < 1 ? f(1 / wd, 3) : '—',
        ratio: x > 0 && x < 1 ? f(Math.exp(TAU * x / wd), 3) : '—',
        ncyc: x > 0 && x < 1 ? f(1 / (TAU * x), 2) : '—',
        t5: x === 0 ? '永不' : (t5 >= 4.99 ? '> 5' : f(t5 + 0.002, 2))
      });
    }
    draw();
  })();

  /* ---------------- W3 旋轉向量 ---------------- */
  (() => {
    const body = SD.widget('w-phasor'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'b', label: '頻率比 β', min: 0.1, max: 3, step: 0.01, value: 0.7, fmt: v => v.toFixed(2) },
      { k: 'xi', label: '阻尼比 ξ', min: 0.01, max: 1, step: 0.01, value: 0.2, fmt: v => v.toFixed(2) },
      { k: 'comp', label: '顯示 G₁、G₂ 分量', type: 'check', value: true }
    ], () => { render(); });
    let playing = true;
    const bt = SD.buttons(body, '預設', [
      { label: 'β = 0.5', on: () => ctl.set({ b: 0.5 }) },
      { label: 'β = 1', on: () => ctl.set({ b: 1 }) },
      { label: 'β = 2', on: () => ctl.set({ b: 2 }) },
      { label: '❚❚ 暫停', primary: true, on: b => { playing = !playing; b.textContent = playing ? '❚❚ 暫停' : '▶ 播放'; if (playing) loop(); } }
    ]);
    const W = 660, H = 300;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot-svg' });
    body.appendChild(svg);
    const g = svgEl('g', null, svg);
    const st = SD.stats(body, [
      { k: 'D', label: '放大因子 D = ρ_d/ρ_s', hl: true }, { k: 'phi', label: '相位落後 φ (°)', hl: true },
      { k: 'g1', label: 'G₁/ρ_s（沿外力）' }, { k: 'g2', label: 'G₂/ρ_s（落後 90°）' }, { k: 'reg', label: '頻率區' }
    ]);
    SD.div(body, 'w-note', '左：灰色向量是外力（長度以 ρs = p0/k 為 1），藍色向量是位移 ρd；兩者都以 ω̄ 旋轉，實際的量是它們在水平軸上的投影。下方方塊就是質量，跟著藍色向量的投影走。右：外力與位移各自正規化後的歷時，看位移晚了多少。');
    let theta = 0, last = null;
    function render() {
      clearG(g);
      const b = ctl.b, x = ctl.xi, D = S.D(b, x), ph = S.phi(b, x);
      const cx = 165, cy = 130, R = 105, s = R / Math.max(1, D);
      // 座標軸與軌跡圓
      svgEl('line', { x1: cx - R - 25, x2: cx + R + 25, y1: cy, y2: cy, class: 's-gray', 'stroke-width': 1 }, g);
      svgEl('line', { x1: cx, x2: cx, y1: cy - R - 15, y2: cy + R + 12, class: 's-gray', 'stroke-width': 1 }, g);
      svgEl('circle', { cx, cy, r: D * s, fill: 'none', class: 's-gray', 'stroke-dasharray': '3 4', 'stroke-width': 1 }, g);
      const aF = theta - Math.PI / 2, aU = aF - ph;
      const P = (len, a) => [cx + len * s * Math.cos(a), cy - len * s * Math.sin(a)];
      const [fx, fy] = P(1, aF), [ux, uy] = P(D, aU);
      // 標籤沿「外力的左側法向」（G2 在右側）偏移，避免重疊
      const nL = [-Math.sin(aF), -Math.cos(aF)];            // 螢幕座標中，外力方向的左法向
      const lab = (x, y, off, s2, cls) => svgText(g, x + nL[0] * off, y + nL[1] * off + 4, s2, cls, 'middle');
      if (ctl.comp) {
        const g1 = D * Math.cos(ph), g2 = D * Math.sin(ph);
        const [ax, ay] = P(g1, aF);
        arrow(g, cx, cy, ax, ay, 'ink', 1.4);
        arrow(g, ax, ay, ax + g2 * s * Math.cos(aF - Math.PI / 2), ay - g2 * s * Math.sin(aF - Math.PI / 2), 'red', 1.4);
        lab(ax, ay, 16, 'G_1', 'txt');
        // G2 垂直於外力，往三角形外側（沿外力方向）推開
        const sg = g1 >= 0 ? 1 : -1;
        svgText(g, (ax + ux) / 2 + sg * Math.cos(aF) * 16, (ay + uy) / 2 - sg * Math.sin(aF) * 16 + 4, 'G_2', 'txt-red', 'middle');
      }
      arrow(g, cx, cy, fx, fy, 'gray', 3);
      arrow(g, cx, cy, ux, uy, 'blue', 2.6);
      lab((cx + fx) / 2, (cy + fy) / 2, 16, 'p_0/k', 'txt-muted');
      const ua = [(ux - cx), (uy - cy)], ul = Math.hypot(...ua) || 1;
      svgText(g, ux + ua[0] / ul * 16, uy + ua[1] / ul * 16 + 4, 'ρ_d', 'txt-blue', 'middle');
      // 投影與質量方塊
      const yb = cy + R + 34;
      svgEl('line', { x1: ux, x2: ux, y1: uy, y2: yb - 10, class: 'guide' }, g);
      svgEl('line', { x1: cx - R - 25, x2: cx + R + 25, y1: yb + 10, y2: yb + 10, class: 's-gray', 'stroke-width': 1 }, g);
      svgEl('rect', { x: ux - 14, y: yb - 10, width: 28, height: 20, rx: 2, class: 's-ink f-surface', 'stroke-width': 1.4 }, g);
      svgEl('path', { d: `M${fx} ${yb - 18}l-5 -8h10z`, class: 'f-gray' }, g);
      // 右側歷時
      const x0 = 360, x1 = 645, ym = 130, A = 80;
      svgEl('rect', { x: x0, y: ym - A - 12, width: x1 - x0, height: 2 * A + 24, class: 'axis-frame' }, g);
      svgEl('line', { x1: x0, x2: x1, y1: ym, y2: ym, class: 's-gray', 'stroke-width': .6 }, g);
      let dF = '', dU = '';
      for (let i = 0; i <= 200; i++) {
        const tt = theta - 2 * TAU + 2 * TAU * i / 200;
        const X = x0 + (x1 - x0) * i / 200;
        dF += (i ? 'L' : 'M') + X.toFixed(1) + ' ' + (ym - A * Math.sin(tt)).toFixed(1);
        dU += (i ? 'L' : 'M') + X.toFixed(1) + ' ' + (ym - A * Math.sin(tt - ph)).toFixed(1);
      }
      svgEl('path', { d: dF, fill: 'none', class: 's-gray', 'stroke-width': 2 }, g);
      svgEl('path', { d: dU, fill: 'none', class: 's-blue', 'stroke-width': 2 }, g);
      svgEl('circle', { cx: x1, cy: ym - A * Math.sin(theta), r: 4, class: 'f-gray' }, g);
      svgEl('circle', { cx: x1, cy: ym - A * Math.sin(theta - ph), r: 4, class: 'f-blue' }, g);
      svgText(g, x0 + 6, ym - A - 16, '外力 p/p_0（灰）　位移 u/ρ_d（藍）', 'txt-muted');
      svgText(g, (x0 + x1) / 2, ym + A + 30, '← 兩個外力週期　　現在 →', 'txt-muted', 'middle');
      const reg = b < 0.8 ? '勁度控制' : b <= 1.2 ? '接近共振（阻尼控制）' : '慣性控制';
      st({ D: f(D, 3), phi: f(deg(ph), 1), g1: f(D * Math.cos(ph), 3), g2: f(-D * Math.sin(ph), 3), reg });
    }
    function loop(ts) {
      if (!playing) { last = null; return; }
      if (ts !== undefined) {
        if (last !== null) theta += (ts - last) / 1000 * TAU / 2.5;
        last = ts;
        render();
      }
      requestAnimationFrame(loop);
    }
    render(); loop();
  })();

  /* ---------------- W4 D–β 與 φ–β ---------------- */
  (() => {
    const body = SD.widget('w-curves'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'xi', label: '阻尼比 ξ', min: 0.01, max: 1, step: 0.005, value: 0.1, fmt: v => v.toFixed(3) },
      { k: 'b', label: '讀值位置 β', min: 0, max: 3, step: 0.01, value: 1, fmt: v => v.toFixed(2) }
    ], draw);
    const row = SD.div(body, 'plot-row');
    const pD = new Plot(row, { w: 420, h: 300, xlim: [0, 3], ylim: [0, 6], xlabel: '頻率比 β', ylabel: '放大因子 D' });
    const pP = new Plot(row, { w: 420, h: 300, xlim: [0, 3], ylim: [0, 180], yticks: [0, 45, 90, 135, 180], xlabel: '頻率比 β', ylabel: '相位落後 φ (°)' });
    pD.legend([{ label: '參考 ξ = 0、0.1、0.2、0.5、1.0', c: 'gray' }, { label: '目前 ξ', c: 'blue', w: 2.2 }, { label: '共振峰', c: 'red', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [
      { k: 'D', label: '此 β 的 D', hl: true }, { k: 'phi', label: '此 β 的 φ (°)', hl: true },
      { k: 'pb', label: '峰的位置 √(1−2ξ²)' }, { k: 'pD', label: '峰高 D_{max}' }, { k: 'half', label: '1/(2ξ)' }, { k: 'cross', label: 'D = 1 的交點' }
    ]);
    const refs = [0, 0.1, 0.2, 0.5, 1.0];
    function draw() {
      const x = ctl.xi, b = ctl.b, pk = S.peak(x);
      const ymax = Math.min(30, Math.max(6, pk ? pk.D * 1.15 : 6));
      pD.setLim(null, [0, ymax]);
      pD.clear(); pP.clear();
      refs.forEach(r => {
        pD.fn(bb => S.D(bb, r), { c: 'gray', w: 1, n: 600 });
        pP.fn(bb => deg(r === 0 ? (bb < 1 ? 0 : 180) : S.phi(bb, r)), { c: 'gray', w: 1, n: 600 });
      });
      pD.hline(1, { dash: '2 3' });
      pP.hline(90, { dash: '2 3' });
      pD.vline(b); pP.vline(b);
      pD.fn(bb => S.D(bb, x), { c: 'blue', w: 2.2, n: 800 });
      pP.fn(bb => deg(S.phi(bb, x)), { c: 'blue', w: 2.2, n: 600 });
      if (pk) pD.markers([pk.b], [pk.D], { c: 'red', fill: true, over: true });
      pD.markers([b], [S.D(b, x)], { c: 'blue', fill: true, over: true, r: 4.5 });
      pP.markers([b], [deg(S.phi(b, x))], { c: 'blue', fill: true, over: true, r: 4.5 });
      const cross = 2 - 4 * x * x > 0 ? Math.sqrt(2 - 4 * x * x) : NaN;
      st({
        D: f(S.D(b, x), 3), phi: f(deg(S.phi(b, x)), 1),
        pb: pk ? f(pk.b, 4) : '無峰（ξ ≥ 0.707）', pD: pk ? f(pk.D, 3) : '—', half: f(1 / (2 * x), 3),
        cross: isFinite(cross) ? 'β = ' + f(cross, 3) : '—'
      });
    }
    pD.readout(bb => `β = ${bb.toFixed(3)}　D = ${S.D(bb, ctl.xi).toFixed(3)}`);
    pP.readout(bb => `β = ${bb.toFixed(3)}　φ = ${deg(S.phi(bb, ctl.xi)).toFixed(1)}°`);
    draw();
  })();

  /* ---------------- W5 力多邊形 ---------------- */
  (() => {
    const body = SD.widget('w-poly'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'b', label: '頻率比 β', min: 0.05, max: 2.5, step: 0.01, value: 0.6, fmt: v => v.toFixed(2) },
      { k: 'xi', label: '阻尼比 ξ', min: 0.02, max: 1, step: 0.01, value: 0.2, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '預設', [0.2, 1, 2].map(v => ({ label: 'β = ' + v, on: () => ctl.set({ b: v }) })));
    const W = 660, H = 230;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot-svg' });
    body.appendChild(svg);
    const g = svgEl('g', null, svg);
    const st = SD.stats(body, [
      { k: 'fs', label: '彈簧力（單位 p_0D）' }, { k: 'fd', label: '阻尼力 2ξβ' }, { k: 'fi', label: '慣性力 β²' },
      { k: 'p0', label: '外力 p_0 = 1/D' }, { k: 'phi', label: '外力領先位移 φ (°)' }, { k: 'reg', label: '誰在扛外力', hl: true }
    ]);
    function draw() {
      clearG(g);
      const b = ctl.b, x = ctl.xi, fd = 2 * x * b, fi = b * b;
      const xs = [0, 1, 1, 1 - fi], ys = [0, 0, fd, fd];
      const xmin = Math.min(...xs, 0), xmax = Math.max(...xs, 1), ymax = Math.max(fd, 0.3);
      const s = Math.min(480 / (xmax - xmin), 150 / ymax, 260);
      const ox = 90 + (-xmin) * s + (480 - (xmax - xmin) * s) / 2, oy = H / 2 + fd * s / 2 + 8;
      const X = v => ox + v * s, Y = v => oy - v * s;
      svgEl('line', { x1: 30, x2: 640, y1: oy, y2: oy, class: 'guide' }, g);
      svgText(g, 636, oy + 16, '位移方向', 'txt-muted', 'end');
      arrow(g, X(0), Y(0), X(1), Y(0), 'ink', 2.4);
      arrow(g, X(1), Y(0), X(1), Y(fd), 'blue', 2.4);
      arrow(g, X(1), Y(fd), X(1 - fi), Y(fd), 'red', 2.4);
      arrow(g, X(0), Y(0), X(1 - fi), Y(fd), 'gray', 3.2);
      svgText(g, X(0.5), Y(0) + 18, '彈簧力 1', 'txt', 'middle');
      svgText(g, X(1) + 8, Y(fd / 2) + 4, '阻尼力 ' + fd.toFixed(2), 'txt-blue');
      svgText(g, X(1 - fi / 2), Y(fd) - 8, '慣性力 ' + fi.toFixed(2), 'txt-red', 'middle');
      const pm = Math.hypot(1 - fi, fd);
      svgText(g, X((1 - fi) / 2) - 10, Y(fd / 2) - 8, 'p_0 = ' + pm.toFixed(2), 'txt-muted', 'end');
      const ph = S.phi(b, x), r = 34;
      svgEl('path', { d: `M${X(0) + r} ${Y(0)}A${r} ${r} 0 0 0 ${X(0) + r * Math.cos(ph)} ${Y(0) - r * Math.sin(ph)}`, fill: 'none', class: 's-ink', 'stroke-width': 1 }, g);
      svgText(g, X(0) + r * Math.cos(ph / 2) + 6, Y(0) - r * Math.sin(ph / 2), 'φ', 'txt');
      const reg = b < 0.7 ? '彈簧（勁度控制）' : b <= 1.3 ? '阻尼器（阻尼控制）' : '質量（慣性控制）';
      st({ fs: '1', fd: f(fd, 3), fi: f(fi, 3), p0: f(pm, 3), phi: f(deg(ph), 1), reg });
    }
    draw();
  })();

  /* ---------------- W6 共振成長 ---------------- */
  (() => {
    const body = SD.widget('w-growth'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'xi', label: '阻尼比 ξ', min: 0, max: 0.3, step: 0.005, value: 0.05, fmt: v => v.toFixed(3) },
      { k: 'b', label: '頻率比 β', min: 0.5, max: 1.5, step: 0.01, value: 1, fmt: v => v.toFixed(2) },
      { k: 'N', label: '顯示圈數', min: 5, max: 40, step: 1, value: 20, fmt: v => v.toFixed(0) }
    ], draw);
    SD.buttons(body, 'ξ =', [0, 0.02, 0.05, 0.15].map(v => ({ label: String(v), on: () => ctl.set({ xi: v, b: 1 }) })));
    const p = new Plot(body, { h: 300, xlim: [0, 20], ylim: [-12, 12], xlabel: '圈數 n = t/T', ylabel: 'R(t) = u/ρ_s' });
    p.legend([{ label: '反應 R(t)', c: 'blue', w: 1.6 }, { label: '包絡線', c: 'ink', dash: true }, { label: '穩態振幅 ±D', c: 'gray', dash: '2 3' }]);
    const st = SD.stats(body, [
      { k: 'D', label: '穩態振幅 D', hl: true }, { k: 'n63', label: '到 63 % 的圈數 1/(2πξ)', hl: true },
      { k: 'p2', label: '2 圈時達到穩態的' }, { k: 'p10', label: '10 圈時達到穩態的' }
    ]);
    const note = SD.div(body, 'w-note');
    function draw() {
      const x = ctl.xi, b = ctl.b, N = ctl.N, w = TAU;
      const n = 4000, ts = linspace(0, N, n), R = ts.map(t => S.fromRest(t, w, x, b));
      const D = S.D(b, x);
      let mx = 0; R.forEach(v => { if (Math.abs(v) > mx) mx = Math.abs(v); });
      const yl = Math.max(1.5, isFinite(D) ? Math.min(D, 400) : 0, mx) * 1.12;
      p.setLim([0, N], [-yl, yl]);
      p.clear();
      p.hline(0, { dash: null, w: 0.6 });
      const res = Math.abs(b - 1) < 1e-9;
      if (isFinite(D) && D < 1e6) { p.hline(D, { dash: '2 3' }); p.hline(-D, { dash: '2 3' }); }
      if (res) {
        const env = t => x > 0 ? (1 - Math.exp(-x * w * t)) / (2 * x) : Math.PI * t;
        p.fn(env, { c: 'ink', dash: '6 4', w: 1.1 });
        p.fn(t => -env(t), { c: 'ink', dash: '6 4', w: 1.1 });
      }
      p.line(ts, R, { c: 'blue', w: 1.4 });
      const frac = nn => x > 0 ? (1 - Math.exp(-TAU * x * nn)) : NaN;
      st({
        D: x === 0 && res ? '∞' : f(D, 2),
        n63: x > 0 ? f(1 / (TAU * x), 2) + ' 圈' : '∞',
        p2: res && x > 0 ? f(100 * frac(2), 0) + ' %' : '—',
        p10: res && x > 0 ? f(100 * frac(10), 0) + ' %' : '—'
      });
      note.innerHTML = res
        ? (x === 0 ? '無阻尼：包絡線是直線 ωt/2，每圈增加 π 倍靜態變位，振幅沒有上限。' : '共振：振幅從零以指數方式趨近 1/(2ξ)。阻尼越小，最終越大，但要累積越久。')
        : 'β ≠ 1：暫態（以 ω<sub>D</sub> 振動）與穩態（以 ω̄ 振動）疊加，兩個頻率接近時出現「拍」，暫態消失後才剩穩態。';
    }
    draw();
  })();

  /* ---------------- W7 能量迴圈（w-energy：單元 4；w-hw3q4c：HW3 Q4(c)） ---------------- */
  [['w-energy', false], ['w-hw3q4c', true]].forEach(([wid, HWM]) => {
    const body = SD.widget(wid); if (!body || (HWM && !HW3)) return;
    const ctl = SD.controls(body, [
      { k: 'b', label: '頻率比 β', min: 0.1, max: 2, step: 0.01, value: 1, fmt: v => v.toFixed(2) },
      { k: 'xi', label: '阻尼比 ξ', min: 0.01, max: 0.5, step: 0.0001, value: HWM ? +HW3.xi.toFixed(4) : 0.1, fmt: v => v.toFixed(4) },
      ...(HWM ? [{ k: 'hw', label: '以隔振墊的數值計算焦耳', type: 'check', value: true }] : [])
    ], draw);
    if (HWM) SD.buttons(body, 'Q4(c)', [
      { label: 'β = 1', on: () => ctl.set({ b: 1, xi: +HW3.xi.toFixed(4), hw: true }) },
      { label: 'β = 0.5', on: () => ctl.set({ b: 0.5, xi: +HW3.xi.toFixed(4), hw: true }) }
    ]);
    const p = new Plot(body, { w: 560, h: 340, xlim: [-1.3, 1.3], ylim: [-1.6, 1.6], xlabel: 'u / ρ_d', ylabel: '力 / (kρ_d)' });
    p.legend([{ label: '彈簧力 ku', c: 'gray', dash: true }, { label: '阻尼力 cu̇（橢圓）', c: 'red', w: 1.6 }, { label: '總力 ku + cu̇（傾斜橢圓）', c: 'blue', w: 2 }]);
    const st = SD.stats(body, [
      { k: 'wd', label: 'W_D' }, { k: 'ws', label: 'W_S' }, { k: 'ratio', label: 'W_D / (4πW_S) = ξβ', hl: true }, { k: 'xi', label: '阻尼比 ξ' }
    ]);
    function draw() {
      const b = ctl.b, x = ctl.xi, e = 2 * x * b;
      const th = linspace(0, TAU, 361);
      const us = th.map(Math.sin);
      const yl = Math.max(1.6, Math.hypot(1, e) * 1.2);
      p.setLim(null, [-yl, yl]);
      p.clear();
      p.hline(0, { dash: null, w: 0.6 }); p.vline(0, { dash: null, w: 0.6 });
      p.line([-1.3, 1.3], [-1.3, 1.3], { c: 'gray', dash: '6 4', w: 1.2 });
      p.line(us, th.map(t => e * Math.cos(t)), { c: 'red', w: 1.6 });
      p.line(us, th.map((t, i) => us[i] + e * Math.cos(t)), { c: 'blue', w: 2 });
      let wd, ws;
      if (ctl.hw) {
        const rd = HW3.rs / 1000 * S.D(b, x), wb = b * TAU * HW3.fn;
        const c = 2 * x * Math.sqrt(HW3.k * HW3.M);
        wd = f(Math.PI * c * wb * rd * rd, 4) + ' J'; ws = f(0.5 * HW3.k * rd * rd, 4) + ' J';
      } else {
        wd = f(Math.PI * e, 4) + ' kρ_d²'; ws = '0.5 kρ_d²';
      }
      st({ wd, ws, ratio: f(x * b, 4), xi: f(x, 4) });
    }
    draw();
  });

  /* ---------------- W8 阻尼量測實驗室 ---------------- */
  (() => {
    const body = SD.widget('w-damp'); if (!body) return;
    const panes = SD.tabs(body, ['a. 自由振動衰減', 'b、c. 掃頻（共振放大、半功率）']);

    // --- a. 自由衰減
    (() => {
      const pane = panes[0];
      let seed = 7;
      const ctl = SD.controls(pane, [
        { k: 'xi', label: '真實阻尼比 ξ', min: 0.005, max: 0.15, step: 0.001, value: 0.03, fmt: v => (v * 100).toFixed(1) + ' %' },
        { k: 'noise', label: '峰值讀數誤差', min: 0, max: 0.1, step: 0.005, value: 0.02, fmt: v => '±' + (v * 100).toFixed(1) + ' %' },
        { k: 'm', label: '兩點法隔幾圈 m', min: 1, max: 10, step: 1, value: 5, fmt: v => v.toFixed(0) }
      ], draw);
      SD.buttons(pane, null, [{ label: '重新量一次（新的讀數誤差）', on: () => { seed++; draw(); } }]);
      const row = SD.div(pane, 'plot-row');
      const p1 = new Plot(row, { w: 420, h: 280, xlim: [0, 12], ylim: [-1.1, 1.1], xlabel: '圈數', ylabel: 'u / u_1' });
      const p2 = new Plot(row, { w: 420, h: 280, xlim: [0, 12], ylim: [-2, 0.2], xlabel: '峰的序號 n', ylabel: 'ln(u_n / u_1)' });
      p1.legend([{ label: '模擬歷時', c: 'blue' }, { label: '量到的峰值', c: 'ink', line: false, marker: 'circle' }]);
      p2.legend([{ label: '量測', c: 'ink', line: false, marker: 'circle' }, { label: '迴歸直線', c: 'blue' }, { label: '兩點法', c: 'red', dash: true }]);
      const st = SD.stats(pane, [
        { k: 'tr', label: '真值 ξ' }, { k: 'two', label: '兩點法 ξ', hl: true }, { k: 'e2', label: '兩點法誤差' },
        { k: 'reg', label: '迴歸 ξ', hl: true }, { k: 'er', label: '迴歸誤差' }
      ]);
      SD.div(pane, 'w-note', '把 m 調成 1，再按幾次「重新量一次」：相鄰兩峰幾乎一樣高，讀數誤差直接吃掉差值，估計值亂跳。隔越多圈、或用所有峰做迴歸，結果越穩定。');
      function draw() {
        const x = ctl.xi, wd = Math.sqrt(1 - x * x), dl = TAU * x / wd, rnd = SD.rng(seed * 7919);
        const K = 12;
        const pk = Array.from({ length: K }, (_, k) => Math.exp(-dl * k) * (1 + ctl.noise * rnd.normal()));
        p1.clear(); p2.clear();
        p1.hline(0, { dash: null, w: 0.6 });
        // 時間以阻尼週期 T_D 為單位：ξω t = δ·(t/T_D)
        p1.fn(t => Math.exp(-dl * t) * Math.cos(TAU * t), { c: 'blue', w: 1.3, n: 1500, range: [0, 12] });
        const tk = pk.map((_, k) => k);
        p1.markers(tk, pk, { c: 'ink', over: true });
        const ln = pk.map(v => Math.log(v / pk[0]));
        const ylo = Math.min(-0.2, Math.min(...ln) * 1.15);
        p2.setLim(null, [ylo, 0.15 * Math.abs(ylo)]);
        p2.clear();
        // 迴歸
        const n = K, mx = tk.reduce((a, b) => a + b) / n, my = ln.reduce((a, b) => a + b) / n;
        let sxy = 0, sxx = 0; tk.forEach((t, i) => { sxy += (t - mx) * (ln[i] - my); sxx += (t - mx) ** 2; });
        const slope = sxy / sxx, icpt = my - slope * mx;
        p2.line([0, 11], [icpt, icpt + slope * 11], { c: 'blue', w: 1.6 });
        const m = ctl.m;
        p2.line([0, m], [ln[0], ln[m]], { c: 'red', dash: '6 4', w: 1.6 });
        p2.markers(tk, ln, { c: 'ink', over: true });
        p2.markers([0, m], [ln[0], ln[m]], { c: 'red', fill: true, over: true, r: 3.2 });
        const toXi = d => d / Math.sqrt(TAU * TAU + d * d);
        const x2 = toXi(Math.log(pk[0] / pk[m]) / m), xr = toXi(-slope);
        st({ tr: f(x * 100, 2) + ' %', two: f(x2 * 100, 2) + ' %', e2: f(100 * (x2 / x - 1), 1) + ' %', reg: f(xr * 100, 2) + ' %', er: f(100 * (xr / x - 1), 1) + ' %' },
          { e2: Math.abs(x2 / x - 1) > 0.1, er: Math.abs(xr / x - 1) > 0.1 });
      }
      draw();
    })();

    // --- b、c. 掃頻
    (() => {
      const pane = panes[1];
      const ctl = SD.controls(pane, [
        { k: 'xi', label: '真實阻尼比 ξ', min: 0.005, max: 0.1, step: 0.001, value: 0.02, fmt: v => (v * 100).toFixed(1) + ' %' },
        { k: 'step', label: '頻率間隔 Δβ', min: 0.001, max: 0.05, log: true, value: 0.01, snap: 2, fmt: v => v.toFixed(4) },
        { k: 'off', label: '掃頻格點的位置偏移', min: 0, max: 1, step: 0.01, value: 0.5, fmt: v => (v * 100).toFixed(0) + ' %' }
      ], draw);
      const p = new Plot(pane, { h: 320, xlim: [0.85, 1.15], ylim: [0, 30], xlabel: '頻率比 β = f / f_n', ylabel: '放大因子 D' });
      p.legend([{ label: '真實曲線', c: 'blue' }, { label: '掃頻量到的點', c: 'ink', dash: true, marker: 'circle' }, { label: '峰值 / √2', c: 'gray', dash: true }, { label: '半功率點 f_1、f_2', c: 'red', line: false, marker: 'circle', fill: true }]);
      const st = SD.stats(pane, [
        { k: 'tr', label: '真值 ξ' }, { k: 'dm', label: '量到的峰值 D_max' },
        { k: 'amp', label: 'b. 共振放大 ξ', hl: true }, { k: 'ea', label: '誤差' },
        { k: 'hp', label: 'c. 半功率 ξ', hl: true }, { k: 'eh', label: '誤差' }
      ]);
      SD.div(pane, 'w-note', '共振放大法假設已知靜態變位 ρs。把 Δβ 調大、拖動格點偏移：錯過峰頂時 D_max 偏小，共振放大法立刻高估阻尼；半功率法也受影響，因為 1/√2 的高度由峰值決定。阻尼越小，峰越窄，需要的間隔越細。');
      function draw() {
        const x = ctl.xi, h = Math.max(6 * x, 0.04);
        const lo = 1 - h, hi = 1 + h, dS = ctl.step;
        const bs = []; for (let b = lo + ctl.off * dS; b <= hi; b += dS) bs.push(b);
        const Ds = bs.map(b => S.D(b, x));
        const pk = S.peak(x);
        p.setLim([lo, hi], [0, pk.D * 1.12]);
        p.clear();
        p.fn(b => S.D(b, x), { c: 'blue', w: 1.6, n: 800 });
        p.line(bs, Ds, { c: 'ink', dash: '5 4', w: 1 });
        p.markers(bs, Ds, { c: 'ink', over: true, r: 3.4 });
        let im = 0; Ds.forEach((v, i) => { if (v > Ds[im]) im = i; });
        const Dm = Ds[im], lvl = Dm / Math.SQRT2;
        p.hline(lvl, { dash: '5 4' });
        const b1 = S.crossing(bs, Ds, lvl, im, -1), b2 = S.crossing(bs, Ds, lvl, im, 1);
        if (isFinite(b1) && isFinite(b2)) p.markers([b1, b2], [lvl, lvl], { c: 'red', fill: true, over: true });
        const xa = S.xiFromDmax(Dm), xh = (b2 - b1) / (b2 + b1);
        st({
          tr: f(x * 100, 2) + ' %', dm: f(Dm, 2) + '（真：' + f(pk.D, 2) + '）',
          amp: f(xa * 100, 2) + ' %', ea: f(100 * (xa / x - 1), 1) + ' %',
          hp: isFinite(xh) ? f(xh * 100, 2) + ' %' : '點太少', eh: isFinite(xh) ? f(100 * (xh / x - 1), 1) + ' %' : '—'
        }, { ea: Math.abs(xa / x - 1) > 0.05, eh: !(Math.abs(xh / x - 1) <= 0.05) });
      }
      draw();
    })();

  })();

  /* ---------------- HW3 Q4(a) 實測資料（作業頁） ---------------- */
  (() => {
    const host = SD.widget('w-hw3q4a'); if (!host || !HW3) return;
    (() => {
      const pane = host, d = HW3.d;
      SD.div(pane, 'w-note', '儀器 m = 2,000 kg 放在隔振墊上，振動機施加 p0 = 0.5 kN，量 4–12 Hz 的穩態位移振幅（81 點，Homework 3/!Question/data/sdof_frequency_response.csv）。以下數值都由網頁即時計算，與 HW3 報告相同。');
      const p = new Plot(pane, { h: 330, xlim: [4, 12], ylim: [0, 0.7], xlabel: '頻率 f (Hz)', ylabel: '位移振幅 (mm)' });
      p.legend([{ label: '實測', c: 'ink', dash: true, marker: 'circle' }, { label: 'SDOF（採用值）', c: 'blue', w: 2 }, { label: '峰（三點拋物線）', c: 'red', line: false, marker: 'circle', fill: true }, { label: '峰值 / √2', c: 'gray', dash: true }]);
      p.line(d.f, d.u, { c: 'ink', dash: '5 4', w: 1 });
      p.markers(d.f, d.u, { c: 'ink', r: 3.2 });
      p.fn(ff => HW3.rs * S.D(ff / HW3.fn, HW3.xi), { c: 'blue', w: 2 });
      p.hline(HW3.m2.lvl, { dash: '5 4' });
      p.vline(HW3.m2.f1, { c: 'red', dash: '2 3' }); p.vline(HW3.m2.f2, { c: 'red', dash: '2 3' });
      p.markers([HW3.pk.x], [HW3.pk.y], { c: 'red', fill: true, over: true });
      p.readout(ff => {
        let j = 0; d.f.forEach((v, i) => { if (Math.abs(v - ff) < Math.abs(d.f[j] - ff)) j = i; });
        return `f = ${d.f[j].toFixed(1)} Hz　實測 ${d.u[j].toFixed(4)} mm　模型 ${(HW3.rs * S.D(d.f[j] / HW3.fn, HW3.xi)).toFixed(4)} mm`;
      });
      const t = SD.div(pane, 'tablewrap');
      t.innerHTML = `<table class="data">
        <tr><th>步驟</th><th>(i) 共振放大法</th><th>(ii) 半功率頻寬法</th><th>採用（兩法平均）</th></tr>
        <tr><td>峰</td><td colspan="2">最大點與左右兩點做拋物線內插：f<sub>peak</sub> = ${HW3.pk.x.toFixed(3)} Hz、U<sub>max</sub> = ${HW3.pk.y.toFixed(4)} mm</td><td></td></tr>
        <tr><td>做法</td><td>k = m(2πf<sub>n</sub>)²、ρ<sub>s</sub> = p<sub>0</sub>/k、D<sub>max</sub> = U<sub>max</sub>/ρ<sub>s</sub> = 1/(2ξ√(1−ξ²))，再以 f<sub>peak</sub> = f<sub>n</sub>√(1−2ξ²) 修正 f<sub>n</sub>，疊代 ${HW3.log.length} 次收斂</td>
          <td>U = U<sub>max</sub>/√2 處線性內插：f<sub>1</sub> = ${HW3.m2.f1.toFixed(3)}、f<sub>2</sub> = ${HW3.m2.f2.toFixed(3)} Hz；ξ = (f<sub>2</sub>−f<sub>1</sub>)/(f<sub>2</sub>+f<sub>1</sub>)</td><td></td></tr>
        <tr class="hl"><td>f<sub>n</sub> (Hz)</td><td class="num">${HW3.m1.fn.toFixed(3)}</td><td class="num">${HW3.m2.fn.toFixed(3)}</td><td class="num">${HW3.fn.toFixed(3)}</td></tr>
        <tr class="hl"><td>ξ (%)</td><td class="num">${(HW3.m1.xi * 100).toFixed(2)}</td><td class="num">${(HW3.m2.xi * 100).toFixed(2)}</td><td class="num">${(HW3.xi * 100).toFixed(2)}</td></tr>
        <tr><td>k (kN/m)</td><td class="num">${(HW3.m1.k / 1000).toFixed(0)}</td><td class="num">${(HW3.m2.k / 1000).toFixed(0)}</td><td class="num">${(HW3.k / 1000).toFixed(0)}</td></tr>
        <tr><td>ρ<sub>s</sub> = p<sub>0</sub>/k (mm)</td><td colspan="2"></td><td class="num">${HW3.rs.toFixed(4)}</td></tr>
        <tr><td>c = 2ξ√(km) (kN·s/m)</td><td class="num">${(2 * HW3.m1.xi * Math.sqrt(HW3.m1.k * HW3.M) / 1000).toFixed(2)}</td><td class="num">${(2 * HW3.m2.xi * Math.sqrt(HW3.m2.k * HW3.M) / 1000).toFixed(2)}</td><td class="num">${(HW3.c / 1000).toFixed(2)}</td></tr>
      </table>`;
      const lg = SD.div(pane, '');
      lg.innerHTML = `<details class="derive"><summary>共振放大法的疊代過程</summary><div class="body"><div class="tablewrap"><table class="data">
        <tr><th>次</th><th class="num">f<sub>n</sub> (Hz)</th><th class="num">ρ<sub>s</sub> (mm)</th><th class="num">D<sub>max</sub></th><th class="num">ξ (%)</th></tr>
        ${HW3.log.map(r => `<tr><td>${r.it}</td><td class="num">${r.fn.toFixed(4)}</td><td class="num">${r.rs.toFixed(5)}</td><td class="num">${r.Dm.toFixed(4)}</td><td class="num">${(r.xi * 100).toFixed(4)}</td></tr>`).join('')}
        </table></div><p>第一次以 f<sub>peak</sub> 當 f<sub>n</sub>；之後用算出的 ξ 把峰的位置修正回自然頻率。兩法差約 0.2 %，取平均作為 (b)、(c) 的採用值。</p></div></details>`;
    })();
  })();

  /* ---------------- W9 傳遞率 ---------------- */
  (() => {
    const body = SD.widget('w-tr'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'xi', label: '阻尼比 ξ', min: 0.01, max: 0.7, step: 0.01, value: 0.05, fmt: v => v.toFixed(2) },
      { k: 'b', label: '運轉時的 β', min: 1.2, max: 6, step: 0.05, value: 3.5, fmt: v => v.toFixed(2) },
      { k: 'rel', label: '同時畫相對位移 β²D（地動情況）', type: 'check', value: false },
      { k: 'rpm', label: '機器轉速 (rpm)', type: 'number', value: 1500 },
      { k: 'm', label: '機器＋機架質量 (kg)', type: 'number', value: 500 },
      { k: 'p0', label: '運轉時的不平衡力 (N)', type: 'number', value: 300 }
    ], draw);
    SD.buttons(body, '範例', [
      { label: '風扇隔振（講義例）', on: () => ctl.set({ rpm: 1500, m: 500, p0: 300, xi: 0.05, b: 3.5 }) },
      { label: '同一台改 ξ = 0.2', on: () => ctl.set({ rpm: 1500, m: 500, p0: 300, xi: 0.2, b: 3.5 }) },
      { label: '1800 rpm、TR ≈ 10 %', on: () => ctl.set({ rpm: 1800, xi: 0.01, b: 3.32 }) }
    ]);
    const p = new Plot(body, { h: 320, xlim: [0, 5], ylim: [0.01, 30], logy: true, xlabel: '頻率比 β', ylabel: '傳遞率 TR' });
    p.legend([{ label: '參考 ξ = 0、0.1、0.2、0.5', c: 'gray' }, { label: '目前 ξ', c: 'blue', w: 2.2 }, { label: '相對位移 β²D', c: 'red', dash: true }, { label: '運轉點', c: 'blue', line: false, marker: 'circle', fill: true }]);
    p.readout(bb => `β = ${bb.toFixed(2)}　TR = ${S.TR(bb, ctl.xi).toFixed(4)}　β²D = ${(bb * bb * S.D(bb, ctl.xi)).toFixed(3)}`);
    const st = SD.stats(body, [
      { k: 'fop', label: '運轉頻率 (Hz)' }, { k: 'fn', label: '支承自然頻率 f_n (Hz)' }, { k: 'k', label: '支承總勁度 k (kN/m)' },
      { k: 'dst', label: '靜態沉陷 (mm)' }, { k: 'tr', label: '運轉時 TR', hl: true }, { k: 'F', label: '運轉時傳到地面的力', hl: true },
      { k: 'trr', label: '啟動通過共振時 TR' }, { k: 'Fr', label: '通過共振時傳遞的力' }
    ]);
    SD.div(body, 'w-note', '不平衡力與轉速平方成正比，所以通過共振（轉速較低）時的不平衡力只有運轉時的 1/β²。比較 ξ = 0.05 與 0.2：運轉時 0.2 傳得較多，但通過共振時少很多，這就是阻尼的折衷。');
    function draw() {
      const x = ctl.xi, b = ctl.b;
      p.clear();
      p.vline(Math.SQRT2, { dash: '4 3' });
      p.hline(1, { dash: '2 3' });
      p.text(Math.SQRT2 + 0.04, 0.016, '√2', { cls: 'txt-muted' });
      [0, 0.1, 0.2, 0.5].forEach(r => p.fn(bb => S.TR(bb, r), { c: 'gray', w: 1, n: 900 }));
      if (ctl.rel) p.fn(bb => bb * bb * S.D(bb, x), { c: 'red', dash: '6 4', w: 1.6, n: 900 });
      p.fn(bb => S.TR(bb, x), { c: 'blue', w: 2.2, n: 900 });
      p.markers([b], [S.TR(b, x)], { c: 'blue', fill: true, over: true, r: 4.5 });
      const fop = ctl.rpm / 60, fn = fop / b, w = TAU * fn, k = ctl.m * w * w;
      const tr = S.TR(b, x), trr = Math.sqrt(1 + 4 * x * x) / (2 * x);
      st({
        fop: f(fop, 2), fn: f(fn, 2), k: f(k / 1000, 0), dst: f(9.81 / (w * w) * 1000, 2),
        tr: f(tr, 4), F: f(tr * ctl.p0, 1) + ' N', trr: f(trr, 2), Fr: f(trr * ctl.p0 / (b * b), 0) + ' N'
      });
    }
    draw();
  })();

  /* ---------------- W10 HW3 Q4(b) ---------------- */
  (() => {
    const body = SD.widget('w-hw3b'); if (!body || !HW3) return;
    const ug0 = 0.05;
    SD.div(body, 'w-note', `隔振墊採用值：f<sub>n</sub> = ${HW3.fn.toFixed(3)} Hz、ξ = ${(HW3.xi * 100).toFixed(2)} %（由 Q4(a) 的實測資料算出）。附近施工使地板振動 u<sub>g0</sub> = 0.05 mm。`);
    const ctl = SD.controls(body, [
      { k: 'fq', label: '地板振動頻率 (Hz)', min: 2, max: 20, step: 0.05, value: 15, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '題目要求', [
      { label: '15 Hz', on: () => ctl.set({ fq: 15 }) },
      { label: '8 Hz', on: () => ctl.set({ fq: 8 }) },
      { label: '√2 f_n', on: () => ctl.set({ fq: +(Math.SQRT2 * HW3.fn).toFixed(2) }) }
    ]);
    const p = new Plot(body, { h: 290, xlim: [2, 20], ylim: [0.1, 10], logy: true, xlabel: '地板振動頻率 f (Hz)', ylabel: 'TR_2 = U_{total}/u_{g0}' });
    p.legend([{ label: 'TR₂（採用值）', c: 'blue', w: 2 }, { label: '目前頻率', c: 'blue', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [
      { k: 'b', label: 'β = f / f_n' }, { k: 'tr', label: 'TR₂', hl: true }, { k: 'u', label: '儀器振幅 U_total (mm)', hl: true }, { k: 'j', label: '判斷' }
    ]);
    function draw() {
      const fq = ctl.fq, b = fq / HW3.fn, tr = S.TR(b, HW3.xi);
      p.clear();
      p.hline(1, { dash: '2 3' });
      p.vline(Math.SQRT2 * HW3.fn, { dash: '4 3' });
      p.text(Math.SQRT2 * HW3.fn + 0.2, 0.14, '√2 f_n = ' + (Math.SQRT2 * HW3.fn).toFixed(2) + ' Hz', { cls: 'txt-muted' });
      p.fn(ff => S.TR(ff / HW3.fn, HW3.xi), { c: 'blue', w: 2, n: 600 });
      p.markers([fq], [tr], { c: 'blue', fill: true, over: true, r: 4.5 });
      st({ b: f(b, 3), tr: f(tr, 3), u: f(ug0 * tr, 4), j: tr < 1 ? '減振' : '放大' }, { j: tr >= 1, tr: tr >= 1 });
    }
    draw();
  })();

  /* ---------------- W11 反應譜 ---------------- */
  (() => {
    const body = SD.widget('w-spec'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'site', label: '地盤', type: 'select', value: 'stiff', options: [['stiff', '硬地盤（卓越頻率約 4 Hz）'], ['soft', '軟弱土層（卓越頻率約 1 Hz）']] },
      { k: 'xi', label: '阻尼比', type: 'select', value: '0.05', options: [['0.02', '2 %'], ['0.05', '5 %'], ['0.1', '10 %'], ['0.2', '20 %']] },
      { k: 'T', label: '讀值週期 T (s)', min: 0.05, max: 4, step: 0.01, value: 0.5, fmt: v => v.toFixed(2) }
    ], draw);
    SD.div(body, 'w-note', '模擬地震：同一段白噪音分別通過兩種地盤濾波（Kanai–Tajimi），加上漸強、持續、衰減的包絡，PGA 都調成 0.3 g。這是教學用的合成紀錄，不是真實測站資料。');
    const dt = 0.005, N = 5000, g0 = 9.81;
    const motions = {};
    function makeMotion(fg, zg) {
      const rnd = SD.rng(2026);
      const wg = TAU * fg;
      let x = 0, v = 0;
      const out = new Array(N);
      for (let i = 0; i < N; i++) {
        const t = i * dt;
        const env = t < 2 ? (t / 2) ** 2 : t < 9 ? 1 : Math.exp(-0.35 * (t - 9));
        const w = rnd.normal() * env;
        const a = -w - 2 * zg * wg * v - wg * wg * x;
        v += a * dt; x += v * dt;
        out[i] = 2 * zg * wg * v + wg * wg * x;
      }
      const pga = Math.max(...out.map(Math.abs));
      return out.map(a => a / pga * 0.3 * g0);
    }
    motions.stiff = makeMotion(4, 0.6);
    motions.soft = makeMotion(1.0, 0.3);
    const Ts = Array.from({ length: 70 }, (_, i) => 0.04 * Math.pow(4 / 0.04, i / 69));
    const cache = {};
    function spectrum(site, xi) {
      const key = site + xi;
      if (!cache[key]) cache[key] = Ts.map(T => { const sd = S.peakDisp(motions[site], dt, T, xi); return { sd, sa: (TAU / T) ** 2 * sd / g0 }; });
      return cache[key];
    }
    const pa = new Plot(body, { h: 200, xlim: [0, 25], ylim: [-0.35, 0.35], xlabel: '時間 t (s)', ylabel: '地表加速度 (g)', m: { b: 44 } });
    const pS = new Plot(body, { h: 300, xlim: [0, 4], ylim: [0, 1.2], xlabel: '自然週期 T (s)', ylabel: '擬加速度 S_a (g)' });
    pS.legend([{ label: '目前地盤', c: 'blue', w: 2 }, { label: '另一種地盤', c: 'gray' }, { label: 'PGA', c: 'ink', dash: '2 3' }]);
    pS.readout(T => {
      const sp = spectrum(ctl.site, +ctl.xi);
      let j = 0; Ts.forEach((v, i) => { if (Math.abs(v - T) < Math.abs(Ts[j] - T)) j = i; });
      return `T = ${Ts[j].toFixed(2)} s　S_a = ${sp[j].sa.toFixed(3)} g　S_d = ${(sp[j].sd * 100).toFixed(2)} cm`;
    });
    const st = SD.stats(body, [
      { k: 'T', label: '週期 T (s)' }, { k: 'sa', label: '擬加速度 S_a (g)', hl: true }, { k: 'amp', label: 'S_a / PGA', hl: true },
      { k: 'sd', label: '最大相對位移 S_d (cm)' }, { k: 'pk', label: '譜的峰值週期 (s)' }
    ]);
    function draw() {
      const site = ctl.site, other = site === 'stiff' ? 'soft' : 'stiff', xi = +ctl.xi;
      const ag = motions[site];
      pa.clear();
      const ts = [], as = [];
      for (let i = 0; i < N; i += 2) { ts.push(i * dt); as.push(ag[i] / g0); }
      pa.line(ts, as, { c: 'ink', w: 0.8 });
      const sp = spectrum(site, xi), so = spectrum(other, xi);
      const ymax = Math.max(1.2, ...sp.map(s => s.sa), ...so.map(s => s.sa)) * 1.08;
      pS.setLim(null, [0, ymax]);
      pS.clear();
      pS.hline(0.3, { c: 'ink', dash: '2 3' });
      pS.line(Ts, so.map(s => s.sa), { c: 'gray', w: 1.2 });
      pS.line(Ts, sp.map(s => s.sa), { c: 'blue', w: 2 });
      const T = ctl.T, sd = S.peakDisp(ag, dt, T, xi), sa = (TAU / T) ** 2 * sd / g0;
      pS.vline(T);
      pS.markers([T], [sa], { c: 'blue', fill: true, over: true, r: 4.5 });
      let jp = 0; sp.forEach((s, i) => { if (s.sa > sp[jp].sa) jp = i; });
      st({ T: f(T, 2), sa: f(sa, 3), amp: f(sa / 0.3, 2), sd: f(sd * 100, 2), pk: f(Ts[jp], 2) });
    }
    draw();
  })();
});
