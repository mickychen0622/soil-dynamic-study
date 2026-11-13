/* =========================================================
   Ch5 動態土壤性質：互動元件
   應變單位：圖上一律用 %；只有算能量（E_D、E_So）時換成小數
   ========================================================= */
document.addEventListener('sd:ready', () => {
  const f = SD.f, PI = Math.PI;
  const logspace = (a, b, n) => Array.from({ length: n }, (_, i) => a * Math.pow(b / a, i / (n - 1)));
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  const fx = (v, d) => isFinite(v) ? v.toFixed(d) : '—';
  const pctFmt = v => v >= 0.01 ? (+v.toPrecision(3)).toString() : v.toExponential(1).replace('e', '×10^');
  const sfmt = v => { // 應變 (%) 的顯示
    if (v >= 0.01) return (+v.toPrecision(3)) + ' %';
    const e = Math.floor(Math.log10(v)), m = v / Math.pow(10, e);
    return (Math.abs(m - 1) < 0.005 ? '' : m.toFixed(2) + ' × ') + '10<sup>' + e + '</sup> %';
  };
  const decFmt = v => { const e = Math.floor(Math.log10(v)), m = v / Math.pow(10, e); return (Math.abs(m - 1) < 0.005 ? '' : m.toFixed(2) + ' × ') + '10<sup>' + e + '</sup>'; };

  /* ---------- 小工具 ---------- */
  function fillPoly(p, xs, ys, cls, op) {
    return svgEl('path', { d: p.path(xs, ys) + 'Z', class: cls || 'f-band', opacity: op, stroke: 'none' }, p.bandG);
  }
  function rectBand(p, x0, x1, y0, y1, cls, op) {
    const X0 = p.sx(x0), X1 = p.sx(x1), Y0 = p.sy(y0), Y1 = p.sy(y1);
    return svgEl('rect', { x: Math.min(X0, X1), y: Math.min(Y0, Y1), width: Math.abs(X1 - X0), height: Math.abs(Y1 - Y0), class: cls, opacity: op }, p.bandG);
  }
  function bisect(fn, a, b) {
    let fa = fn(a);
    for (let i = 0; i < 200; i++) {
      const m = (a + b) / 2, fm = fn(m);
      if ((fm > 0) === (fa > 0)) { a = m; fa = fm; } else b = m;
    }
    return (a + b) / 2;
  }
  function golden(fn, a, b) {
    const g = (Math.sqrt(5) - 1) / 2;
    let c = b - g * (b - a), d = a + g * (b - a);
    for (let i = 0; i < 120; i++) {
      if (fn(c) < fn(d)) b = d; else a = c;
      c = b - g * (b - a); d = a + g * (b - a);
    }
    return (a + b) / 2;
  }
  function nelderMead(fn, x0, step) {
    let s = [x0.slice(), [x0[0] + step[0], x0[1]], [x0[0], x0[1] + step[1]]].map(p => ({ p, v: fn(p) }));
    for (let it = 0; it < 2000; it++) {
      s.sort((a, b) => a.v - b.v);
      if (Math.abs(s[2].v - s[0].v) <= 1e-14 * (Math.abs(s[0].v) + 1e-30)) break;
      const c = [(s[0].p[0] + s[1].p[0]) / 2, (s[0].p[1] + s[1].p[1]) / 2];
      const at = t => [c[0] + t * (s[2].p[0] - c[0]), c[1] + t * (s[2].p[1] - c[1])];
      const r = at(-1), vr = fn(r);
      if (vr < s[0].v) { const e = at(-2), ve = fn(e); s[2] = ve < vr ? { p: e, v: ve } : { p: r, v: vr }; }
      else if (vr < s[1].v) s[2] = { p: r, v: vr };
      else {
        const k = vr < s[2].v ? at(-0.5) : at(0.5), vk = fn(k);
        if (vk < Math.min(vr, s[2].v)) s[2] = { p: k, v: vk };
        else s = s.map((q, i) => i === 0 ? q : (pp => ({ p: pp, v: fn(pp) }))([(q.p[0] + s[0].p[0]) / 2, (q.p[1] + s[0].p[1]) / 2]));
      }
    }
    s.sort((a, b) => a.v - b.v);
    return s[0].p;
  }
  function linreg(xs, ys) {
    const n = xs.length, mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
    let sxy = 0, sxx = 0; xs.forEach((x, i) => { sxy += (x - mx) * (ys[i] - my); sxx += (x - mx) ** 2; });
    const slope = sxy / sxx;
    return { slope, icpt: my - slope * mx };
  }
  function interpLog(x, xs, ys) {
    const lx = Math.log10(x);
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    for (let i = 1; i < xs.length; i++) if (x <= xs[i]) {
      const a = Math.log10(xs[i - 1]), b = Math.log10(xs[i]);
      return ys[i - 1] + (ys[i] - ys[i - 1]) * (lx - a) / (b - a);
    }
    return NaN;
  }

  /* ---------- Darendeli (2001) 平均模型（與 darendeli_curves.py 相同的式子） ---------- */
  const PHI = { p1: 0.0352, p2: 0.0010, p3: 0.3246, p4: 0.3483, a: 0.919, p6: 0.8005, p7: 0.0129, p8: -0.1069, p9: -0.2889, p10: 0.2919, p11: 0.6329, p12: -0.0057 };
  /** 雙曲線（a = 1）Masing 阻尼，x = γ/γr，回傳小數 */
  function masingA1(x) {
    if (x <= 0) return 0;
    const s = x < 0.01 ? x * x / 2 - x ** 3 / 3 + x ** 4 / 4 - x ** 5 / 5 + x ** 6 / 6 : x - Math.log1p(x);
    return (4 * (1 + x) * s / (x * x) - 2) / PI;
  }
  /** 任意曲率 a 的修正雙曲線骨架之 Masing 阻尼（數值積分），回傳小數 */
  function masingNum(x, a) {
    if (x <= 0) return 0;
    const n = 600, L1 = Math.log(x), L0 = L1 - 32, h = (L1 - L0) / n;
    let s = 0;
    for (let i = 0; i <= n; i++) {
      const u = Math.exp(L0 + i * h), v = u * u / (1 + Math.pow(u, a));
      s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * v;
    }
    const A = s * h / 3, fx_ = x / (1 + Math.pow(x, a));
    return 2 / PI * (2 * A / (x * fx_) - 1);
  }
  function DR(g, o) {
    o = Object.assign({ PI: 0, OCR: 1, sm: 1, N: 10, f: 1 }, o || {});
    const gr = (PHI.p1 + PHI.p2 * o.PI * Math.pow(o.OCR, PHI.p3)) * Math.pow(o.sm, PHI.p4);
    const a = PHI.a, gg = 1 / (1 + Math.pow(g / gr, a));
    const d1 = 100 * masingA1(g / gr);
    const c1 = -1.1143 * a * a + 1.8618 * a + 0.2523, c2 = 0.0805 * a * a - 0.0710 * a - 0.0095, c3 = -0.0005 * a * a + 0.0002 * a + 0.0003;
    const dm = c1 * d1 + c2 * d1 * d1 + c3 * d1 ** 3;
    const dmin = (PHI.p6 + PHI.p7 * o.PI * Math.pow(o.OCR, PHI.p8)) * Math.pow(o.sm, PHI.p9) * (1 + PHI.p10 * Math.log(o.f));
    const b = PHI.p11 + PHI.p12 * Math.log(o.N);
    return { gr, gg, d1, dm, dmin, b, D: b * Math.pow(gg, 0.1) * dm + dmin };
  }

  /* ---------- 雙曲線骨架與 Masing 路徑 ---------- */
  /** τ (kPa)，γ (%)，Gmax (kPa)，γr (%) */
  const hyp = (Gmax, gr) => g => Gmax * (g / 100) / (1 + Math.abs(g / gr));
  /** 擴充 Masing rules（反轉點堆疊）；targets 為應變反轉點序列，從 0 開始 */
  function masingExt(targets, fb, nPer) {
    const G = [targets[0]], T = [fb(targets[0])];
    const stack = []; let dir = 0, gc = targets[0], tc = T[0];
    for (let k = 1; k < targets.length; k++) {
      const a = targets[k - 1], b = targets[k];
      if (b === a) continue;
      const d = Math.sign(b - a);
      if (dir !== 0 && d !== dir) stack.push({ g: gc, t: tc });
      dir = d;
      const n = Math.max(20, Math.round(nPer * Math.abs(b - a) / (Math.abs(targets[1] - targets[0]) || 1)));
      for (let i = 1; i <= n; i++) {
        const g = a + (b - a) * i / n;
        for (;;) {
          const L = stack.length;
          if (L >= 2) { const tg = stack[L - 2].g; if (d > 0 ? g >= tg : g <= tg) { stack.length = L - 2; continue; } }
          else if (L === 1 && d * g >= Math.abs(stack[0].g)) { stack.length = 0; continue; }
          break;
        }
        const top = stack[stack.length - 1];
        const t = top ? top.t + 2 * fb((g - top.g) / 2) : fb(g);
        G.push(g); T.push(t); gc = g; tc = t;
      }
    }
    return { g: G, t: T, open: stack.length };
  }
  /** HW3 網頁工具的簡化做法：只記最後一個反轉點；|γ| 超過歷史最大值時回到骨架 */
  function masingSimple(targets, fb, nPer, sameN) {
    const G = [], T = [], S = [];
    let rev = { g: targets[0], t: fb(targets[0]) }, gmp = 0;
    for (let k = 1; k < targets.length; k++) {
      const a = targets[k - 1], b = targets[k];
      const n = sameN ? nPer : Math.max(20, Math.round(nPer * Math.abs(b - a) / (Math.abs(targets[1] - targets[0]) || 1)));
      let last = rev;
      for (let i = 0; i < n; i++) {
        const g = a + (b - a) * i / (n - 1);
        let t;
        if (Math.abs(g) > gmp) { gmp = Math.abs(g); t = fb(g); }
        else t = rev.t + 2 * fb((g - rev.g) / 2);
        G.push(g); T.push(t); S.push(k); last = { g, t };
      }
      rev = last;
    }
    return { g: G, t: T, seg: S };
  }

  /* ---------- HW3 資料處理（Q1、Q3；與 A1、C1 相同的演算法） ---------- */
  const HW = (() => {
    const d = window.CH5_HW3; if (!d) return null;
    const sp = d.spec, TAU = 2 * PI;
    const beta = bisect(b => b * Math.tan(b) - sp.IratioI0, 1e-9, PI / 2 - 1e-9);
    // Q1(a)：stage 1 掃頻，扭矩激振的加速度振幅 A = A0 r²/√[(1−r²)²+(2ξr)²]
    const sw = d.sweep;
    const shape = (q, fn, xi) => { const r = q / fn; return r * r / Math.sqrt((1 - r * r) ** 2 + (2 * xi * r) ** 2); };
    const a0of = (fn, xi) => { let sgy = 0, sgg = 0; sw.f.forEach((q, i) => { const v = shape(q, fn, xi); sgy += v * sw.a[i]; sgg += v * v; }); return sgy / sgg; };
    const sse = p => { if (p[1] <= 1e-5 || p[0] <= 0) return 1e99; const a0 = a0of(p[0], p[1]); let s = 0; sw.f.forEach((q, i) => { s += (sw.a[i] - a0 * shape(q, p[0], p[1])) ** 2; }); return s; };
    let ipk = 0; sw.a.forEach((v, i) => { if (v > sw.a[ipk]) ipk = i; });
    const fit = nelderMead(sse, [sw.f[ipk], 0.01], [0.4, 0.004]);
    const A0 = a0of(fit[0], fit[1]);
    const [x1, x2, x3] = [sw.f[ipk - 1], sw.f[ipk], sw.f[ipk + 1]], [y1, y2, y3] = [sw.a[ipk - 1], sw.a[ipk], sw.a[ipk + 1]];
    const den = (x1 - x2) * (x1 - x3) * (x2 - x3);
    const pa_ = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / den, pb_ = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / den;
    const fpar = -pb_ / (2 * pa_);
    const lvl = sw.a[ipk] / Math.SQRT2;
    let j = ipk; while (j > 0 && sw.a[j] >= lvl) j--;
    const f1 = sw.f[j] + (lvl - sw.a[j]) * (sw.f[j + 1] - sw.f[j]) / (sw.a[j + 1] - sw.a[j]);
    let k = ipk; while (k < sw.f.length - 1 && sw.a[k] >= lvl) k++;
    const f2 = sw.f[k - 1] + (lvl - sw.a[k - 1]) * (sw.f[k] - sw.f[k - 1]) / (sw.a[k] - sw.a[k - 1]);
    const q1a = { fpk: sw.f[ipk], apk: sw.a[ipk], fpar, fn: fit[0], xi: fit[1], A0, lvl, f1, f2, Dhp: (f2 - f1) / (2 * fit[0]) };
    // Q1(b)
    const ln = Math.log;
    function stages(fn1) {
      const rows = d.stages.map((s, i) => {
        const fn = i === 0 ? fn1 : s.fn, vs = TAU * fn * sp.H / beta, G = sp.rho * vs * vs;
        const del = ln(s.X[0] / s.X[5]) / 5;
        const lr = linreg([1, 2, 3, 4, 5, 6], s.X.map(v => ln(v)));
        const dr = DR(s.g);
        return { g: s.g, fn, vs, G, del, delr: -lr.slope, icpt: lr.icpt, D: 100 * del / TAU, Dx: 100 * del / Math.sqrt(TAU * TAU + del * del), drg: dr.gg, drD: dr.D, X: s.X };
      });
      const Gmax = rows[0].G;
      rows.forEach(r => { r.gg = r.G / Gmax; });
      const gs = rows.map(r => r.g), gg = rows.map(r => r.gg);
      let jj = gg.findIndex(v => v < 0.5);
      const grInt = Math.pow(10, Math.log10(gs[jj - 1]) + (0.5 - gg[jj - 1]) * (Math.log10(gs[jj]) - Math.log10(gs[jj - 1])) / (gg[jj] - gg[jj - 1]));
      const fitA = a => Math.pow(10, golden(L => gs.reduce((s, g, i) => s + (gg[i] - 1 / (1 + Math.pow(g / Math.pow(10, L), a))) ** 2, 0), -3, 0));
      return { rows, Gmax, grInt, grH: fitA(1), grM: fitA(PHI.a) };
    }
    const fit1 = stages(q1a.fn), tab1 = stages(d.stages[0].fn);
    // Q3 下孔
    const dh = d.downhole, z = [0, ...dh.z], t = [0, ...dh.t];
    const stress = zz => { const sv = dh.rho * 9.81 * zz, u = 9.81 * Math.max(zz - dh.gwt, 0), s1 = sv - u; return { sv, u, s1, sm: (1 + 2 * dh.K0) / 3 * s1 }; };
    return { d, beta, q1a, stages, fit1, tab1, z, t, stress, TAU };
  })();

  /* =====================================================================
     W1 應變範圍地圖
     ===================================================================== */
  (() => {
    const body = SD.widget('w-range'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'g', label: '剪應變 γ (%)', min: 1e-5, max: 5, log: true, snap: 2, value: 0.05, fmt: v => pctFmt(v) + ' %' }
    ], draw);
    SD.buttons(body, '例', [
      { label: '機器基礎 3×10⁻⁴ %', on: () => ctl.set({ g: 3e-4 }) },
      { label: '例：0.05 %', on: () => ctl.set({ g: 0.05 }) },
      { label: '強震軟弱層 0.3 %', on: () => ctl.set({ g: 0.3 }) }
    ]);
    const p = new Plot(body, { h: 250, logx: true, xlim: [1e-5, 10], ylim: [0, 1], yticks: [], m: { l: 16, r: 16, t: 10, b: 46 }, xlabel: '剪應變 γ (%)　［小數 = % ÷ 100］' });
    const st = SD.stats(body, [
      { k: 'pc', label: '應變（%）' }, { k: 'dec', label: '應變（小數）' }, { k: 'reg', label: '行為', hl: true },
      { k: 'mod', label: '模型' }, { k: 'ana', label: '分析方法', hl: true }, { k: 'thr', label: '相對砂的門檻應變' }
    ]);
    SD.div(body, 'w-note', '灰底三段是 Ishihara 的三個應變範圍（以 % 表示）。虛線是砂的線性門檻 γ<sub>tl</sub> ≈ 10<sup>−3</sup> % 與體積門檻 γ<sub>tv</sub> ≈ 10<sup>−2</sup> %；高塑性黏土的 γ<sub>tv</sub> 可到約 0.1 %。下方橫條是講義提到的典型問題，範圍為約略值。');
    function draw() {
      const g = ctl.g;
      p.clear();
      rectBand(p, 1e-5, 1e-3, 0, 1, 'f-blue', 0.08);
      rectBand(p, 1e-3, 1e-1, 0, 1, 'f-gray', 0.16);
      rectBand(p, 1e-1, 10, 0, 1, 'f-red', 0.10);
      p.text(1e-4, 0.88, '小應變：線彈性', { anchor: 'middle', cls: 'txt' });
      p.text(1e-2, 0.88, '中應變：彈塑性', { anchor: 'middle', cls: 'txt' });
      p.text(1, 0.88, '大應變：可能破壞', { anchor: 'middle', cls: 'txt' });
      p.vline(1e-3, { c: 'ink', dash: '4 3' }); p.vline(1e-2, { c: 'ink', dash: '4 3' });
      p.text(1e-3, 0.74, 'γ_{tl}', { anchor: 'middle', cls: 'txt-muted', dx: -14 });
      p.text(1e-2, 0.74, 'γ_{tv}（砂）', { anchor: 'middle', cls: 'txt-muted', dx: 30 });
      const bars = [[2e-5, 1e-3, 0.55, '機器基礎、交通振動'], [1e-3, 5, 0.38, '地震載重（中到大應變）'], [0.1, 5, 0.21, '強震下的軟弱層']];
      bars.forEach(([a, b, y, s]) => {
        p.line([a, b], [y, y], { c: 'gray', w: 7, op: 0.75 });
        p.text(Math.sqrt(a * b), y + 0.07, s, { anchor: 'middle', cls: 'txt-muted' });
      });
      p.vline(g, { c: 'blue', w: 2.2, dash: null });
      p.markers([g], [0.05], { c: 'blue', fill: true, r: 5 });
      const reg = g < 1e-3 ? ['近乎彈性', '線彈性', '線性分析（直接用 Gmax）'] : g <= 0.1 ? ['非線性、遲滯，但不破壞', '黏彈性（非退化）', '等值線性法'] : ['可能破壞；圈數、速率效應明顯', '追蹤載重歷史（退化）', '時間域逐步積分'];
      const thr = g < 1e-3 ? '< γ<sub>tl</sub>：G ≈ G<sub>max</sub>' : g < 1e-2 ? 'γ<sub>tl</sub>–γ<sub>tv</sub>：無永久體積變化' : '> γ<sub>tv</sub>：體積變化或孔隙水壓';
      st({ pc: sfmt(g), dec: decFmt(g / 100), reg: reg[0], mod: reg[1], ana: reg[2], thr });
    }
    draw();
  })();

  /* =====================================================================
     W2 一個迴圈、兩種勁度
     ===================================================================== */
  (() => {
    const body = SD.widget('w-loop'); if (!body) return;
    const Gmax = 94000, gr = 0.04, fb = hyp(Gmax, gr);
    const fp = g => Gmax / (1 + Math.abs(g) / gr) ** 2; // dτ/dγ (kPa per decimal)
    SD.div(body, 'w-note', '骨架：雙曲線，G<sub>max</sub> = 94 MPa、γ<sub>r</sub> = 0.04 %；迴圈由 Masing rules 產生（8.7）。');
    const ctl = SD.controls(body, [
      { k: 'gc', label: '應變幅度 γ_c (%)', min: 0.001, max: 1, log: true, snap: 2, value: 0.05, fmt: v => pctFmt(v) + ' %' },
      { k: 's', label: '迴圈上的位置', min: 0, max: 1, step: 0.005, value: 0.12, fmt: v => v < 0.5 ? '卸載 ' + (v * 200).toFixed(0) + ' %' : '再載入 ' + ((v - 0.5) * 200).toFixed(0) + ' %' }
    ], draw);
    SD.buttons(body, 'γ_c =', [0.003, 0.01, 0.05, 0.1, 0.5].map(v => ({ label: v + ' %', on: () => ctl.set({ gc: v }) })));
    const p = new Plot(body, { w: 600, h: 360, xlim: [-1, 1], ylim: [-1, 1], xlabel: '剪應變 γ (%)', ylabel: '剪應力 τ (kPa)' });
    p.legend([{ label: '骨架', c: 'gray', dash: true }, { label: '遲滯迴圈', c: 'blue', w: 2 }, { label: 'G_{sec}（原點到尖端）', c: 'red', dash: true }, { label: 'G_{tan}（目前位置的切線）', c: 'ink' }, { label: 'E_{So} 三角形', c: 'blue', line: false, marker: 'square', fill: true }]);
    const st = SD.stats(body, [
      { k: 'tc', label: 'τ_c (kPa)' }, { k: 'gs', label: 'G_{sec}/G_{max}', hl: true }, { k: 'gt', label: 'G_{tan}/G_{max}（目前位置）', hl: true },
      { k: 'ed', label: 'E_D 迴圈面積 (J/m³)' }, { k: 'es', label: 'E_{So} (J/m³)' }, { k: 'D', label: 'D = E_D/(4πE_{So})' }
    ]);
    function draw() {
      const gc = ctl.gc, tc = fb(gc), X = 1.3 * gc, Y = 1.25 * fb(X);
      p.setLim([-X, X], [-Y, Y]);
      p.clear();
      p.hline(0, { dash: null, w: 0.6 }); p.vline(0, { dash: null, w: 0.6 });
      const gu = linspace(gc, -gc, 200), gl = linspace(-gc, gc, 200);
      const tu = gu.map(g => tc + 2 * fb((g - gc) / 2)), tl = gl.map(g => -tc + 2 * fb((g + gc) / 2));
      fillPoly(p, gu.concat(gl), tu.concat(tl), 'f-band');
      fillPoly(p, [0, gc, gc], [0, 0, tc], 'f-blue', 0.22);
      p.fn(fb, { c: 'gray', dash: '6 4', w: 1.2, range: [-X, X] });
      p.line(gu, tu, { c: 'blue', w: 2 }); p.line(gl, tl, { c: 'blue', w: 2 });
      p.line([-X, X], [-X * tc / gc, X * tc / gc], { c: 'red', dash: '6 4', w: 1.4 });
      // 目前位置
      const s = ctl.s;
      let g, t, gt;
      if (s < 0.5) { g = gc - 4 * gc * s; t = tc + 2 * fb((g - gc) / 2); gt = fp((g - gc) / 2); }
      else { g = -gc + 4 * gc * (s - 0.5); t = -tc + 2 * fb((g + gc) / 2); gt = fp((g + gc) / 2); }
      const L = 0.35 * gc;
      p.line([g - L, g + L], [t - gt * L / 100, t + gt * L / 100], { c: 'ink', w: 1.6 });
      p.markers([gc, -gc], [tc, -tc], { c: 'ink', r: 3.6, over: true });
      p.markers([g], [t], { c: 'red', fill: true, over: true, r: 4.5 });
      p.text(gc, tc, 'A', { dx: 7, dy: -6 }); p.text(-gc, -tc, 'B', { dx: -16, dy: 14 });
      const x = gc / gr, Dm = masingA1(x);
      const eso = 0.5 * tc * gc / 100, ed = Dm * 4 * PI * eso;
      st({ tc: fx(tc, 2), gs: fx(tc / (gc / 100) / Gmax, 3), gt: fx(gt / Gmax, 3), ed: f(ed * 1000, 3), es: f(eso * 1000, 3), D: fx(Dm * 100, 2) + ' %' });
    }
    draw();
  })();

  /* =====================================================================
     W3 由迴圈面積算阻尼比
     ===================================================================== */
  (() => {
    const body = SD.widget('w-loopdamp'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'tc', label: 'τ_c (kPa)', type: 'number', value: 40 },
      { k: 'gc', label: 'γ_c (%)', type: 'number', value: 0.1 },
      { k: 'ed', label: '迴圈面積 E_D (kJ/m³)', type: 'number', value: 0.03 },
      { k: 'bad', label: '錯誤示範：γ 用 % 直接代入', type: 'check', value: false }
    ], draw);
    SD.buttons(body, '範例', [
      { label: '講義例：40 kPa、0.1 %', on: () => ctl.set({ tc: 40, gc: 0.1, ed: 0.03 }) },
      { label: '例：Masing γc = 0.1 %', on: () => ctl.set({ tc: 26.86, gc: 0.1, ed: 0.04264 }) },
      { label: '例：Masing γc = 0.01 %', on: () => ctl.set({ tc: 7.52, gc: 0.01, ed: 0.000223 }) }
    ]);
    const p = new Plot(body, { w: 600, h: 340, xlim: [-1, 1], ylim: [-1, 1], xlabel: '剪應變 γ (%)', ylabel: '剪應力 τ (kPa)' });
    p.legend([{ label: '等面積橢圓（線性黏彈）', c: 'blue', w: 2 }, { label: 'E_D（灰）', c: 'gray', line: false, marker: 'square', fill: true }, { label: 'E_{So} = ½τ_cγ_c（淡藍）', c: 'blue', line: false, marker: 'square', fill: true }, { label: 'G_{sec}', c: 'red', dash: true }]);
    const st = SD.stats(body, [
      { k: 'gd', label: 'γ_c（小數）' }, { k: 'es', label: 'E_{So} (kJ/m³)' }, { k: 'ed', label: 'E_D (kJ/m³)' },
      { k: 'D', label: 'D = E_D/(4πE_{So})', hl: true }, { k: 'gs', label: 'G_{sec} (MPa)' }
    ]);
    const note = SD.div(body, 'w-note');
    function draw() {
      const tc = Math.max(1e-6, ctl.tc), gc = Math.max(1e-6, ctl.gc), ed = Math.max(0, ctl.ed);
      const gd = gc / 100, eso = 0.5 * tc * gd, D = ed / (4 * PI * eso);
      const b = ed / (PI * gd); // 橢圓的 cos 振幅 (kPa)
      const th = linspace(0, 2 * PI, 241), xs = th.map(t => gc * Math.sin(t)), ys = th.map(t => tc * Math.sin(t) + b * Math.cos(t));
      const Y = 1.25 * Math.max(tc, Math.hypot(tc, b));
      p.setLim([-1.3 * gc, 1.3 * gc], [-Y, Y]);
      p.clear();
      p.hline(0, { dash: null, w: 0.6 }); p.vline(0, { dash: null, w: 0.6 });
      fillPoly(p, xs, ys, 'f-band');
      fillPoly(p, [0, gc, gc], [0, 0, tc], 'f-blue', 0.25);
      p.line(xs, ys, { c: 'blue', w: 2 });
      p.line([-1.3 * gc, 1.3 * gc], [-1.3 * tc, 1.3 * tc], { c: 'red', dash: '6 4', w: 1.3 });
      p.markers([gc], [tc], { c: 'ink', fill: true, over: true, r: 4 });
      const Dbad = ed / (4 * PI * 0.5 * tc * gc);
      st({ gd: decFmt(gd), es: f(eso, 4), ed: f(ed, 4), D: ctl.bad ? '<s>' + fx(D * 100, 2) + ' %</s> → ' + f(Dbad * 100, 3) + ' %' : fx(D * 100, 2) + ' %', gs: fx(tc / gd / 1000, 2) }, { D: ctl.bad });
      note.innerHTML = ctl.bad
        ? `錯誤示範：把 γ<sub>c</sub> = ${gc} 直接代入，E<sub>So</sub> = ${f(0.5 * tc * gc, 3)} 大了 100 倍，算出的阻尼比只有 ${f(Dbad * 100, 3)} %，比真值小 100 倍。`
        : '橢圓通過尖端 (γ<sub>c</sub>, τ<sub>c</sub>)，面積等於輸入的 E<sub>D</sub>：這就是等值線性法用來代表真實迴圈的線性黏彈性迴圈。分子是整個迴圈（灰），分母只有第一象限的三角形（淡藍）。';
    }
    draw();
  })();

  /* =====================================================================
     W4 MRD 曲線實驗室
     ===================================================================== */
  (() => {
    const body = SD.widget('w-mrd'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'PI', label: '塑性指數 PI', min: 0, max: 200, step: 1, value: 15, fmt: v => v.toFixed(0) },
      { k: 'OCR', label: '過壓密比 OCR', min: 1, max: 8, step: 0.1, value: 1, fmt: v => v.toFixed(1) },
      { k: 'sm', label: 'σ′_m (atm)', min: 0.25, max: 16, log: true, snap: 2, value: 1, fmt: v => v.toPrecision(2) },
      { k: 'N', label: '反覆次數 N', min: 1, max: 1000, log: true, snap: 2, value: 10, fmt: v => v.toPrecision(2) },
      { k: 'fq', label: '頻率 f (Hz)', min: 0.1, max: 100, log: true, snap: 2, value: 1, fmt: v => v.toPrecision(2) },
      { k: 'g', label: '讀值應變 γ (%)', min: 1e-4, max: 10, log: true, snap: 2, value: 0.1, fmt: v => pctFmt(v) + ' %' }
    ], draw);
    SD.buttons(body, '預設', [
      { label: '講義例：PI = 15、0.1 %', on: () => ctl.set({ PI: 15, OCR: 1, sm: 1, N: 10, fq: 1, g: 0.1 }) },
      { label: '例：砂 PI = 0、0.05 %', on: () => ctl.set({ PI: 0, OCR: 1, sm: 1, N: 10, fq: 1, g: 0.05 }) },
      { label: '墨西哥城黏土：PI = 200', on: () => ctl.set({ PI: 200, OCR: 1, sm: 1, N: 10, fq: 1, g: 0.1 }) }
    ]);
    const row = SD.div(body, 'plot-row');
    const xl = [1e-4, 10];
    const pG = new Plot(row, { w: 420, h: 300, logx: true, xlim: xl, ylim: [0, 1.05], xlabel: '剪應變 γ (%)', ylabel: 'G/G_{max}' });
    const pD = new Plot(row, { w: 420, h: 300, logx: true, xlim: xl, ylim: [0, 30], xlabel: '剪應變 γ (%)', ylabel: '阻尼比 D (%)' });
    pG.legend([{ label: '目前參數', c: 'blue', w: 2 }, { label: '參考：PI = 0、1 atm', c: 'gray' }, { label: 'γ_{tl}–γ_{tv}（砂）', c: 'gray', line: false, marker: 'square', fill: true }]);
    const st = SD.stats(body, [
      { k: 'gr', label: '參考應變 γ_r (%)' }, { k: 'dmin', label: 'D_{min} (%)' }, { k: 'gg', label: '此應變 G/G_{max}', hl: true },
      { k: 'D', label: '此應變 D (%)', hl: true }, { k: 'vs', label: 'V_{s,eq}/V_s = √(G/G_{max})' }, { k: 'reg', label: '範圍（以砂的門檻）' }
    ]);
    const gs = logspace(xl[0], xl[1], 300);
    const ref = gs.map(g => DR(g));
    pG.readout(g => { const r = DR(g, opt()); return `γ = ${pctFmt(g)} %　G/G<sub>max</sub> = ${r.gg.toFixed(3)}`; });
    pD.readout(g => { const r = DR(g, opt()); return `γ = ${pctFmt(g)} %　D = ${r.D.toFixed(2)} %`; });
    function opt() { return { PI: ctl.PI, OCR: ctl.OCR, sm: ctl.sm, N: ctl.N, f: ctl.fq }; }
    function draw() {
      const o = opt(), cur = gs.map(g => DR(g, o)), r = DR(ctl.g, o);
      [pG, pD].forEach(p => { p.clear(); p.band(1e-3, 1e-2); p.vline(ctl.g, { c: 'blue', dash: '4 3' }); });
      pG.line(gs, ref.map(q => q.gg), { c: 'gray', w: 1.2 }); pD.line(gs, ref.map(q => q.D), { c: 'gray', w: 1.2 });
      pG.line(gs, cur.map(q => q.gg), { c: 'blue', w: 2 }); pD.line(gs, cur.map(q => q.D), { c: 'blue', w: 2 });
      pG.markers([ctl.g], [r.gg], { c: 'blue', fill: true, over: true, r: 4.5 });
      pD.markers([ctl.g], [r.D], { c: 'blue', fill: true, over: true, r: 4.5 });
      pG.markers([r.gr], [0.5], { c: 'red', fill: true, over: true, r: 3.6 });
      pG.text(r.gr, 0.5, 'γ_r', { dx: 7, dy: -6, cls: 'txt-red' });
      const reg = ctl.g < 1e-3 ? '< γ<sub>tl</sub>：線性' : ctl.g < 1e-2 ? 'γ<sub>tl</sub>–γ<sub>tv</sub>：等值線性' : '> γ<sub>tv</sub>：體積變化／孔隙水壓';
      st({ gr: fx(r.gr, 4), dmin: fx(r.dmin, 2), gg: fx(r.gg, 3), D: fx(r.D, 2), vs: fx(Math.sqrt(r.gg), 3), reg });
    }
    draw();
  })();

  /* =====================================================================
     W5 骨架曲線
     ===================================================================== */
  (() => {
    const body = SD.widget('w-backbone'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'gr', label: '參考應變 γ_r (%)', min: 0.005, max: 1, log: true, snap: 2, value: 0.04, fmt: v => pctFmt(v) + ' %' },
      { k: 'a', label: '曲率 a', min: 0.6, max: 1.3, step: 0.001, value: 0.919, fmt: v => v.toFixed(3) },
      { k: 'G', label: 'G_max (MPa)', min: 10, max: 300, step: 1, value: 94, fmt: v => v.toFixed(0) }
    ], draw);
    SD.buttons(body, '預設', [
      { label: '雙曲線 a = 1', on: () => ctl.set({ a: 1 }) },
      { label: 'Darendeli a = 0.919', on: () => ctl.set({ a: 0.919 }) },
      { label: '例：94 MPa、0.04 %', on: () => ctl.set({ G: 94, gr: 0.04, a: 1 }) },
      { label: '講義：黏土 50 MPa、強度推 γr = 0.1 %', on: () => ctl.set({ G: 50, gr: 0.1, a: 1 }) }
    ]);
    const row = SD.div(body, 'plot-row');
    const pT = new Plot(row, { w: 420, h: 300, xlim: [0, 1], ylim: [0, 1], xlabel: '剪應變 γ (%)', ylabel: '剪應力 τ (kPa)' });
    const pR = new Plot(row, { w: 420, h: 300, logx: true, xlim: [0.01, 100], ylim: [0, 1.05], xlabel: 'γ/γ_r', ylabel: 'G/G_{max}' });
    pT.legend([{ label: '目前 a', c: 'blue', w: 2 }, { label: '雙曲線 a = 1', c: 'gray' }, { label: 'τ_{max} = G_{max}γ_r', c: 'gray', dash: true }, { label: '初始斜率 G_{max}', c: 'ink', dash: '2 3' }]);
    const st = SD.stats(body, [
      { k: 'tm', label: 'τ_{max} = G_{max}γ_r (kPa)' }, { k: 'g01', label: 'G/G_{max} @ 0.1γ_r' }, { k: 'g1', label: 'G/G_{max} @ γ_r', hl: true },
      { k: 'g3', label: 'G/G_{max} @ 3γ_r' }, { k: 'g10', label: 'G/G_{max} @ 10γ_r' }, { k: 't10', label: 'τ(10γ_r)/τ_{max}' }
    ]);
    SD.div(body, 'w-note', '右圖橫軸以 γ/γ<sub>r</sub> 正規化：任何 γ<sub>r</sub> 都疊在同一條 S 形上，並以 γ = γ<sub>r</sub>（G/G<sub>max</sub> = ½）為中心對稱。a &lt; 1 時大應變段較平、小應變段較早開始下降；a &lt; 1 的骨架沒有水平漸近線，τ 會慢慢超過 G<sub>max</sub>γ<sub>r</sub>。');
    function draw() {
      const G = ctl.G * 1000, gr = ctl.gr, a = ctl.a;
      const tau = (g, aa) => G * (g / 100) / (1 + Math.pow(g / gr, aa));
      const X = 10 * gr, tm = G * gr / 100;
      const ymax = 1.15 * Math.max(tm, tau(X, a));
      pT.setLim([0, X], [0, ymax]); pT.clear(); pR.clear();
      pT.hline(tm, { dash: '6 4' });
      pT.line([0, tm / (G / 100)], [0, tm], { c: 'ink', dash: '2 3', w: 1.2 });
      pT.fn(g => tau(g, 1), { c: 'gray', w: 1.2 });
      pT.fn(g => tau(g, a), { c: 'blue', w: 2 });
      pT.markers([gr], [tau(gr, a)], { c: 'red', fill: true, over: true });
      pT.text(gr, tau(gr, a), 'γ = γ_r', { dx: 8, dy: 16, cls: 'txt-red' });
      pR.vline(1, { dash: '4 3' }); pR.hline(0.5, { dash: '4 3' });
      pR.fn(x => 1 / (1 + x), { c: 'gray', w: 1.2 });
      pR.fn(x => 1 / (1 + Math.pow(x, a)), { c: 'blue', w: 2 });
      pR.markers([0.1, 1, 10], [0.1, 1, 10].map(x => 1 / (1 + Math.pow(x, a))), { c: 'blue', fill: true, over: true, r: 3.8 });
      const gg = x => 1 / (1 + Math.pow(x, a));
      st({ tm: fx(tm, 1), g01: fx(gg(0.1), 3), g1: fx(gg(1), 3), g3: fx(gg(3), 3), g10: fx(gg(10), 3), t10: fx(tau(X, a) / tm, 3) });
    }
    draw();
  })();

  /* =====================================================================
     W6 Masing 迴圈產生器
     ===================================================================== */
  (() => {
    const body = SD.widget('w-masing'); if (!body) return;
    const rnd = SD.rng(20261105);
    const eqPath = (() => { const out = [0]; let s = 1; for (let i = 0; i < 18; i++) { const env = 0.12 * Math.exp(-Math.pow((i - 5) / 6, 2)); out.push(s * env * (0.35 + 0.65 * rnd())); s = -s; } out.push(0); return out.map(v => +v.toFixed(4)); })();
    const PATHS = {
      hw: { n: '例：0 → γc → −γc → γc（γc = 0.1 %）', p: [0, 0.1, -0.1, 0.1] },
      grow: { n: '漸增振幅', p: [0, 0.02, -0.04, 0.06, -0.08, 0.1, -0.1, 0.1] },
      inner: { n: '內迴圈：0.1 → −0.02 → 0.05 → −0.1 …', p: [0, 0.1, -0.02, 0.05, -0.1, 0.03, -0.04, 0.1] },
      eq: { n: '不規則（合成的地震應變歷時）', p: eqPath }
    };
    const ctl = SD.controls(body, [
      { k: 'path', label: '應變路徑', type: 'select', value: 'inner', options: Object.keys(PATHS).map(k => [k, PATHS[k].n]) },
      { k: 'prog', label: '播放進度', min: 0.02, max: 1, step: 0.005, value: 1, fmt: v => (v * 100).toFixed(0) + ' %' },
      { k: 'cmp', label: '同時畫 HW3 網頁工具的簡化做法', type: 'check', value: true }
    ], draw);
    let playing = false, raf = null;
    SD.buttons(body, null, [{ label: '▶ 從頭播放', primary: true, on: () => { ctl.set({ prog: 0.02 }); playing = true; step(); } }]);
    const fb = hyp(94000, 0.04);
    const row = SD.div(body, 'plot-row');
    const p = new Plot(row, { w: 420, h: 330, xlim: [-0.12, 0.12], ylim: [-40, 40], xlabel: '剪應變 γ (%)', ylabel: '剪應力 τ (kPa)' });
    const pg = new Plot(row, { w: 420, h: 330, xlim: [0, 1], ylim: [-0.12, 0.12], xlabel: '載重步驟（正規化）', ylabel: '剪應變 γ (%)' });
    p.legend([{ label: '骨架', c: 'gray', dash: true }, { label: '擴充 Masing rules', c: 'blue', w: 2 }, { label: '簡化（只記最後反轉點）', c: 'red', dash: true }, { label: '反轉點', c: 'ink', line: false, marker: 'circle' }]);
    const st = SD.stats(body, [
      { k: 'nr', label: '反轉點數' }, { k: 'open', label: '尚未閉合的反轉點（堆疊深度）' }, { k: 'g', label: '目前 γ (%)' },
      { k: 'te', label: '目前 τ：擴充 (kPa)', hl: true }, { k: 'ts', label: '目前 τ：簡化 (kPa)' }, { k: 'dm', label: '兩者最大差 |Δτ| (kPa)' }
    ]);
    SD.div(body, 'w-note', '骨架與 W2 相同（94 MPa、γ<sub>r</sub> = 0.04 %）。「內迴圈」路徑在 −0.02 % 反轉、上到 0.05 % 再往下：越過 −0.02 % 時內迴圈閉合，擴充規則改沿 A 出發的外分支走，到 −0.1 % 正好回到骨架；簡化做法則沿內分支一路往下，直到超過歷史最大應變才跳回骨架。合成路徑由固定亂數種子產生，不是實測紀錄。');
    function step() {
      if (!playing) return;
      const v = Math.min(1, ctl.prog + 0.006);
      ctl.set({ prog: v });
      if (v >= 1) { playing = false; return; }
      raf = requestAnimationFrame(step);
    }
    function draw(k) {
      if (k === 'prog' && playing === true && false) return;
      if (k && k !== 'prog' && k !== '*') playing = false;
      const P = PATHS[ctl.path].p;
      const ext = masingExt(P, fb, 160);
      const sim = masingSimple(P, fb, 160, false);
      const gm = Math.max(...P.map(Math.abs)) * 1.2, tm = Math.max(...ext.t.map(Math.abs), ...sim.t.map(Math.abs)) * 1.15;
      p.setLim([-gm, gm], [-tm, tm]); pg.setLim(null, [-gm, gm]);
      p.clear(); pg.clear();
      p.hline(0, { dash: null, w: 0.6 }); p.vline(0, { dash: null, w: 0.6 });
      p.fn(fb, { c: 'gray', dash: '6 4', w: 1.1, range: [-gm, gm] });
      const nE = Math.max(2, Math.round(ctl.prog * ext.g.length)), nS = Math.max(2, Math.round(ctl.prog * sim.g.length));
      if (ctl.cmp) p.line(sim.g.slice(0, nS), sim.t.slice(0, nS), { c: 'red', dash: '6 4', w: 1.5 });
      p.line(ext.g.slice(0, nE), ext.t.slice(0, nE), { c: 'blue', w: 2 });
      // 反轉點：在擴充路徑上找方向改變處
      const rg = [], rt = [];
      for (let i = 1; i < nE - 1; i++) if ((ext.g[i] - ext.g[i - 1]) * (ext.g[i + 1] - ext.g[i]) < 0) { rg.push(ext.g[i]); rt.push(ext.t[i]); }
      p.markers(rg, rt, { c: 'ink', r: 3.6, over: true });
      p.markers([ext.g[nE - 1]], [ext.t[nE - 1]], { c: 'blue', fill: true, r: 4.5, over: true });
      if (ctl.cmp) p.markers([sim.g[nS - 1]], [sim.t[nS - 1]], { c: 'red', fill: true, r: 3.6, over: true });
      const idx = P.map((_, i) => i / (P.length - 1));
      pg.hline(0, { dash: null, w: 0.6 });
      pg.line(idx, P, { c: 'ink', w: 1.4 });
      pg.markers(idx, P, { c: 'ink', r: 3.4 });
      pg.vline(ctl.prog, { c: 'blue', dash: '4 3' });
      // 簡化 vs 擴充：以 γ 序列的位置比較（兩者以相同比例取點）
      let dm = 0;
      const M = Math.min(ext.g.length, sim.g.length);
      for (let i = 0; i < 400; i++) {
        const q = i / 399, a = Math.round(q * (ext.g.length - 1)), b = Math.round(q * (sim.g.length - 1));
        if (Math.abs(ext.g[a] - sim.g[b]) < 0.002 * gm) dm = Math.max(dm, Math.abs(ext.t[a] - sim.t[b]));
      }
      const openNow = masingExt(P.slice(0, 1 + Math.ceil(ctl.prog * (P.length - 1) - 1e-9)), fb, 40).open;
      st({ nr: String(rg.length), open: String(openNow), g: fx(ext.g[nE - 1], 4), te: fx(ext.t[nE - 1], 2), ts: ctl.cmp ? fx(sim.t[nS - 1], 2) : '—', dm: ctl.cmp ? fx(dm, 2) : '—' }, { dm: ctl.cmp && dm > 0.5 });
      void M;
    }
    draw();
  })();

  /* =====================================================================
     W7 Masing 阻尼 vs 修正後的阻尼
     ===================================================================== */
  (() => {
    const body = SD.widget('w-masingD'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'PI', label: '塑性指數 PI', min: 0, max: 100, step: 1, value: 0, fmt: v => v.toFixed(0) },
      { k: 'sm', label: 'σ′_m (atm)', min: 0.25, max: 16, log: true, snap: 2, value: 1, fmt: v => v.toPrecision(2) },
      { k: 'N', label: '反覆次數 N', min: 1, max: 1000, log: true, snap: 2, value: 10, fmt: v => v.toPrecision(2) },
      { k: 'g', label: '讀值應變 γ (%)', min: 1e-4, max: 10, log: true, snap: 2, value: 0.1, fmt: v => pctFmt(v) + ' %' },
      { k: 'rc', label: '加上共振柱實測資料（PI = 0、1 atm）', type: 'check', value: true }
    ], draw);
    const p = new Plot(body, { w: 640, h: 340, logx: true, xlim: [1e-4, 10], ylim: [0, 70], xlabel: '剪應變 γ (%)', ylabel: '阻尼比 D (%)' });
    p.legend([{ label: 'Masing，雙曲線 a = 1', c: 'gray', dash: true }, { label: 'Masing，修正雙曲線 a = 0.919', c: 'red', w: 1.8 }, { label: 'Darendeli：b(G/G_{max})^{0.1}D_{Masing} + D_{min}', c: 'blue', w: 2 }, { label: 'D_{min}', c: 'gray', dash: '2 3' }, { label: '共振柱實測', c: 'ink', dash: true, marker: 'circle' }]);
    const st = SD.stats(body, [
      { k: 'x', label: 'x = γ/γ_r' }, { k: 'dm', label: 'D_{Masing}（a = 0.919，數值積分）' }, { k: 'dmc', label: 'D_{Masing}（Darendeli 三次近似）' },
      { k: 'b', label: 'b = 0.6329 − 0.0057 ln N' }, { k: 'gg', label: '(G/G_{max})^{0.1}' }, { k: 'dmin', label: 'D_{min} (%)' },
      { k: 'D', label: '修正後 D (%)', hl: true }, { k: 'ratio', label: 'Masing / 修正後', hl: true }
    ]);
    const gs = logspace(1e-4, 10, 220);
    function draw() {
      const o = { PI: ctl.PI, sm: ctl.sm, N: ctl.N }, gr = DR(1, o).gr;
      p.clear();
      p.hline(200 / PI, { dash: '4 3' });
      p.text(2e-4, 200 / PI, '2/π ≈ 63.7 %（a = 1 的上限）', { dy: -6, cls: 'txt-muted' });
      p.line(gs, gs.map(g => 100 * masingA1(g / gr)), { c: 'gray', dash: '6 4', w: 1.3 });
      p.line(gs, gs.map(g => 100 * masingNum(g / gr, PHI.a)), { c: 'red', w: 1.8 });
      p.line(gs, gs.map(g => DR(g, o).D), { c: 'blue', w: 2 });
      p.hline(DR(1e-6, o).dmin, { c: 'gray', dash: '2 3' });
      if (ctl.rc && HW) {
        const r = HW.fit1.rows;
        p.line(r.map(q => q.g), r.map(q => q.D), { c: 'ink', dash: '5 4', w: 1 });
        p.markers(r.map(q => q.g), r.map(q => q.D), { c: 'ink', r: 3.8 });
      }
      p.vline(ctl.g, { c: 'blue', dash: '4 3' });
      const r = DR(ctl.g, o), dmN = 100 * masingNum(ctl.g / gr, PHI.a);
      p.markers([ctl.g, ctl.g], [dmN, r.D], { c: 'ink', fill: true, over: true, r: 3.6 });
      st({ x: f(ctl.g / gr, 3), dm: fx(dmN, 2) + ' %', dmc: fx(r.dm, 2) + ' %', b: fx(r.b, 4), gg: fx(Math.pow(r.gg, 0.1), 4), dmin: fx(r.dmin, 3), D: fx(r.D, 2), ratio: fx(dmN / r.D, 2) });
    }
    draw();
  })();

  /* =====================================================================
     W8 影響因素
     ===================================================================== */
  (() => {
    const body = SD.widget('w-factors'); if (!body) return;
    const panes = SD.tabs(body, ['塑性 PI', '圍壓 σ′m', '礫石級配 Cu（Menq 2003）']);
    const gs = logspace(1e-4, 10, 220);
    function famPane(pane, defs, famVals, famKey, label, fixedOpt) {
      const ctl = SD.controls(pane, defs, draw);
      const row = SD.div(pane, 'plot-row');
      const pG = new Plot(row, { w: 420, h: 290, logx: true, xlim: [1e-4, 10], ylim: [0, 1.05], xlabel: '剪應變 γ (%)', ylabel: 'G/G_{max}' });
      const pD = new Plot(row, { w: 420, h: 290, logx: true, xlim: [1e-4, 10], ylim: [0, 30], xlabel: '剪應變 γ (%)', ylabel: '阻尼比 D (%)' });
      pG.legend([{ label: label + ' = ' + famVals.join('、'), c: 'gray' }, { label: '選定的一條', c: 'blue', w: 2.2 }]);
      const st = SD.stats(pane, [{ k: 'gr', label: '選定曲線 γ_r (%)' }, { k: 'gg', label: 'G/G_{max} @ 0.1 %', hl: true }, { k: 'D', label: 'D @ 0.1 % (%)', hl: true }, { k: 'dmin', label: 'D_{min} (%)' }]);
      pG.readout(g => { const r = DR(g, opt(ctl[famKey])); return `γ = ${pctFmt(g)} %　選定曲線 G/G<sub>max</sub> = ${r.gg.toFixed(3)}`; });
      pD.readout(g => { const r = DR(g, opt(ctl[famKey])); return `γ = ${pctFmt(g)} %　選定曲線 D = ${r.D.toFixed(2)} %`; });
      function opt(v) { const o = fixedOpt(ctl); o[famKey] = +v; return o; }
      function draw() {
        pG.clear(); pD.clear();
        famVals.forEach(v => {
          const c = gs.map(g => DR(g, opt(v)));
          pG.line(gs, c.map(q => q.gg), { c: 'gray', w: 1 }); pD.line(gs, c.map(q => q.D), { c: 'gray', w: 1 });
        });
        const c = gs.map(g => DR(g, opt(ctl[famKey])));
        pG.line(gs, c.map(q => q.gg), { c: 'blue', w: 2.2 }); pD.line(gs, c.map(q => q.D), { c: 'blue', w: 2.2 });
        const r = DR(0.1, opt(ctl[famKey]));
        pG.vline(0.1, { dash: '4 3' }); pD.vline(0.1, { dash: '4 3' });
        st({ gr: fx(r.gr, 4), gg: fx(r.gg, 3), D: fx(r.D, 2), dmin: fx(r.dmin, 2) });
      }
      draw();
    }
    famPane(panes[0], [
      { k: 'PI', label: '選定 PI', type: 'select', value: '0', options: [0, 15, 30, 50, 100, 200].map(v => [String(v), 'PI = ' + v]) },
      { k: 'sm', label: 'σ′_m (atm)', min: 0.25, max: 16, log: true, snap: 2, value: 1, fmt: v => v.toPrecision(2) }
    ], [0, 15, 30, 50, 100, 200], 'PI', 'PI', c => ({ sm: c.sm }));
    famPane(panes[1], [
      { k: 'sm', label: '選定 σ′m', type: 'select', value: '1', options: [0.25, 1, 4, 16].map(v => [String(v), v + ' atm']) },
      { k: 'PI', label: 'PI', min: 0, max: 100, step: 1, value: 0, fmt: v => v.toFixed(0) }
    ], [0.25, 1, 4, 16], 'sm', 'σ′_m (atm)', c => ({ PI: c.PI }));

    // 礫石
    (() => {
      const pane = panes[2];
      const ctl = SD.controls(pane, [
        { k: 'cu', label: '均勻係數 C_u', min: 1.2, max: 50, log: true, snap: 2, value: 2, fmt: v => v.toPrecision(2) },
        { k: 'sm', label: 'σ′_m (atm)', min: 0.25, max: 16, log: true, snap: 2, value: 1, fmt: v => v.toPrecision(2) }
      ], draw);
      SD.buttons(pane, '講義例', [{ label: 'C_u = 2（均勻級配）', on: () => ctl.set({ cu: 2, sm: 1 }) }, { label: 'C_u = 20（優良級配）', on: () => ctl.set({ cu: 20, sm: 1 }) }]);
      const p = new Plot(pane, { w: 640, h: 300, logx: true, xlim: [1e-4, 10], ylim: [0, 1.05], xlabel: '剪應變 γ (%)', ylabel: 'G/G_{max}' });
      p.legend([{ label: '礫石（Menq 的 γ_r；形狀借用 a = 0.919，示意）', c: 'blue', w: 2.2 }, { label: '乾淨砂 Darendeli（同 σ′_m）', c: 'gray' }]);
      const st = SD.stats(pane, [{ k: 'gr', label: 'γ_r = 0.12C_u^{−0.6}(σ′_m/p_a)^{n} (%)', hl: true }, { k: 'n', label: '圍壓指數 n = 0.5C_u^{−0.15}' }, { k: 'grs', label: '同 σ′m 乾淨砂 γ_r (%)' }, { k: 'r', label: '礫石 / 砂' }]);
      SD.div(pane, 'w-note', '講義只給 Menq (2003) 的參考應變式；這裡為了畫出整條曲線，形狀暫借 Darendeli 的平均曲率 a = 0.919，所以曲線形狀只是示意，γ<sub>r</sub> 的數值才是重點：C<sub>u</sub> 越大 γ<sub>r</sub> 越小、越早進入非線性，而且圍壓的影響變弱。');
      function draw() {
        const n = 0.5 * Math.pow(ctl.cu, -0.15), gr = 0.12 * Math.pow(ctl.cu, -0.6) * Math.pow(ctl.sm, n), grs = DR(1, { sm: ctl.sm }).gr;
        p.clear();
        p.line(gs, gs.map(g => DR(g, { sm: ctl.sm }).gg), { c: 'gray', w: 1.2 });
        p.line(gs, gs.map(g => 1 / (1 + Math.pow(g / gr, PHI.a))), { c: 'blue', w: 2.2 });
        p.markers([gr], [0.5], { c: 'blue', fill: true, over: true });
        p.markers([grs], [0.5], { c: 'gray', fill: true, over: true, r: 3.4 });
        st({ gr: fx(gr, 4), n: fx(n, 3), grs: fx(grs, 4), r: fx(gr / grs, 2) });
      }
      draw();
    })();
  })();

  /* =====================================================================
     W9 Gmax 計算器
     ===================================================================== */
  (() => {
    const body = SD.widget('w-gmax'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'mode', label: 'σ′m 的來源', type: 'select', value: 'depth', options: [['depth', '由深度、地下水位、K0 計算'], ['direct', '直接輸入 σ′m']] },
      { k: 'z', label: '深度 z (m)', min: 0.5, max: 40, step: 0.5, value: 15, fmt: v => v.toFixed(1) },
      { k: 'gwt', label: '地下水位 (m)', min: 0, max: 10, step: 0.5, value: 1, fmt: v => v.toFixed(1) },
      { k: 'rho', label: '密度 ρ (t/m³)', min: 1.5, max: 2.3, step: 0.01, value: 1.95, fmt: v => v.toFixed(2) },
      { k: 'K0', label: 'K_0', min: 0.3, max: 2, step: 0.05, value: 0.5, fmt: v => v.toFixed(2) },
      { k: 'smd', label: '直接輸入 σ′m (kPa)', type: 'number', value: 100 },
      { k: 'e', label: '孔隙比 e', min: 0.35, max: 1.3, step: 0.01, value: 0.7, fmt: v => v.toFixed(2) },
      { k: 'OCR', label: 'OCR', min: 1, max: 8, step: 0.1, value: 1, fmt: v => v.toFixed(1) },
      { k: 'PI', label: 'PI（決定 k）', min: 0, max: 100, step: 1, value: 0, fmt: v => v.toFixed(0) },
      { k: 'N1', label: '(N_1)_{60}（Seed 式）', min: 2, max: 60, step: 1, value: 15, fmt: v => v.toFixed(0) },
      { k: 'pa', label: 'p_a', type: 'select', value: '101.325', options: [['101.325', '101.325 kPa'], ['100', '100 kPa']] }
    ], draw);
    SD.buttons(body, '範例', [
      { label: '講義例：e = 0.6、σ′m = 100 kPa', on: () => ctl.set({ mode: 'direct', smd: 100, e: 0.6, OCR: 1, PI: 0, rho: 2.0, N1: 15, pa: '100' }) },
      { label: '例：深 15 m、e = 0.70', on: () => ctl.set({ mode: 'depth', z: 15, gwt: 1, rho: 1.95, K0: 0.5, e: 0.7, OCR: 1, PI: 0, pa: '101.325' }) }
    ]);
    const p = new Plot(body, { w: 640, h: 300, xlim: [0, 400], ylim: [0, 250], xlabel: '平均有效應力 σ′_m (kPa)', ylabel: 'G_{max} (MPa)' });
    p.legend([{ label: 'Hardin (1978)', c: 'blue', w: 2 }, { label: 'Hardin & Drnevich (1972)', c: 'gray' }, { label: 'Seed et al. (1986)', c: 'red', dash: true }, { label: '目前 σ′_m', c: 'blue', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [
      { k: 'sv', label: 'σ′_v (kPa)' }, { k: 'sm', label: 'σ′_m (kPa)', hl: true }, { k: 'Fe', label: 'F(e) = 1/(0.3+0.7e²)' }, { k: 'k', label: 'OCR 指數 k' },
      { k: 'h78', label: 'Hardin 1978 G_{max} (MPa)', hl: true }, { k: 'vh', label: '對應 V_s (m/s)' },
      { k: 'hd', label: 'Hardin–Drnevich G_{max} (MPa)' }, { k: 'sd', label: 'Seed (N_1)_{60} G_{max} (MPa)' }
    ]);
    const kOfPI = pi => { const P = [0, 20, 40, 60, 80, 100], K = [0, 0.18, 0.30, 0.41, 0.48, 0.5]; if (pi >= 100) return 0.5; for (let i = 1; i < P.length; i++) if (pi <= P[i]) return K[i - 1] + (K[i] - K[i - 1]) * (pi - P[i - 1]) / (P[i] - P[i - 1]); return 0.5; };
    SD.div(body, 'w-note', 'k 依 PI 線性內插（PI = 0、20、40、60、80、≥100 時 0、0.18、0.30、0.41、0.48、0.50）。Seed 式內部先把 kPa 換成 psf（1 psf = 0.04788 kPa）再換回 MPa。三條經驗式的差別就是「經驗式離散」的縮影，現地量測才是依據。');
    function draw() {
      const pa = +ctl.pa, F = 1 / (0.3 + 0.7 * ctl.e * ctl.e), k = kOfPI(ctl.PI), oc = Math.pow(ctl.OCR, k);
      let sv1 = NaN, sm;
      if (ctl.mode === 'depth') {
        const sv = ctl.rho * 9.81 * ctl.z, u = 9.81 * Math.max(ctl.z - ctl.gwt, 0);
        sv1 = sv - u; sm = (1 + 2 * ctl.K0) / 3 * sv1;
      } else sm = Math.max(1, ctl.smd);
      const h78 = s => 625 * F * oc * Math.sqrt(pa * s) / 1000;
      const hd = s => 3230 * (2.97 - ctl.e) ** 2 / (1 + ctl.e) * oc * Math.sqrt(s) / 1000;
      const K2 = 20 * Math.cbrt(ctl.N1), psf = 0.04788026;
      const sd = s => 1000 * K2 * Math.sqrt(s / psf) * psf / 1000;
      const ymax = 1.15 * Math.max(h78(400), hd(400), sd(400));
      p.setLim(null, [0, ymax]); p.clear();
      p.fn(hd, { c: 'gray', w: 1.3 }); p.fn(sd, { c: 'red', dash: '6 4', w: 1.5 }); p.fn(h78, { c: 'blue', w: 2 });
      p.vline(sm);
      p.markers([sm], [h78(sm)], { c: 'blue', fill: true, over: true });
      const G = h78(sm);
      st({ sv: fx(sv1, 1), sm: fx(sm, 1), Fe: fx(F, 4), k: fx(k, 2), h78: fx(G, 1), vh: fx(Math.sqrt(G * 1000 / ctl.rho), 0), hd: fx(hd(sm), 1), sd: fx(sd(sm), 1) });
    }
    draw();
  })();

  /* =====================================================================
     W10 下孔試驗：斜射路徑
     ===================================================================== */
  (() => {
    const body = SD.widget('w-downhole'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'v1', label: '上層 V_{s1} (m/s)', min: 80, max: 400, step: 5, value: 150, fmt: v => v.toFixed(0) },
      { k: 'h1', label: '上層厚度 (m)', min: 2, max: 15, step: 0.5, value: 6, fmt: v => v.toFixed(1) },
      { k: 'v2', label: '下層 V_{s2} (m/s)', min: 100, max: 600, step: 5, value: 300, fmt: v => v.toFixed(0) },
      { k: 'x', label: '震源離孔口 x (m)', min: 0, max: 10, step: 0.25, value: 3, fmt: v => v.toFixed(2) },
      { k: 'cor', label: '做斜距修正 t_v = t·z/√(z²+x²)', type: 'check', value: true }
    ], draw);
    SD.buttons(body, '範例', [{ label: '講義例：150／300 m/s、6 m、震源緊鄰', on: () => ctl.set({ v1: 150, h1: 6, v2: 300, x: 0 }) }, { label: '同上、震源離孔 5 m', on: () => ctl.set({ v1: 150, h1: 6, v2: 300, x: 5 }) }]);
    const row = SD.div(body, 'plot-row');
    const pt = new Plot(row, { w: 420, h: 300, xlim: [0, 20], ylim: [0, 100], xlabel: '受波器深度 z (m)', ylabel: '走時 t (ms)' });
    const pv = new Plot(row, { w: 420, h: 300, xlim: [0, 20], ylim: [0, 400], xlabel: '深度 z (m)', ylabel: '區間 V_s (m/s)' });
    pt.legend([{ label: '真實垂直走時', c: 'gray' }, { label: '合成量測（斜射）', c: 'ink', dash: true, marker: 'circle' }, { label: '修正後', c: 'blue', marker: 'circle', fill: true, line: false }]);
    pv.legend([{ label: '真實 V_s', c: 'gray', w: 2 }, { label: '由量測直接算', c: 'ink', dash: true, marker: 'circle' }, { label: '修正後', c: 'blue', marker: 'circle', fill: true, line: false }]);
    const st = SD.stats(body, [
      { k: 't12', label: '12 m 的垂直走時 (ms)' }, { k: 'raw1', label: '0–1 m 區間 V_s（未修正）' }, { k: 'cor1', label: '0–1 m 區間 V_s（修正）' },
      { k: 'reg', label: '上層迴歸 V_s（未修正）' }, { k: 'regc', label: '上層迴歸 V_s（修正）', hl: true }
    ]);
    SD.div(body, 'w-note', '合成資料：受波器每 1 m 一個，走時用直線射線近似（忽略層界折射），只為示範斜射造成的偏差。震源越遠、深度越淺，未修正的走時越長，淺層 V<sub>s</sub> 被低估；深處射線接近垂直，偏差消失。');
    function draw() {
      const { v1, h1, v2, x } = ctl;
      const zs = linspace(1, 20, 20);
      const tv = z => 1000 * (Math.min(z, h1) / v1 + Math.max(z - h1, 0) / v2);
      const tr = zs.map(z => tv(z) * Math.sqrt(z * z + x * x) / z);
      const tc = zs.map((z, i) => tr[i] * z / Math.sqrt(z * z + x * x));
      const use = ctl.cor ? tc : tr;
      const ymax = 1.1 * Math.max(...tr, tv(20));
      pt.setLim(null, [0, ymax]); pt.clear(); pv.clear();
      pt.fn(tv, { c: 'gray', w: 1.4, range: [0, 20] });
      pt.line([0, ...zs], [0, ...tr], { c: 'ink', dash: '5 4', w: 1 });
      pt.markers(zs, tr, { c: 'ink', r: 3.4 });
      if (ctl.cor) pt.markers(zs, tc, { c: 'blue', fill: true, r: 3, over: true });
      pt.vline(h1);
      const iv = arr => zs.map((z, i) => 1000 * (i ? 1 : z) / (arr[i] - (i ? arr[i - 1] : 0)));
      const mid = zs.map(z => z - 0.5);
      const vr = iv(tr), vc = iv(tc);
      pv.setLim(null, [0, 1.25 * Math.max(v1, v2, ...vr.filter(isFinite))]);
      pv.clear();
      pv.line([0, h1, h1, 20], [v1, v1, v2, v2], { c: 'gray', w: 2 });
      pv.line(mid, vr, { c: 'ink', dash: '5 4', w: 1 }); pv.markers(mid, vr, { c: 'ink', r: 3.4 });
      if (ctl.cor) pv.markers(mid, vc, { c: 'blue', fill: true, r: 3, over: true });
      pv.vline(h1);
      const sel = zs.map((z, i) => i).filter(i => zs[i] <= h1);
      const regOf = arr => 1000 / linreg([0, ...sel.map(i => zs[i])], [0, ...sel.map(i => arr[i])]).slope;
      st({ t12: fx(tv(12), 1), raw1: fx(vr[0], 0) + ' m/s', cor1: fx(vc[0], 0) + ' m/s', reg: fx(regOf(tr), 0) + ' m/s', regc: fx(regOf(tc), 0) + ' m/s' }, { raw1: Math.abs(vr[0] / v1 - 1) > 0.1, reg: Math.abs(regOf(tr) / v1 - 1) > 0.05 });
      void use;
    }
    draw();
  })();

  /* =====================================================================
     W11 試驗應變範圍、共振柱頻率方程式、反例
     ===================================================================== */
  (() => {
    const body = SD.widget('w-lab'); if (!body) return;
    const panes = SD.tabs(body, ['試驗的應變範圍', '共振柱：β tan β = I/I₀', '反例：講義 p. 52–53']);

    // --- 應變範圍
    (() => {
      const pane = panes[0];
      const TESTS = [
        ['現地震測（下孔、跨孔、MASW）', 1e-5, 1e-3, 'Gmax'],
        ['彎曲元件', 1e-5, 1e-3, 'Gmax'],
        ['共振柱', 1e-4, 1e-1, '小到中應變的 G/Gmax、D'],
        ['反覆扭剪', 1e-3, 1, '中應變到約 1 % 的迴圈'],
        ['反覆三軸', 1e-2, 5, '大應變 G、D、孔隙水壓'],
        ['反覆單剪', 1e-2, 5, '大應變；應力路徑最接近現地']
      ];
      const ctl = SD.controls(pane, [{ k: 'g', label: '分析需要的應變 γ (%)', min: 1e-5, max: 5, log: true, snap: 2, value: 0.05, fmt: v => pctFmt(v) + ' %' }], draw);
      const p = new Plot(pane, { w: 680, h: 300, logx: true, xlim: [1e-5, 10], ylim: [0, TESTS.length + 0.6], yticks: [], m: { l: 168, r: 16, t: 10, b: 48 }, xlabel: '剪應變 γ (%)' });
      const st = SD.stats(pane, [{ k: 'n', label: '涵蓋此應變的試驗', hl: true }, { k: 'what', label: '它們提供' }]);
      SD.div(pane, 'w-note', '範圍為約略值（講義 3.7 的示意圖），隨儀器與土壤而變；反覆扭剪的下限為估計。虛線為砂的 γ<sub>tl</sub>、γ<sub>tv</sub>。');
      function draw() {
        p.clear();
        p.vline(1e-3, { c: 'ink', dash: '4 3' }); p.vline(1e-2, { c: 'ink', dash: '4 3' });
        const hit = [];
        TESTS.forEach(([n, a, b, w], i) => {
          const y = TESTS.length - i;
          const on = ctl.g >= a && ctl.g <= b;
          if (on) hit.push([n, w]);
          const X0 = p.sx(a), X1 = p.sx(b), Y = p.sy(y);
          svgEl('rect', { x: X0, y: Y - 9, width: X1 - X0, height: 18, rx: 4, class: on ? 'f-blue' : 'f-gray', opacity: on ? 0.85 : 0.35 }, p.bandG);
          const t = svgEl('text', { x: p.x0 - 8, y: Y + 4.5, 'text-anchor': 'end', class: on ? 'txt-blue' : 'txt' }, p.overG);
          t.textContent = n;
        });
        p.vline(ctl.g, { c: 'red', w: 2, dash: null });
        st({ n: hit.length ? hit.map(h => h[0]).join('、') : '沒有', what: hit.length ? hit.map(h => h[1]).join('；') : '—' });
      }
      draw();
    })();

    // --- 共振柱頻率方程式
    (() => {
      const pane = panes[1];
      const ctl = SD.controls(pane, [
        { k: 'r', label: 'I/I_0（試體 / 頂部質量）', min: 0.01, max: 100, log: true, snap: 3, value: 0.1, fmt: v => v.toPrecision(3) },
        { k: 'H', label: '試體高 H (m)', type: 'number', value: 0.1 },
        { k: 'fn', label: '共振頻率 f_n (Hz)', type: 'number', value: 108.22 },
        { k: 'rho', label: 'ρ (t/m³)', type: 'number', value: 1.95 }
      ], draw);
      SD.buttons(pane, '範例', [
        { label: '例：共振柱 stage 1', on: () => ctl.set({ r: 0.1, H: 0.1, fn: 108.22, rho: 1.95 }) },
        { label: '講義備註例：β = 0.5、100 Hz', on: () => ctl.set({ r: +(0.5 * Math.tan(0.5)).toFixed(4), H: 0.1, fn: 100, rho: 1.9 }) }
      ]);
      const p = new Plot(pane, { w: 640, h: 300, logx: true, xlim: [0.01, 100], ylim: [0, 1.7], xlabel: 'I/I_0', ylabel: 'β = ω_nH/V_s' });
      p.legend([{ label: 'β tan β = I/I_0 的第一個根', c: 'blue', w: 2 }, { label: 'π/2（頂部無質量）', c: 'gray', dash: true }, { label: '√(I/I_0)（頂部很重的 SDOF 近似）', c: 'red', dash: true }]);
      const st = SD.stats(pane, [{ k: 'b', label: 'β (rad)', hl: true }, { k: 'ba', label: '√(I/I_0) 近似' }, { k: 'vs', label: 'V_s = 2πf_nH/β (m/s)', hl: true }, { k: 'G', label: 'G = ρV_s² (MPa)', hl: true }]);
      const root = r => bisect(b => b * Math.tan(b) - r, 1e-9, PI / 2 - 1e-9);
      const rs = logspace(0.01, 100, 200), bs = rs.map(root);
      function draw() {
        p.clear();
        p.hline(PI / 2, { dash: '6 4' });
        p.line(rs, rs.map(Math.sqrt), { c: 'red', dash: '6 4', w: 1.3 });
        p.line(rs, bs, { c: 'blue', w: 2 });
        const b = root(ctl.r), vs = 2 * PI * ctl.fn * ctl.H / b;
        p.markers([ctl.r], [b], { c: 'blue', fill: true, over: true });
        st({ b: fx(b, 4), ba: fx(Math.sqrt(ctl.r), 4), vs: fx(vs, 1), G: fx(ctl.rho * vs * vs / 1000, 2) });
      }
      draw();
    })();

    // --- 反例（合成資料）
    (() => {
      const pane = panes[2];
      const ctl = SD.controls(pane, [
        { k: 'mode', label: '畫法', type: 'select', value: 'good', options: [['good', '正確：對數應變軸、G_max 取小應變平台'], ['bad', '錯誤：線性應變軸、G_max 取離群值']] },
        { k: 'mpa', label: '把 G 軸誤標成 MPa', type: 'check', value: false }
      ], draw);
      const rr = SD.rng(52);
      const gs = logspace(1e-5, 6e-3, 14);
      const Gm = 120000, Gd = gs.map(g => Gm / (1 + Math.pow(g / 0.05, 0.919)) * (1 + 0.025 * rr.normal()));
      Gd[9] = 470000; // 一筆離群值（例如輸入錯誤）
      const row = SD.div(pane, 'plot-row');
      const pG = new Plot(row, { w: 420, h: 300, logx: true, xlim: [1e-5, 1e-2], ylim: [0, 500000], xlabel: '剪應變 γ (%)', ylabel: 'G (kPa)', yfmt: v => (v / 1000).toFixed(0) + 'k', m: { l: 58 } });
      const pR = new Plot(row, { w: 420, h: 300, logx: true, xlim: [1e-5, 1e-2], ylim: [0, 1.2], xlabel: '剪應變 γ (%)', ylabel: 'G/G_{max}' });
      pR.legend([{ label: '合成試驗點', c: 'ink', dash: true, marker: 'circle' }]);
      const st = SD.stats(pane, [{ k: 'gmax', label: '採用的 G_{max}' }, { k: 'vs', label: '對應 V_s（ρ = 1.9）' }, { k: 'chk', label: '合理性檢查', hl: true }]);
      SD.div(pane, 'w-note', '合成資料，不是講義原始數據：一個 G<sub>max</sub> ≈ 120 MPa 的砂，在 0.0006 % 附近混入一筆 470 000 kPa 的離群值。用離群值當 G<sub>max</sub>，整條 G/G<sub>max</sub> 都掉到 0.3 以下；線性應變軸又把小應變的點全擠在左邊，看不出曲線形狀。');
      function draw() {
        const bad = ctl.mode === 'bad';
        [pG, pR].forEach(p => { p.o.logx = !bad; });
        pG.o.ylabel = ctl.mpa ? 'G (MPa)　← 標錯' : 'G (kPa)';
        pG.o.yfmt = ctl.mpa ? (v => v.toFixed(0)) : (v => (v / 1000).toFixed(0) + 'k');
        pG.o.m.l = ctl.mpa ? 70 : 58;
        const xl = bad ? [0, 0.007] : [1e-5, 1e-2];
        pG.setLim(xl, null); pR.setLim(xl, null);
        pG.clear(); pR.clear();
        pG.line(gs, Gd, { c: 'ink', dash: '5 4', w: 1 }); pG.markers(gs, Gd, { c: 'ink', r: 3.6 });
        const gmax = bad ? Math.max(...Gd) : (Gd[0] + Gd[1] + Gd[2]) / 3;
        const R = Gd.map(v => v / gmax);
        pR.hline(1, { dash: '4 3' });
        pR.line(gs, R, { c: 'ink', dash: '5 4', w: 1 }); pR.markers(gs, R, { c: 'ink', r: 3.6 });
        const vs = Math.sqrt(gmax / 1.9);
        st({ gmax: ctl.mpa ? f(gmax, 0) + ' MPa？' : f(gmax / 1000, 1) + ' MPa', vs: ctl.mpa ? f(Math.sqrt(gmax * 1000 / 1.9), 0) + ' m/s' : f(vs, 0) + ' m/s',
          chk: ctl.mpa ? 'V<sub>s</sub> 上萬 m/s，比岩石還硬：單位應是 kPa' : bad ? '小應變端 G/G<sub>max</sub> ≈ 0.25，不可能' : '小應變端 ≈ 1，合理' }, { chk: ctl.mpa || bad });
      }
      draw();
    })();
  })();

  /* =====================================================================
     W12 現地 Gmax × 實驗室曲線
     ===================================================================== */
  (() => {
    const body = SD.widget('w-eqlin'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'src', label: '曲線來源', type: 'select', value: 'rc', options: [['rc', '共振柱實測資料（log γ 內插）'], ['dar', 'Darendeli (2001)'], ['man', '手動輸入 G/Gmax 與 D']] },
      { k: 'vs', label: '現地 V_s (m/s)', min: 80, max: 500, step: 1, value: 269, fmt: v => v.toFixed(0) },
      { k: 'rho', label: 'ρ (t/m³)', min: 1.5, max: 2.3, step: 0.01, value: 1.95, fmt: v => v.toFixed(2) },
      { k: 'g', label: '有效應變 γ (%)', min: 5e-4, max: 0.1, log: true, snap: 2, value: 0.05, fmt: v => pctFmt(v) + ' %' },
      { k: 'PI', label: 'PI（Darendeli）', min: 0, max: 100, step: 1, value: 0, fmt: v => v.toFixed(0) },
      { k: 'sm', label: 'σ′m (atm)（Darendeli）', min: 0.25, max: 16, log: true, snap: 2, value: 1, fmt: v => v.toPrecision(2) },
      { k: 'mg', label: '手動 G/Gmax', type: 'number', value: 0.6 },
      { k: 'mD', label: '手動 D (%)', type: 'number', value: 8 }
    ], draw);
    SD.buttons(body, '範例', [
      { label: '講義例：200 m/s、1.9、G/Gmax = 0.6', on: () => ctl.set({ src: 'man', vs: 200, rho: 1.9, g: 0.05, mg: 0.6 }) },
      { label: '例：Vs = 269 m/s、0.05 %', on: () => ctl.set({ src: 'rc', vs: 269, rho: 1.95, g: 0.05 }) }
    ]);
    const row = SD.div(body, 'plot-row');
    const pG = new Plot(row, { w: 420, h: 290, logx: true, xlim: [1e-4, 1], ylim: [0, 160], xlabel: '剪應變 γ (%)', ylabel: 'G (MPa)' });
    const pD = new Plot(row, { w: 420, h: 290, logx: true, xlim: [1e-4, 1], ylim: [0, 25], xlabel: '剪應變 γ (%)', ylabel: 'D (%)' });
    pG.legend([{ label: '現地 G_{max} × 曲線', c: 'blue', w: 2 }, { label: '實驗室 G_{max} × 曲線（僅共振柱實測資料）', c: 'gray', dash: true }, { label: '共振柱資料點', c: 'ink', dash: true, marker: 'circle' }]);
    const st = SD.stats(body, [
      { k: 'gm', label: '現地 G_{max} = ρV_s² (MPa)' }, { k: 'gg', label: 'G/G_{max}', hl: true }, { k: 'G', label: '現地 G (MPa)', hl: true },
      { k: 'D', label: 'D (%)', hl: true }, { k: 'vse', label: '等效 V_s (m/s)' }, { k: 'Gl', label: '若只用實驗室 G_{max} (MPa)' }
    ]);
    const gs = logspace(1e-4, 1, 200);
    function curve(g) {
      if (ctl.src === 'rc' && HW) { const r = HW.fit1.rows; return { gg: interpLog(g, r.map(q => q.g), r.map(q => q.gg)), D: interpLog(g, r.map(q => q.g), r.map(q => q.D)) }; }
      if (ctl.src === 'man') return { gg: ctl.mg, D: ctl.mD };
      const r = DR(g, { PI: ctl.PI, sm: ctl.sm }); return { gg: r.gg, D: r.D };
    }
    function draw() {
      const gm = ctl.rho * ctl.vs * ctl.vs / 1000, c = curve(ctl.g), lab = HW ? HW.fit1.Gmax / 1000 : NaN;
      pG.setLim(null, [0, 1.15 * Math.max(gm, ctl.src === 'rc' ? lab : 0)]);
      pG.clear(); pD.clear();
      pG.vline(ctl.g); pD.vline(ctl.g);
      if (ctl.src !== 'man') {
        const cs = gs.map(curve);
        pG.line(gs, cs.map(q => q.gg * gm), { c: 'blue', w: 2 });
        pD.line(gs, cs.map(q => q.D), { c: 'blue', w: 2 });
        if (ctl.src === 'rc' && HW) {
          const r = HW.fit1.rows;
          pG.line(gs, cs.map(q => q.gg * lab), { c: 'gray', dash: '6 4', w: 1.4 });
          pG.markers(r.map(q => q.g), r.map(q => q.gg * gm), { c: 'ink', r: 3.6 });
          pD.line(r.map(q => q.g), r.map(q => q.D), { c: 'ink', dash: '5 4', w: 1 });
          pD.markers(r.map(q => q.g), r.map(q => q.D), { c: 'ink', r: 3.6 });
        }
      }
      pG.markers([ctl.g], [c.gg * gm], { c: 'blue', fill: true, over: true, r: 4.5 });
      pD.markers([ctl.g], [c.D], { c: 'blue', fill: true, over: true, r: 4.5 });
      st({ gm: fx(gm, 1), gg: fx(c.gg, 3), G: fx(c.gg * gm, 1), D: fx(c.D, 2), vse: fx(ctl.vs * Math.sqrt(c.gg), 0), Gl: ctl.src === 'rc' ? fx(c.gg * lab, 1) : '—' });
    }
    draw();
  })();

  /* =====================================================================
     W13 HW3 Q1
     ===================================================================== */
  (() => {
    const body = SD.widget('w-hw3q1'); if (!body || !HW) return;
    const panes = SD.tabs(body, ['Q1(a) stage 1 掃頻', 'Q1(b) G/Gmax 與 D', '六個階段的自由衰減']);
    const sp = HW.d.spec, q = HW.q1a;

    // (a)
    (() => {
      const pane = panes[0], sw = HW.d.sweep;
      SD.div(pane, 'w-note', '資料：resonant_column_sweep.csv（96–120 Hz，97 點，加速度計 mV）。以扭矩激振 SDOF 的加速度振幅 A = A<sub>0</sub>r²/√[(1−r²)² + (2ξr)²] 做最小平方擬合（Nelder–Mead，A<sub>0</sub> 線性解出），結果在網頁中即時計算。');
      const p = new Plot(pane, { w: 640, h: 320, xlim: [96, 120], ylim: [0, 66], xlabel: '頻率 f (Hz)', ylabel: '加速度計輸出 (mV)' });
      p.legend([{ label: '量測', c: 'ink', dash: true, marker: 'circle' }, { label: 'SDOF 擬合', c: 'blue', w: 2 }, { label: '峰值 / √2', c: 'gray', dash: true }, { label: '半功率點', c: 'red', line: false, marker: 'circle', fill: true }]);
      p.line(sw.f, sw.a, { c: 'ink', dash: '5 4', w: 1 });
      p.markers(sw.f, sw.a, { c: 'ink', r: 2.8 });
      const shape = ff => { const r = ff / q.fn; return q.A0 * r * r / Math.sqrt((1 - r * r) ** 2 + (2 * q.xi * r) ** 2); };
      p.fn(shape, { c: 'blue', w: 2, n: 800 });
      p.hline(q.lvl, { dash: '5 4' });
      p.markers([q.f1, q.f2], [q.lvl, q.lvl], { c: 'red', fill: true, over: true });
      p.readout(ff => `f = ${ff.toFixed(2)} Hz　擬合 ${shape(ff).toFixed(2)} mV`);
      const st = SD.stats(pane, [
        { k: 'pk', label: '量測最大點 (Hz)' }, { k: 'par', label: '三點拋物線峰 (Hz)' }, { k: 'fn', label: '擬合 f_n (Hz)', hl: true }, { k: 'xi', label: '擬合 ξ (%)' },
        { k: 'hp', label: '半功率 D = (f_2−f_1)/2f_n (%)' }, { k: 'b', label: 'β（β tan β = 0.10）' }, { k: 'vs', label: 'V_s = 2πf_nH/β (m/s)', hl: true }, { k: 'G', label: 'G_{max} = ρV_s² (MPa)', hl: true }
      ]);
      const r1 = HW.fit1.rows[0];
      st({ pk: fx(q.fpk, 2), par: fx(q.fpar, 2), fn: fx(q.fn, 2), xi: fx(q.xi * 100, 2), hp: fx(q.Dhp * 100, 2), b: fx(HW.beta, 4), vs: fx(r1.vs, 1), G: fx(r1.G / 1000, 1) });
      SD.div(pane, 'w-note', `stage 1 對數遞減：δ = (1/5) ln(${HW.d.stages[0].X[0]}/${HW.d.stages[0].X[5]}) = ${fx(r1.del, 4)}，D = δ/2π = ${fx(r1.D, 2)} %，與半功率頻寬 ${fx(q.Dhp * 100, 2)} % 一致。`);
    })();

    // (b)
    (() => {
      const pane = panes[1];
      const ctl = SD.controls(pane, [
        { k: 'src', label: 'stage 1 的 f_n', type: 'select', value: 'fit', options: [['fit', '掃頻擬合 ' + q.fn.toFixed(2) + ' Hz（學生採用）'], ['tab', 'Table 1 的 108.3 Hz']] },
        { k: 'norm', label: 'Darendeli 也除以它在 stage 1 應變（0.0005 %）的值', type: 'check', value: false },
        { k: 'hyp', label: '畫雙曲線擬合（γr 由六點最小平方）', type: 'check', value: true }
      ], draw);
      const row = SD.div(pane, 'plot-row');
      const pG = new Plot(row, { w: 420, h: 300, logx: true, xlim: [1e-4, 1], ylim: [0, 1.1], xlabel: '剪應變 γ (%)', ylabel: 'G/G_{max}' });
      const pD = new Plot(row, { w: 420, h: 300, logx: true, xlim: [1e-4, 1], ylim: [0, 22], xlabel: '剪應變 γ (%)', ylabel: '阻尼比 D (%)' });
      pG.legend([{ label: '共振柱（HW3）', c: 'ink', dash: true, marker: 'circle' }, { label: 'Darendeli PI = 0、1 atm', c: 'blue', w: 2 }, { label: '雙曲線擬合', c: 'red', dash: true }, { label: 'γ_{tl}、γ_{tv}', c: 'gray', dash: '4 3' }]);
      const tbl = SD.div(pane, 'tablewrap');
      const st = SD.stats(pane, [
        { k: 'gmax', label: 'G_{max}（stage 1）(MPa)' }, { k: 'gi', label: 'γ_r：G/G_{max} = 0.5 內插 (%)', hl: true }, { k: 'gh', label: 'γ_r：雙曲線擬合 (%)', hl: true },
        { k: 'gm', label: 'γ_r：修正雙曲線 a = 0.919 (%)' }, { k: 'gd', label: 'γ_r：Darendeli (%)' }, { k: 'dd', label: '|G/G_{max} − Darendeli| 範圍' }
      ]);
      const gs = logspace(1e-4, 1, 200);
      function draw() {
        const S = ctl.src === 'fit' ? HW.fit1 : HW.tab1, rows = S.rows;
        const n0 = ctl.norm ? DR(rows[0].g).gg : 1;
        [pG, pD].forEach(p => { p.clear(); p.vline(1e-3, { dash: '4 3' }); p.vline(1e-2, { dash: '4 3' }); });
        pG.text(1e-3, 0.08, 'γ_{tl}', { dx: -22, cls: 'txt-muted' }); pG.text(1e-2, 0.08, 'γ_{tv}', { dx: 4, cls: 'txt-muted' });
        pG.line(gs, gs.map(g => DR(g).gg / n0), { c: 'blue', w: 2 });
        pD.line(gs, gs.map(g => DR(g).D), { c: 'blue', w: 2 });
        if (ctl.hyp) pG.line(gs, gs.map(g => 1 / (1 + g / S.grH)), { c: 'red', dash: '6 4', w: 1.4 });
        const g = rows.map(r => r.g);
        pG.line(g, rows.map(r => r.gg), { c: 'ink', dash: '5 4', w: 1 }); pG.markers(g, rows.map(r => r.gg), { c: 'ink', r: 4 });
        pD.line(g, rows.map(r => r.D), { c: 'ink', dash: '5 4', w: 1 }); pD.markers(g, rows.map(r => r.D), { c: 'ink', r: 4 });
        pG.hline(0.5, { dash: '2 3' });
        pG.markers([S.grInt], [0.5], { c: 'red', fill: true, over: true, r: 3.6 });
        const zone = gg => gg < 0.999e-3 ? '&lt; γ<sub>tl</sub>' : gg <= 1.001e-3 ? '≈ γ<sub>tl</sub>' : gg < 0.999e-2 ? 'γ<sub>tl</sub>–γ<sub>tv</sub>' : gg <= 1.001e-2 ? '≈ γ<sub>tv</sub>' : '&gt; γ<sub>tv</sub>';
        tbl.innerHTML = `<table class="data"><tr><th>Stage</th><th class="num">γ (%)</th><th class="num">f<sub>n</sub> (Hz)</th><th class="num">V<sub>s</sub> (m/s)</th><th class="num">G (MPa)</th><th class="num">G/G<sub>max</sub></th><th class="num">δ</th><th class="num">D (%)</th><th class="num">Darendeli G/G<sub>max</sub></th><th class="num">Darendeli D (%)</th><th>範圍</th></tr>` +
          rows.map((r, i) => `<tr${i === 5 ? ' class="hl"' : ''}><td>${i + 1}</td><td class="num">${r.g}</td><td class="num">${r.fn.toFixed(2)}</td><td class="num">${r.vs.toFixed(1)}</td><td class="num">${(r.G / 1000).toFixed(1)}</td><td class="num">${r.gg.toFixed(3)}</td><td class="num">${r.del.toFixed(4)}</td><td class="num">${r.D.toFixed(2)}</td><td class="num">${(r.drg / n0).toFixed(3)}</td><td class="num">${r.drD.toFixed(2)}</td><td>${zone(r.g)}</td></tr>`).join('') + '</table>';
        const diffs = rows.map(r => r.gg - r.drg / n0);
        st({ gmax: fx(S.Gmax / 1000, 2), gi: fx(S.grInt, 4), gh: fx(S.grH, 4), gm: fx(S.grM, 4), gd: fx(DR(1).gr, 4), dd: fx(Math.min(...diffs), 3) + ' ～ ' + fx(Math.max(...diffs), 3) });
      }
      draw();
    })();

    // (c) 自由衰減
    (() => {
      const pane = panes[2];
      const p = new Plot(pane, { w: 640, h: 320, xlim: [0.5, 6.5], ylim: [-4.5, 0.3], xticks: [1, 2, 3, 4, 5, 6], xlabel: '峰的序號 i', ylabel: 'ln(X_i/X_1)' });
      const shapes = [['circle', false], ['tri', false], ['square', false], ['circle', true], ['tri', true], ['square', true]];
      p.legend(HW.fit1.rows.map((r, i) => ({ label: 'stage ' + (i + 1) + '（' + r.g + ' %）', c: 'ink', line: false, marker: shapes[i][0], fill: shapes[i][1] })).concat([{ label: '迴歸直線', c: 'blue' }]));
      HW.fit1.rows.forEach((r, i) => {
        const xs = [1, 2, 3, 4, 5, 6], ys = r.X.map(v => Math.log(v / r.X[0]));
        const b0 = r.icpt - Math.log(r.X[0]);
        p.line([1, 6], [b0 - r.delr, b0 - 6 * r.delr], { c: 'blue', w: 1.2 });
        p.markers(xs, ys, { c: 'ink', shape: shapes[i][0], fill: shapes[i][1], r: 3.8, over: true });
      });
      const t = SD.div(pane, 'tablewrap');
      t.innerHTML = `<table class="data"><tr><th>Stage</th><th class="num">X<sub>1</sub>/X<sub>6</sub></th><th class="num">δ（兩點）</th><th class="num">δ（迴歸）</th><th class="num">D = δ/2π (%)</th><th class="num">精確式 (%)</th></tr>` +
        HW.fit1.rows.map((r, i) => `<tr><td>${i + 1}</td><td class="num">${(r.X[0] / r.X[5]).toFixed(3)}</td><td class="num">${r.del.toFixed(4)}</td><td class="num">${r.delr.toFixed(4)}</td><td class="num">${r.D.toFixed(2)}</td><td class="num">${r.Dx.toFixed(2)}</td></tr>`).join('') + '</table>';
      SD.div(pane, 'w-note', '斜率就是 −δ。小應變的 stage 1–3 幾乎是直線（等比例衰減，線性黏彈）；stage 6 斜率最陡，D = 13 %。兩點公式與迴歸差不到 2 %。自由衰減時振幅一路變小，量到的是一段應變範圍的平均阻尼。');
    })();
  })();

  /* =====================================================================
     W14 HW3 Q2 Masing
     ===================================================================== */
  (() => {
    const body = SD.widget('w-hw3q2'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'gc', label: 'γ_c (%)', min: 0.001, max: 1, log: true, snap: 3, value: 0.1, fmt: v => pctFmt(v) + ' %' },
      { k: 'G', label: 'G_max (kPa)', type: 'number', value: 94000 },
      { k: 'gref', label: 'γ_ref（小數，與網頁工具相同）', type: 'number', value: 0.0004 }
    ], draw);
    SD.buttons(body, '題目', [{ label: 'γc = 0.01 %', on: () => ctl.set({ gc: 0.01 }) }, { label: 'γc = 0.1 %', on: () => ctl.set({ gc: 0.1 }) }]);
    const row = SD.div(body, 'plot-row');
    const p = new Plot(row, { w: 420, h: 320, xlim: [-1, 1], ylim: [-1, 1], xlabel: '剪應變 γ (%)', ylabel: '剪應力 τ (kPa)' });
    const pd = new Plot(row, { w: 420, h: 320, logx: true, xlim: [1e-4, 1], ylim: [0, 50], xlabel: '剪應變 γ (%)', ylabel: '阻尼比 D (%)' });
    p.legend([{ label: '骨架', c: 'gray', dash: true }, { label: 'Masing 迴圈', c: 'blue', w: 2 }, { label: 'E_D（灰）、E_{So}（淡藍）', c: 'gray', line: false, marker: 'square', fill: true }]);
    pd.legend([{ label: 'Masing（此骨架）', c: 'red', w: 1.6 }, { label: 'Darendeli', c: 'blue', w: 2 }, { label: '共振柱量測', c: 'ink', dash: true, marker: 'circle' }]);
    const st = SD.stats(body, [
      { k: 'tc', label: 'τ_c (kPa)' }, { k: 'gs', label: 'G_{sec}/G_{max}' }, { k: 'ed', label: 'E_D (J/m³)' }, { k: 'es', label: 'E_{So} (J/m³)' },
      { k: 'D', label: 'D Masing（梯形積分）', hl: true }, { k: 'Dc', label: 'D 封閉解' }, { k: 'Dm', label: 'D 量測（Q1，log 內插）', hl: true }, { k: 'Dd', label: 'D Darendeli' }
    ]);
    const note = SD.div(body, 'w-note');
    const gs = logspace(1e-4, 1, 200);
    function draw() {
      const gc = ctl.gc, Gm = ctl.G, grP = ctl.gref * 100, fb = hyp(Gm, grP);
      const path = masingSimple([0, gc, -gc, gc], fb, 200, true);
      const seg = path.g.map((g, i) => [g, path.t[i], path.seg[i]]);
      const s23 = seg.filter(s => s[2] >= 2);
      let ed = 0; for (let i = 1; i < s23.length; i++) ed += (s23[i][0] - s23[i - 1][0]) / 100 * (s23[i][1] + s23[i - 1][1]) / 2;
      ed = Math.abs(ed);
      const A = seg.filter(s => s[2] === 1).pop(), tc = A[1], eso = 0.5 * tc * gc / 100, D = ed / (4 * PI * eso);
      const X = 1.3 * gc, Y = 1.3 * fb(X);
      p.setLim([-X, X], [-Y, Y]); p.clear(); pd.clear();
      p.hline(0, { dash: null, w: 0.6 }); p.vline(0, { dash: null, w: 0.6 });
      fillPoly(p, s23.map(s => s[0]), s23.map(s => s[1]), 'f-band');
      fillPoly(p, [0, gc, gc], [0, 0, tc], 'f-blue', 0.25);
      p.fn(fb, { c: 'gray', dash: '6 4', w: 1.1, range: [-X, X] });
      p.line(path.g, path.t, { c: 'blue', w: 2 });
      p.markers([gc, -gc], [tc, -tc], { c: 'ink', r: 3.8, over: true });
      p.text(gc, tc, 'A', { dx: 7, dy: -4 }); p.text(-gc, -tc, 'B', { dx: -16, dy: 12 });
      pd.line(gs, gs.map(g => 100 * masingA1(g / grP)), { c: 'red', w: 1.6 });
      pd.line(gs, gs.map(g => DR(g).D), { c: 'blue', w: 2 });
      let Dm = NaN;
      if (HW) {
        const r = HW.fit1.rows;
        pd.line(r.map(q => q.g), r.map(q => q.D), { c: 'ink', dash: '5 4', w: 1 }); pd.markers(r.map(q => q.g), r.map(q => q.D), { c: 'ink', r: 3.8 });
        if (gc >= r[0].g && gc <= r[5].g) Dm = interpLog(gc, r.map(q => q.g), r.map(q => q.D));
      }
      pd.vline(gc, { c: 'blue', dash: '4 3' });
      pd.markers([gc], [100 * D], { c: 'red', fill: true, over: true, r: 4.5 });
      const Dd = DR(gc).D;
      st({ tc: fx(tc, 2), gs: fx(tc / (gc / 100) / Gm, 3), ed: f(ed * 1000, 4), es: f(eso * 1000, 4), D: fx(100 * D, 2) + ' %', Dc: fx(100 * masingA1(gc / grP), 2) + ' %', Dm: isFinite(Dm) ? fx(Dm, 2) + ' %' : '超出試驗範圍', Dd: fx(Dd, 2) + ' %' });
      note.innerHTML = isFinite(Dm) ? `Masing／量測 = <b>${fx(100 * D / Dm, 2)}</b>。γ<sub>c</sub>/γ<sub>ref</sub> = ${fx(gc / grP, 2)}：比值越大，迴圈越胖，Masing 高估越多；γ<sub>c</sub> 很小時 Masing 阻尼趨近 0，反而低於量測（試著拉到 0.001 %）。每段 200 點，與網頁工具相同。` : 'γ<sub>c</sub> 超出共振柱試驗範圍（0.0005–0.1 %），沒有量測值可比較。';
    }
    draw();
  })();

  /* =====================================================================
     W15 HW3 Q3 下孔試驗
     ===================================================================== */
  (() => {
    const body = SD.widget('w-hw3q3'); if (!body || !HW) return;
    const dh = HW.d.downhole, Z = HW.z, T = HW.t;
    const ctl = SD.controls(body, [
      { k: 'm', label: '各層 V_s 的求法', type: 'select', value: 'reg', options: [['reg', '層頂、層內、層底受波器迴歸（學生採用）'], ['end', '只用層頂與層底兩點']] },
      { k: 'pa', label: 'Hardin 的 p_a', type: 'select', value: '101.325', options: [['101.325', '101.325 kPa（學生採用）'], ['100', '100 kPa']] },
      { k: 'g', label: '第 3 層有效應變 γ (%)', min: 0.003, max: 0.1, log: true, snap: 2, value: 0.05, fmt: v => pctFmt(v) + ' %' }
    ], draw);
    const row = SD.div(body, 'plot-row');
    const pt = new Plot(row, { w: 420, h: 300, xlim: [0, 20], ylim: [0, 100], xlabel: '深度 z (m)', ylabel: '初達時間 t (ms)' });
    const pv = new Plot(row, { w: 420, h: 300, xlim: [0, 20], ylim: [0, 320], xlabel: '深度 z (m)', ylabel: 'V_s (m/s)' });
    pt.legend([{ label: '量測（含原點）', c: 'ink', dash: true, marker: 'circle' }, { label: '各層直線', c: 'blue', w: 2 }]);
    pv.legend([{ label: '各層 V_s', c: 'blue', w: 2 }, { label: '區間速度', c: 'ink', dash: true, marker: 'circle' }]);
    const t1 = SD.div(body, 'tablewrap'), t2 = SD.div(body, 'tablewrap'), t3 = SD.div(body, 'tablewrap');
    const note = SD.div(body, 'w-note');
    function draw() {
      pt.clear(); pv.clear();
      pt.line(Z, T, { c: 'ink', dash: '5 4', w: 1 }); pt.markers(Z, T, { c: 'ink', r: 4 });
      const lay = dh.layers.map(([a, b, name], i) => {
        const idx = Z.map((z, j) => j).filter(j => Z[j] >= a && Z[j] <= b);
        let slope, icpt;
        if (ctl.m === 'reg') ({ slope, icpt } = linreg(idx.map(j => Z[j]), idx.map(j => T[j])));
        else { const j0 = idx[0], j1 = idx[idx.length - 1]; slope = (T[j1] - T[j0]) / (Z[j1] - Z[j0]); icpt = T[j0] - slope * Z[j0]; }
        const vs = 1000 / slope, zm = (a + b) / 2, s = HW.stress(zm);
        const res = Math.max(...idx.map(j => Math.abs(T[j] - (icpt + slope * Z[j]))));
        pt.line([a, b], [icpt + slope * a, icpt + slope * b], { c: 'blue', w: 2 });
        pv.line([a, b], [vs, vs], { c: 'blue', w: 2.2 });
        if (i) { pt.vline(a); pv.vline(a); }
        return { name, a, b, slope, icpt, vs, G: dh.rho * vs * vs, sm: s.sm, res };
      });
      const iz = [], iv = [];
      for (let j = 1; j < Z.length; j++) { iz.push((Z[j] + Z[j - 1]) / 2); iv.push(1000 * (Z[j] - Z[j - 1]) / (T[j] - T[j - 1])); }
      pv.line(iz, iv, { c: 'ink', dash: '5 4', w: 1 }); pv.markers(iz, iv, { c: 'ink', r: 4 });
      t1.innerHTML = '<table class="data"><tr><th>土層</th><th class="num">斜率 (ms/m)</th><th class="num">截距 (ms)</th><th class="num">最大殘差 (ms)</th><th class="num">V<sub>s</sub> (m/s)</th><th class="num">G<sub>max</sub> (MPa)</th><th class="num">層中點 σ′<sub>m</sub> (kPa)</th></tr>' +
        lay.map((l, i) => `<tr${i === 2 ? ' class="hl"' : ''}><td>${i + 1} ${l.name}（${l.a}–${l.b} m）</td><td class="num">${l.slope.toFixed(3)}</td><td class="num">${l.icpt.toFixed(2)}</td><td class="num">${l.res.toFixed(2)}</td><td class="num">${l.vs.toFixed(1)}</td><td class="num">${(l.G / 1000).toFixed(1)}</td><td class="num">${l.sm.toFixed(1)}</td></tr>`).join('') + '</table>';
      // 第 3 層比較
      const s15 = HW.stress(dh.zSample), pa = +ctl.pa, e = HW.d.spec.e, F = 1 / (0.3 + 0.7 * e * e);
      const gH = 625 * F * Math.sqrt(pa * s15.sm), L3 = lay[2], rc = HW.fit1.rows[0];
      const rowsC = [['下孔試驗（現地）', L3.vs, L3.G], ['共振柱（重模）', rc.vs, rc.G], [`Hardin (1978)，σ′<sub>m</sub> = ${s15.sm.toFixed(1)} kPa`, Math.sqrt(gH / dh.rho), gH]];
      t2.innerHTML = '<table class="data"><tr><th>第 3 層 G<sub>max</sub></th><th class="num">V<sub>s</sub> (m/s)</th><th class="num">G<sub>max</sub> (MPa)</th><th class="num">與 Hardin 之比</th><th class="num">與共振柱之比</th></tr>' +
        rowsC.map((r, i) => `<tr${i === 0 ? ' class="hl"' : ''}><td>${r[0]}</td><td class="num">${r[1].toFixed(1)}</td><td class="num">${(r[2] / 1000).toFixed(1)}</td><td class="num">${(r[2] / gH).toFixed(2)}</td><td class="num">${(r[2] / rc.G).toFixed(2)}</td></tr>`).join('') + '</table>';
      // Q3(b)
      const g = ctl.g, rr = HW.fit1.rows, gsx = rr.map(q => q.g);
      const ggI = interpLog(g, gsx, rr.map(q => q.gg)), DI = interpLog(g, gsx, rr.map(q => q.D));
      const ggH = 1 / (1 + g / HW.fit1.grH), dr = DR(g);
      const rowsB = [['共振柱資料，log γ 內插（採用）', ggI, DI], [`雙曲線擬合（γ<sub>r</sub> = ${HW.fit1.grH.toFixed(4)} %）`, ggH, NaN], ['Darendeli (2001)，1 atm', dr.gg, dr.D]];
      t3.innerHTML = `<table class="data"><tr><th>γ = ${pctFmt(g)} % 的估計</th><th class="num">G/G<sub>max</sub></th><th class="num">D (%)</th><th class="num">G，現地 G<sub>max</sub> (MPa)</th><th class="num">G，實驗室 G<sub>max</sub> (MPa)</th><th class="num">等效 V<sub>s</sub> (m/s)</th></tr>` +
        rowsB.map((r, i) => `<tr${i === 0 ? ' class="hl"' : ''}><td>${r[0]}</td><td class="num">${r[1].toFixed(3)}</td><td class="num">${isFinite(r[2]) ? r[2].toFixed(2) : '—'}</td><td class="num">${(r[1] * L3.G / 1000).toFixed(1)}</td><td class="num">${(r[1] * rc.G / 1000).toFixed(1)}</td><td class="num">${Math.sqrt(r[1] * L3.G / dh.rho).toFixed(0)}</td></tr>`).join('') + '</table>';
      note.innerHTML = `15 m：σ<sub>v</sub> = ${s15.sv.toFixed(1)}、u = ${s15.u.toFixed(1)}、σ′<sub>v</sub> = ${s15.s1.toFixed(1)}、σ′<sub>m</sub> = ${s15.sm.toFixed(1)} kPa；F(0.70) = ${F.toFixed(4)}。改用「只用兩點」時第 3 層 V<sub>s</sub> 變成 ${(1000 * 10 / (T[6] - T[4])).toFixed(1)} m/s，差異小於 1 %，結論不變。γ = ${pctFmt(g)} % ${g > 0.01 ? '高於' : g < 0.001 ? '低於 γ<sub>tl</sub>，也低於' : '低於'} 砂的 γ<sub>tv</sub> ≈ 0.01 %${g > 0.01 ? '：飽和砂會累積孔隙水壓，等值線性的總應力分析算不出來。' : '。'}`;
    }
    draw();
  })();
});
