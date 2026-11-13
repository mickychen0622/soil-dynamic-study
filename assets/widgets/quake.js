/* =========================================================
   Ch0–2 導論、地震與強地動：互動元件
   資料：assets/data/ch02_*.js（HW1 手機紀錄、PEER 集集紀錄、B1 結果表）
   ========================================================= */

/* ---------- 訊號處理工具（與 HW1 的 A1/B1 程式同一套流程） ---------- */
const DSP = (() => {
  const TAU = 2 * Math.PI, G = 9.80665;
  /** RBJ 雙二階（= 雙線性轉換＋預扭曲的 Butterworth 節） */
  function biquad(type, f0, fs, Q) {
    const w0 = TAU * f0 / fs, c = Math.cos(w0), al = Math.sin(w0) / (2 * Q), a0 = 1 + al;
    const b = type === 'lp' ? [(1 - c) / 2, 1 - c, (1 - c) / 2] : [(1 + c) / 2, -(1 + c), (1 + c) / 2];
    return { b0: b[0] / a0, b1: b[1] / a0, b2: b[2] / a0, a1: -2 * c / a0, a2: (1 - al) / a0 };
  }
  /** 4 階 Butterworth = 兩個雙二階節（Q = 0.5412、1.3066） */
  const butter4 = (type, f0, fs) => [biquad(type, f0, fs, 0.5411961), biquad(type, f0, fs, 1.3065630)];
  function sosfilt(secs, x) {
    const y = Float64Array.from(x);
    secs.forEach(s => {
      let z1 = 0, z2 = 0;
      for (let i = 0; i < y.length; i++) {
        const xi = y[i], o = s.b0 * xi + z1;
        z1 = s.b1 * xi - s.a1 * o + z2; z2 = s.b2 * xi - s.a2 * o; y[i] = o;
      }
    });
    return y;
  }
  function filtfilt(secs, x) { let y = sosfilt(secs, x); y.reverse(); y = sosfilt(secs, y); y.reverse(); return y; }
  /** 漸縮（可選）→ 前後補零 → 零相位濾波 → 刪除補零 */
  function zp(x, secs, npad, nt) {
    const n = x.length, y = new Float64Array(n + 2 * npad);
    for (let i = 0; i < n; i++) {
      let w = 1;
      if (nt) { if (i < nt) w = 0.5 * (1 - Math.cos(Math.PI * i / nt)); else if (i >= n - nt) w = 0.5 * (1 - Math.cos(Math.PI * (n - 1 - i) / nt)); }
      y[npad + i] = x[i] * w;
    }
    return filtfilt(secs, y).slice(npad, npad + n);
  }
  function ctz(x, dt) { const o = new Float64Array(x.length); for (let i = 1; i < x.length; i++) o[i] = o[i - 1] + 0.5 * (x[i - 1] + x[i]) * dt; return o; }
  function mean(x) { let s = 0; for (const v of x) s += v; return s / x.length; }
  function absmax(x) { let m = 0, j = 0; for (let i = 0; i < x.length; i++) { const a = Math.abs(x[i]); if (a > m) { m = a; j = i; } } return { v: m, i: j }; }
  /** 0.1–30 Hz 帶通（4 階高通＋4 階低通串接，零相位） */
  function bandpass(a, dt, lo, hi) {
    const fs = 1 / dt, npad = Math.round(30 / dt), nt = Math.round(1 / dt);
    let y = zp(a, butter4('hp', lo || 0.1, fs), npad, nt);
    return zp(y, butter4('lp', hi || 30, fs), npad, 0);
  }
  /** 積分＋積分後 0.1 Hz 高通（HW1 步驟 5–6） */
  function integrate(a, dt, hp) {
    const fs = 1 / dt, npad = Math.round(30 / dt), sec = butter4('hp', 0.1, fs);
    let v = ctz(a, dt); if (hp) v = zp(v, sec, npad, 0);
    let d = ctz(v, dt); if (hp) d = zp(d, sec, npad, 0);
    return { v, d };
  }
  /** 基 2 FFT（實數輸入，補零到 2 的次方） */
  function fft(x, nfft) {
    let N = 1; while (N < (nfft || x.length)) N <<= 1;
    const re = new Float64Array(N), im = new Float64Array(N);
    re.set(x.length > N ? x.slice(0, N) : x);
    for (let i = 1, j = 0; i < N; i++) {
      let bit = N >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
      if (i < j) { const t = re[i]; re[i] = re[j]; re[j] = t; }
    }
    for (let len = 2; len <= N; len <<= 1) {
      const ang = -TAU / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < N; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2;
          const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
          re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
          const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
        }
      }
    }
    return { re, im, N };
  }
  /** FAS = Δt|X|（單位：加速度 m/s² → m/s） */
  function fas(a, dt, nfft) {
    const { re, im, N } = fft(a, nfft);
    const n2 = N / 2, f = new Float64Array(n2 + 1), A = new Float64Array(n2 + 1);
    for (let k = 0; k <= n2; k++) { f[k] = k / (N * dt); A[k] = dt * Math.hypot(re[k], im[k]); }
    return { f, A, N };
  }
  /** 直接 DFT 在指定頻率的 Δt|X(f)|（用於精修峰值） */
  function dftAt(a, dt, fq) {
    let r = 0, i = 0; const w = TAU * fq * dt;
    for (let n = 0; n < a.length; n++) { r += a[n] * Math.cos(w * n); i -= a[n] * Math.sin(w * n); }
    return dt * Math.hypot(r, i);
  }
  /** Konno–Ohmachi 平滑（b = 40） */
  function koSmooth(f, A, centers, b) {
    b = b || 40;
    return centers.map(fc => {
      const lo = fc * Math.pow(10, -Math.PI / b), hi = fc * Math.pow(10, Math.PI / b);
      let sw = 0, sa = 0;
      for (let k = 1; k < f.length; k++) {
        if (f[k] < lo) continue; if (f[k] > hi) break;
        const x = b * Math.log10(f[k] / fc), w = Math.abs(x) < 1e-10 ? 1 : Math.pow(Math.sin(x) / x, 4);
        sw += w; sa += w * A[k];
      }
      return sw > 0 ? sa / sw : NaN;
    });
  }
  /** 地動參數（與 B1 params() 相同的定義） */
  function params(a, v, d, dt) {
    const n = a.length; let cav = 0, s2 = 0; const H = new Float64Array(n);
    for (let i = 0; i < n; i++) { cav += Math.abs(a[i]) * dt; s2 += a[i] * a[i] * dt; H[i] = s2; }
    for (let i = 0; i < n; i++) H[i] /= s2;
    let i5 = 0; while (H[i5] < 0.05) i5++;
    let i95 = i5; while (H[i95] < 0.95) i95++;
    const F = fas(a, dt, 2 * n);
    let num = 0, den = 0;
    for (let k = 1; k < F.f.length; k++) if (F.f[k] >= 0.25 && F.f[k] <= 20) { const c2 = F.A[k] * F.A[k]; num += c2 / F.f[k]; den += c2; }
    const pa = absmax(a), pv = absmax(v), pd = absmax(d);
    return { PGA: pa.v / G, iA: pa.i, PGV: pv.v * 100, iV: pv.i, PGD: pd.v * 100, iD: pd.i, CAV: cav, AI: Math.PI / (2 * G) * s2,
      t5: i5 * dt, t95: i95 * dt, Ds: (i95 - i5) * dt, Tm: num / den, H };
  }
  return { TAU, G, butter4, sosfilt, filtfilt, zp, ctz, mean, absmax, bandpass, integrate, fft, fas, dftAt, koSmooth, params };
})();

document.addEventListener('sd:ready', () => {
  const S = SDOF, TAU = Math.PI * 2, G = DSP.G, f = SD.f;
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  const logspace = (a, b, n) => Array.from({ length: n }, (_, i) => a * Math.pow(b / a, i / (n - 1)));
  const deg = r => r * 180 / Math.PI, rad = d => d * Math.PI / 180;
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
  /** 畫圖用：每個像素格只留最大與最小值，保留峰值 */
  function thin(ts, ys, nb) {
    nb = nb || 1200;
    const n = ys.length;
    if (n <= 2 * nb) return [Array.from(ts), Array.from(ys)];
    const X = [], Y = [], step = n / nb;
    for (let b = 0; b < nb; b++) {
      const i0 = Math.floor(b * step), i1 = Math.min(n, Math.floor((b + 1) * step));
      let imn = i0, imx = i0;
      for (let i = i0; i < i1; i++) { if (ys[i] < ys[imn]) imn = i; if (ys[i] > ys[imx]) imx = i; }
      const [p, q] = imn < imx ? [imn, imx] : [imx, imn];
      X.push(ts[p], ts[q]); Y.push(ys[p], ys[q]);
    }
    return [X, Y];
  }
  const tvec = (n, dt) => { const t = new Float64Array(n); for (let i = 0; i < n; i++) t[i] = i * dt; return t; };
  const css = name => getComputedStyle(document.body).getPropertyValue(name).trim();

  /* ---------- PEER 紀錄 ---------- */
  const REC = window.CH02_REC || {};
  const STA = ['CHY010', 'HWA035', 'CHY065', 'CHY078'].filter(k => REC[k]);
  const staLabel = k => `${k}（R_rup ${REC[k].rrup.toFixed(0)} km、V_s30 ${REC[k].vs30.toFixed(0)} m/s）`.replace(/_/g, '');
  const compName = c => c === 'V' ? '垂直 U-D' : ({ N: '水平 N-S', S: '水平 N-S', E: '水平 E-W', W: '水平 E-W' }[c] || c);
  const recCache = {};
  function getComp(sta, code) {
    const key = sta + code;
    if (recCache[key]) return recCache[key];
    const r = REC[sta], q = r.comps[code], dt = r.dt;
    const a = Float64Array.from(q, v => v / 1000);
    const { v, d } = DSP.integrate(a, dt, true);
    const P = DSP.params(a, v, d, dt);
    return (recCache[key] = { a, v, d, P, dt, t: tvec(a.length, dt), code, sta });
  }
  const horiz = sta => Object.keys(REC[sta].comps).filter(c => c !== 'V');

  /* ================= W1 Mohr 圓 ================= */
  (() => {
    const body = SD.widget('w-mohr'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 's1', label: '最大主應力 σ_1 (MPa)', min: 20, max: 160, step: 1, value: 100, fmt: v => v.toFixed(0) },
      { k: 's3', label: '最小主應力 σ_3 (MPa)', min: 5, max: 120, step: 1, value: 40, fmt: v => v.toFixed(0) },
      { k: 'u', label: '孔隙水壓 u (MPa)', min: 0, max: 60, step: 0.5, value: 0, fmt: v => v.toFixed(1) },
      { k: 'mu', label: '摩擦係數 μ = tanφ', min: 0.3, max: 0.9, step: 0.01, value: 0.6, fmt: v => v.toFixed(2) }
    ].map(d => Object.assign(d, { label: SD.rich(d.label) })), draw);
    SD.buttons(body, '操作', [
      { label: '注水 +5 MPa', on: () => ctl.set({ u: Math.min(60, ctl.u + 5) }) },
      { label: '剛好觸發', on: () => ctl.set({ u: +crit().toFixed(1) }) },
      { label: '歸零', on: () => ctl.set({ u: 0 }) }
    ]);
    const p = new Plot(body, { w: 640, h: 330, xlim: [0, 170], ylim: [0, 80], xlabel: '有效正向應力 σ′ (MPa)', ylabel: '剪應力 τ (MPa)' });
    p.legend([{ label: '破壞包絡線 τ = μσ′', c: 'ink', w: 1.6 }, { label: '有效應力圓（扣掉 u）', c: 'blue', w: 2 }, { label: '總應力圓', c: 'gray', dash: true }, { label: '最接近破壞的面', c: 'red', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [
      { k: 'c', label: '圓心 (σ′_1+σ′_3)/2' }, { k: 'R', label: '半徑 (σ_1−σ_3)/2' }, { k: 'FS', label: '安全係數', hl: true },
      { k: 'uc', label: '觸發所需 u (MPa)', hl: true }, { k: 'st', label: '狀態' }
    ]);
    SD.div(body, 'w-note', '孔隙水壓不改變圓的大小（剪應力不受水壓影響），只把整個圓往左平移 u。圓碰到包絡線就破壞：對斷層是滑動（誘發地震），對土壤元素就是液化或剪力破壞。');
    function crit() {
      const sphi = ctl.mu / Math.sqrt(1 + ctl.mu * ctl.mu), s1 = Math.max(ctl.s1, ctl.s3), s3 = Math.min(ctl.s1, ctl.s3);
      return (s1 + s3) / 2 - (s1 - s3) / 2 / sphi;
    }
    function draw() {
      const s1 = Math.max(ctl.s1, ctl.s3), s3 = Math.min(ctl.s1, ctl.s3), u = ctl.u, mu = ctl.mu;
      const R = (s1 - s3) / 2, c = (s1 + s3) / 2 - u, ct = (s1 + s3) / 2;
      const sphi = mu / Math.sqrt(1 + mu * mu), cphi = 1 / Math.sqrt(1 + mu * mu);
      p.clear();
      const th = linspace(0, Math.PI, 181);
      p.line(th.map(t => ct + R * Math.cos(t)), th.map(t => R * Math.sin(t)), { c: 'gray', dash: '6 4', w: 1.2 });
      p.fn(x => mu * x, { c: 'ink', w: 1.6 });
      p.line(th.map(t => c + R * Math.cos(t)), th.map(t => R * Math.sin(t)), { c: 'blue', w: 2 });
      if (u > 0) { p.line([ct, c], [R * 0.5, R * 0.5], { c: 'gray', w: 1 }); p.text((ct + c) / 2, R * 0.5, '−u', { anchor: 'middle', dy: -6, cls: 'txt-muted' }); }
      // 最接近破壞的點：圓上法線與包絡線垂直
      const px = c - R * sphi, py = R * cphi;
      p.markers([px], [py], { c: 'red', fill: true, over: true });
      const FS = c > 0 ? c * sphi / R : 0, uc = crit();
      st({ c: f(c, 1), R: f(R, 1), FS: c > 0 ? f(FS, 3) : '—', uc: uc > 0 ? f(uc, 1) : '已破壞', st: FS > 1 ? '穩定' : '滑動！' }, { FS: FS <= 1, st: FS <= 1 });
    }
    draw();
  })();

  /* ================= W2 斷層與海灘球 ================= */
  (() => {
    const body = SD.widget('w-fault'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'phi', label: '走向 strike φ (°)', min: 0, max: 359, step: 1, value: 0, fmt: v => v.toFixed(0) },
      { k: 'dl', label: '傾角 dip δ (°)', min: 5, max: 90, step: 1, value: 45, fmt: v => v.toFixed(0) },
      { k: 'lam', label: '滑移角 rake λ (°)', min: -180, max: 180, step: 1, value: 90, fmt: v => v.toFixed(0) }
    ], draw);
    SD.buttons(body, '預設', [
      { label: '逆斷層', on: () => ctl.set({ phi: 0, dl: 45, lam: 90 }) },
      { label: '正斷層', on: () => ctl.set({ phi: 0, dl: 60, lam: -90 }) },
      { label: '右移', on: () => ctl.set({ phi: 0, dl: 90, lam: 180 }) },
      { label: '左移', on: () => ctl.set({ phi: 0, dl: 90, lam: 0 }) },
      { label: '逆衝＋左移（斜移）', on: () => ctl.set({ phi: 20, dl: 30, lam: 60 }) }
    ]);
    const wrap = SD.div(body, 'bb-wrap');
    const cv = document.createElement('canvas'); cv.width = 400; cv.height = 400;
    const cbox = SD.div(wrap, ''); cbox.appendChild(cv); SD.div(cbox, 'cap', '海灘球（下半球等面積投影）：深色＝壓縮（含 T 軸），淺色＝伸張（含 P 軸）');
    const W = 420, H = 280;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot-svg' }); wrap.appendChild(svg);
    const g = svgEl('g', null, svg);
    const st = SD.stats(body, [
      { k: 'type', label: '斷層類型', hl: true }, { k: 'ds', label: '傾移分量 sinλ（+ 上盤往上）' },
      { k: 'ss', label: '走向分量 cosλ（+ 左移）' }, { k: 'P', label: 'P 軸 傾伏/方位' }, { k: 'T', label: 'T 軸 傾伏/方位' }
    ]);
    SD.div(body, 'w-note', '左：海灘球。右上：地圖（北朝上），粗線是地表斷層線，短刺指向傾斜方向（上盤側）。右下：沿走向看的剖面，箭頭是上盤相對下盤的傾移運動。試著把 λ 從 90° 慢慢轉到 180°：逆斷層的中心陰影會轉成走向滑移的四象限。');
    function vecs(phi, dl, lam) {
      const [p, d, l] = [rad(phi), rad(dl), rad(lam)];
      const n = [-Math.sin(d) * Math.sin(p), Math.sin(d) * Math.cos(p), -Math.cos(d)];
      const s = [Math.cos(l) * Math.cos(p) + Math.cos(d) * Math.sin(l) * Math.sin(p), Math.cos(l) * Math.sin(p) - Math.cos(d) * Math.sin(l) * Math.cos(p), -Math.sin(l) * Math.sin(d)];
      return { n, s };
    }
    function axis(v) {
      let [x, y, z] = v; const L = Math.hypot(x, y, z); x /= L; y /= L; z /= L;
      if (z < 0) { x = -x; y = -y; z = -z; }
      const inc = Math.acos(Math.min(1, z)), az = (deg(Math.atan2(y, x)) + 360) % 360;
      const r = Math.SQRT2 * Math.sin(inc / 2);
      return { plunge: 90 - deg(inc), az, X: r * Math.sin(rad(az)), Y: -r * Math.cos(rad(az)) };
    }
    function hexRGB(s) {
      const m = s.match(/^#([0-9a-f]{6})$/i);
      if (m) { const v = parseInt(m[1], 16); return [v >> 16, (v >> 8) & 255, v & 255]; }
      const r = s.match(/\d+/g); return r ? r.slice(0, 3).map(Number) : [0, 0, 0];
    }
    function draw() {
      const { n, s } = vecs(ctl.phi, ctl.dl, ctl.lam);
      const ctx = cv.getContext('2d'), Wc = cv.width, R = Wc / 2 - 14, cx = Wc / 2, cy = Wc / 2;
      const img = ctx.createImageData(Wc, Wc);
      const ink = hexRGB(css('--ink')), bg = hexRGB(css('--inset') || '#ffffff');
      for (let j = 0; j < Wc; j++) for (let i = 0; i < Wc; i++) {
        const X = (i - cx) / R, Y = (j - cy) / R, rr = X * X + Y * Y, k = 4 * (j * Wc + i);
        if (rr > 1) { img.data[k + 3] = 0; continue; }
        const r = Math.sqrt(rr) * Math.SQRT2, inc = 2 * Math.asin(Math.min(1, r / Math.SQRT2 * 1));
        const incl = 2 * Math.asin(Math.min(1, Math.sqrt(rr) / Math.SQRT2 * Math.SQRT2 / Math.SQRT2));
        const az = Math.atan2(X, -Y);
        const sinI = Math.sin(incl), dir = [sinI * Math.cos(az), sinI * Math.sin(az), Math.cos(incl)];
        const amp = (dir[0] * n[0] + dir[1] * n[1] + dir[2] * n[2]) * (dir[0] * s[0] + dir[1] * s[1] + dir[2] * s[2]);
        const col = amp > 0 ? ink : bg;
        img.data[k] = col[0]; img.data[k + 1] = col[1]; img.data[k + 2] = col[2]; img.data[k + 3] = 255;
        void inc;
      }
      ctx.clearRect(0, 0, Wc, Wc);
      ctx.putImageData(img, 0, 0);
      ctx.strokeStyle = css('--ink'); ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(cx, cy, R, 0, TAU); ctx.stroke();
      const T = axis([n[0] + s[0], n[1] + s[1], n[2] + s[2]]), P = axis([n[0] - s[0], n[1] - s[1], n[2] - s[2]]);
      ctx.font = 'bold 22px ' + css('--font-text');
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      [[T, 'T'], [P, 'P']].forEach(([a, lab]) => {
        const x = cx + a.X * R, y = cy + a.Y * R;
        ctx.fillStyle = css('--red'); ctx.beginPath(); ctx.arc(x, y, 13, 0, TAU); ctx.fill();
        ctx.fillStyle = '#fff'; ctx.fillText(lab, x, y + 1);
      });
      ctx.fillStyle = css('--ink'); ctx.font = '20px ' + css('--font-text'); ctx.fillText('N', cx, 9);
      // 右側：地圖與剖面
      clearG(g);
      const ph = rad(ctl.phi), dl = rad(ctl.dl), lam = rad(ctl.lam), ss = Math.cos(lam), ds = Math.sin(lam);
      // 地圖
      const mx = 105, my = 92, L = 70, ux = Math.sin(ph), uy = -Math.cos(ph), rx = -uy, ry = ux; // 走向方向與右側（上盤）
      svgText(g, mx, 14, '(a) 地圖', 'txt-muted', 'middle');
      svgEl('line', { x1: mx - ux * L, y1: my - uy * L, x2: mx + ux * L, y2: my + uy * L, class: 's-ink', 'stroke-width': 2.4 }, g);
      for (let k = -2; k <= 2; k++) {
        const bx = mx + ux * k * 26, by = my + uy * k * 26;
        svgEl('line', { x1: bx, y1: by, x2: bx + rx * 8, y2: by + ry * 8, class: 's-ink', 'stroke-width': 1.6 }, g);
      }
      arrow(g, mx + ux * 50 + rx * 30, my + uy * 50 + ry * 30, mx + ux * 50 + rx * 30 + ux * 0, my + uy * 50 + ry * 30, 'ink');
      if (Math.abs(ss) > 0.05) {
        const a = 34 * ss;
        arrow(g, mx + rx * 26 - ux * a, my + ry * 26 - uy * a, mx + rx * 26 + ux * a, my + ry * 26 + uy * a, 'blue', 2.4);
        arrow(g, mx - rx * 26 + ux * a, my - ry * 26 + uy * a, mx - rx * 26 - ux * a, my - ry * 26 - uy * a, 'gray', 2.4);
      }
      svgText(g, mx + rx * 52, my + ry * 52 + 4, '上盤', 'txt-blue', 'middle');
      svgText(g, mx - rx * 52, my - ry * 52 + 4, '下盤', 'txt-muted', 'middle');
      svgEl('path', { d: `M${W / 2 - 18} 26 l-5 12 h10z`, class: 'f-ink' }, g);
      svgText(g, W / 2 - 18, 52, 'N', 'txt', 'middle');
      // 剖面（沿走向看，上盤在右）
      const sx0 = 250, sy0 = 40, Ls = 120, dxs = Math.cos(dl), dys = Math.sin(dl);
      svgText(g, 330, 14, '(b) 剖面（沿走向看）', 'txt-muted', 'middle');
      svgEl('line', { x1: 230, y1: sy0, x2: 410, y2: sy0, class: 's-gray', 'stroke-width': 1 }, g);
      svgEl('line', { x1: sx0 + 10, y1: sy0, x2: sx0 + 10 + dxs * Ls, y2: sy0 + dys * Ls, class: 's-ink', 'stroke-width': 2.4 }, g);
      const mxs = sx0 + 10 + dxs * Ls * 0.5, mys = sy0 + dys * Ls * 0.5, nx = -dys, ny = dxs; // 法線（往上盤為 -?）
      // 上盤在斷層上方右側：法線往右上 = (dys, -dxs)
      const hx = dys, hy = -dxs;
      if (Math.abs(ds) > 0.05) {
        const a = 30 * ds; // 正 = 上盤往上（逆）
        arrow(g, mxs + hx * 16 + dxs * a, mys + hy * 16 + dys * a, mxs + hx * 16 - dxs * a, mys + hy * 16 - dys * a, 'blue', 2.4);
        arrow(g, mxs - hx * 16 - dxs * a, mys - hy * 16 - dys * a, mxs - hx * 16 + dxs * a, mys - hy * 16 + dys * a, 'gray', 2.4);
      } else {
        svgText(g, mxs + 40, mys + 4, '⊙ ⊗ 水平錯動', 'txt-muted');
      }
      void nx; void ny;
      svgText(g, 395, 80, '上盤', 'txt-blue', 'end');
      svgText(g, 240, 150, '下盤', 'txt-muted');
      svgText(g, sx0 + 10 + 28, sy0 + 16, 'δ = ' + ctl.dl.toFixed(0) + '°', 'txt');
      svgText(g, 330, H - 50, ds > 0.05 ? '上盤往上 → 逆斷層成分' : ds < -0.05 ? '上盤往下 → 正斷層成分' : '幾乎沒有傾移', 'txt', 'middle');
      svgText(g, 330, H - 28, Math.abs(ss) > 0.05 ? (ss > 0 ? '看對面那塊往左 → 左移成分' : '看對面那塊往右 → 右移成分') : '幾乎沒有走向滑移', 'txt', 'middle');
      const l = ((ctl.lam + 540) % 360) - 180;
      let type;
      if (l >= 60 && l <= 120) type = '逆斷層';
      else if (l <= -60 && l >= -120) type = '正斷層';
      else if (Math.abs(l) <= 30) type = '左移走向滑移';
      else if (Math.abs(l) >= 150) type = '右移走向滑移';
      else type = (l > 0 ? '逆' : '正') + '＋' + (Math.abs(l) < 90 ? '左移' : '右移') + '（斜移）';
      const z0 = v => Math.abs(v) < 1e-9 ? 0 : v;   // 浮點誤差的 6e-17 顯示成 0
      st({ type, ds: f(z0(ds), 2), ss: f(z0(ss), 2), P: `${P.plunge.toFixed(0)}° / ${P.az.toFixed(0)}°`, T: `${T.plunge.toFixed(0)}° / ${T.az.toFixed(0)}°` });
    }
    draw();
    // 深淺色切換時重畫 canvas
    new MutationObserver(draw).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  })();

  /* ================= W3 彈性回跳 ================= */
  (() => {
    const body = SD.widget('w-rebound'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'v', label: '板塊相對速率 (cm/年)', min: 0.5, max: 8, step: 0.1, value: 2, fmt: v => v.toFixed(1) },
      { k: 'rate', label: '應力累積速率 (kPa/年)', min: 2, max: 40, step: 0.5, value: 10, fmt: v => v.toFixed(1) },
      { k: 'drop', label: '應力降 Δτ (MPa)', min: 1, max: 8, step: 0.1, value: 3, fmt: v => v.toFixed(1) },
      { k: 'du', label: '第 700 年注水 Δu (MPa)', min: 0, max: 4, step: 0.1, value: 0, fmt: v => v.toFixed(1) }
    ], draw);
    const row = SD.div(body, 'plot-row');
    const p1 = new Plot(row, { w: 420, h: 260, xlim: [0, 2000], ylim: [0, 1], xlabel: '時間 (年)', ylabel: '斷層剪應力 (MPa，相對值)' });
    const p2 = new Plot(row, { w: 420, h: 260, xlim: [0, 2000], ylim: [0, 40], xlabel: '時間 (年)', ylabel: '累積錯動 (m)' });
    p1.legend([{ label: '剪應力', c: 'blue', w: 2 }, { label: '斷層強度 μ(σ_n − u)', c: 'ink', dash: true }]);
    p2.legend([{ label: '斷層錯動（階梯）', c: 'blue', w: 2 }, { label: '遠處板塊位移 vt', c: 'gray', dash: true }]);
    const st = SD.stats(body, [
      { k: 'T', label: '再現週期 Δτ/τ̇ (年)', hl: true }, { k: 'D', label: '每次錯動 D ≈ vT (m)', hl: true },
      { k: 'n', label: '2000 年內地震次數' }, { k: 'Mw', label: 'M_w（破裂 40 km × 20 km、μ = 30 GPa）' }
    ]);
    SD.div(body, 'w-note', '模型：剪應力以固定速率累積，碰到強度就降一個 Δτ，期間累積的板塊位移一次釋放。拉動 Δu：強度線下移 0.6Δu（μ = 0.6），下一次地震提早，錯動量變小，這就是誘發地震「提早觸發」的效果，並不會憑空多出能量。');
    function draw() {
      const v = ctl.v / 100, rate = ctl.rate / 1000, drop = ctl.drop, du = ctl.du;
      const base = 10; // 強度基準 (MPa，相對)
      let tau = base - drop * 0.6, last = 0, slip = 0;
      const ts = [0], ys = [tau], ss = [0], st2 = [0], ev = [];
      const dt = 1;
      for (let t = 1; t <= 2000; t += dt) {
        tau += rate * dt;
        const str = base - (t >= 700 ? 0.6 * du : 0);
        ts.push(t); ys.push(tau); ss.push(slip); st2.push(t);
        if (tau >= str) {
          tau = str - drop; slip += v * (t - last); ev.push({ t, D: v * (t - last) }); last = t;
          ts.push(t); ys.push(tau); ss.push(slip); st2.push(t);
        }
      }
      const ymin = base - drop - 1.5, ymax = base + 1;
      p1.setLim(null, [ymin, ymax]); p1.clear();
      p1.line([0, 700, 700, 2000], [base, base, base - 0.6 * du, base - 0.6 * du], { c: 'ink', dash: '6 4', w: 1.4 });
      p1.line(ts, ys, { c: 'blue', w: 1.8 });
      if (du > 0) p1.vline(700);
      const smax = Math.max(v * 2000, slip) * 1.08;
      p2.setLim(null, [0, smax]); p2.clear();
      p2.line([0, 2000], [0, v * 2000], { c: 'gray', dash: '6 4', w: 1.4 });
      p2.line(st2, ss, { c: 'blue', w: 2 });
      const T = drop / rate, D = v * T, M0 = 3e10 * 40e3 * 20e3 * D;
      st({ T: f(T, 0), D: f(D, 2), n: String(ev.length), Mw: f((Math.log10(M0) - 9.05) / 1.5, 2) });
    }
    draw();
  })();

  /* ================= W4 震度計 ================= */
  (() => {
    const body = SD.widget('w-intensity'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'pga', label: 'PGA (cm/s²)', min: 0.3, max: 1500, log: true, snap: 3, value: 120, fmt: v => v.toPrecision(3) },
      { k: 'pgv', label: 'PGV (cm/s)', min: 0.1, max: 300, log: true, snap: 3, value: 35, fmt: v => v.toPrecision(3) }
    ], draw);
    SD.buttons(body, '範例', [
      { label: '講義例 120、35', on: () => ctl.set({ pga: 120, pgv: 35 }) },
      { label: 'CHY010（RotD100）', on: () => ctl.set({ pga: 247, pgv: 21.4 }) },
      { label: 'HWA035（RotD100）', on: () => ctl.set({ pga: 81.5, pgv: 9.08 }) },
      { label: '小區域強搖 400、12', on: () => ctl.set({ pga: 400, pgv: 12 }) }
    ]);
    const p = new Plot(body, { w: 640, h: 360, xlim: [0.3, 1500], ylim: [0.1, 300], logx: true, logy: true, xlabel: 'PGA (cm/s²)', ylabel: 'PGV (cm/s)' });
    p.legend([{ label: '分級界線', c: 'ink', w: 1.2 }, { label: '目前的測站', c: 'blue', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [{ k: 'lv', label: '新制震度', hl: true }, { k: 'by', label: '判定依據' }, { k: 'pgaonly', label: '若只看 PGA 對照' }, { k: 'g', label: 'PGA (g)' }]);
    const note = SD.div(body, 'w-note');
    const names = ['0 級', '1 級', '2 級', '3 級', '4 級', '5 弱', '5 強', '6 弱', '6 強', '7 級'];
    const pgaB = [0.8, 2.5, 8, 25, 80, 140, 250, 440, 800], pgvB = [15, 30, 50, 80, 140];
    function level(pga, pgv) {
      if (pga < 80) { let k = 0; while (k < 4 && pga >= pgaB[k]) k++; return { k, by: 'PGA（< 80 cm/s²）' }; }
      if (pgv < 15) return { k: 4, by: 'PGA ≥ 80 但 PGV < 15 → 4 級' };
      let k = 5; for (let j = 1; j < pgvB.length; j++) if (pgv >= pgvB[j]) k = 5 + j;
      return { k, by: 'PGV（PGA ≥ 80）' };
    }
    function draw() {
      p.clear();
      // 左區：PGA 分級（垂直線）；右區：PGV 分級（水平線）
      pgaB.slice(0, 5).forEach(x => p.line([x, x], [0.1, x >= 80 ? 15 : 300], { c: 'ink', w: 1.1 }));
      p.line([80, 80], [0.1, 300], { c: 'ink', w: 1.6 });
      pgvB.forEach(y => p.line([80, 1500], [y, y], { c: 'ink', w: 1.1 }));
      const lab = [[0.45, '0'], [1.4, '1'], [4.5, '2'], [14, '3'], [45, '4']];
      lab.forEach(([x, s]) => p.text(x, 0.16, s, { anchor: 'middle', cls: 'txt-muted' }));
      p.text(300, 4, '4 級（PGV < 15）', { anchor: 'middle', cls: 'txt-muted' });
      [[21, '5 弱'], [39, '5 強'], [63, '6 弱'], [106, '6 強'], [200, '7']].forEach(([y, s]) => p.text(700, y, s, { anchor: 'middle', dy: 4, cls: 'txt-muted' }));
      p.text(90, 220, '≥ 80 改看 PGV →', { cls: 'txt-muted' });
      p.markers([ctl.pga], [ctl.pgv], { c: 'blue', fill: true, over: true, r: 5.5 });
      const L = level(ctl.pga, ctl.pgv);
      let ka = 0; while (ka < 9 && ctl.pga >= pgaB[ka]) ka++;
      st({ lv: names[L.k], by: L.by, pgaonly: names[ka], g: f(ctl.pga / 981, 3) }, { pgaonly: ka !== L.k });
      note.textContent = ka !== L.k ? `只看 PGA 會判成 ${names[ka]}，與新制的 ${names[L.k]} 不同：短週期的加速度尖峰對建物損害的代表性不如速度。` : 'PGA 與 PGV 給出同一級。';
      document.querySelectorAll('#cwa-table th').forEach((th, i) => th.classList.toggle('on', i === L.k + 1));
    }
    draw();
  })();

  /* ================= W5 規模計算機 ================= */
  (() => {
    const body = SD.widget('w-mag'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'mu', label: '剛性模數 μ (GPa)', min: 10, max: 50, step: 1, value: 30, fmt: v => v.toFixed(0) },
      { k: 'L', label: '破裂長度 L (km)', min: 1, max: 1000, log: true, snap: 3, value: 80, fmt: v => v.toPrecision(3) },
      { k: 'W', label: '破裂寬度 W (km)', min: 1, max: 200, log: true, snap: 3, value: 40, fmt: v => v.toPrecision(3) },
      { k: 'D', label: '平均滑移 D (m)', min: 0.01, max: 40, log: true, snap: 3, value: 3, fmt: v => v.toPrecision(3) },
      { k: 'M2', label: '比較的地震 M_w', min: 4, max: 9.5, step: 0.1, value: 7.4, fmt: v => v.toFixed(1) }
    ].map(d => Object.assign(d, { label: SD.rich(d.label) })), draw);
    SD.buttons(body, '範例', [
      { label: '集集（講義例）', on: () => ctl.set({ mu: 30, L: 80, W: 40, D: 3, M2: 7.4 }) },
      { label: '測驗：20×10 km、0.8 m', on: () => ctl.set({ mu: 30, L: 20, W: 10, D: 0.8 }) },
      { label: '東日本 量級 500×200 km、10 m', on: () => ctl.set({ mu: 40, L: 500, W: 200, D: 10, M2: 7.6 }) }
    ]);
    const st = SD.stats(body, [
      { k: 'A', label: '破裂面積 A (m²)' }, { k: 'M0', label: 'M_0 = μAD (N·m)' }, { k: 'Mw', label: 'M_w（9.05）', hl: true },
      { k: 'Mw2', label: 'M_w（9.1，另一套常數）' }, { k: 'E', label: '輻射能量 ≈ M_0/2×10⁴ (J)' }, { k: 'r', label: '與比較地震的 M_0 比', hl: true }
    ]);
    const p = new Plot(body, { w: 640, h: 330, xlim: [4, 9.6], ylim: [4, 9.6], xlabel: '地震矩規模 M_w', ylabel: '量到的規模' });
    p.legend([{ label: 'M_w（不飽和）', c: 'ink', w: 1.4 }, { label: 'm_b（約 1 s P 波）', c: 'blue', w: 2 }, { label: 'M_L（0.1–1 s）', c: 'red', w: 2 }, { label: 'M_s（約 20 s）', c: 'gray', w: 2 }, { label: '臺灣實例 (M_w, M_L)', c: 'ink', line: false, marker: 'circle' }]);
    SD.div(body, 'w-note', '下圖的三條飽和曲線是示意（平滑地趨近 m_b ≈ 6.3、M_L ≈ 7.4、M_s ≈ 8.3），用來看出「量測週期越短，越早飽和」；空心點是集集（7.6, 7.3）與 2024 花蓮（7.4, 7.2）的實際值。');
    const sat = (M, Ms, k) => M - Math.log(1 + Math.exp(k * (M - Ms))) / k;
    function draw() {
      const A = ctl.L * 1e3 * ctl.W * 1e3, M0 = ctl.mu * 1e9 * A * ctl.D;
      const Mw = (Math.log10(M0) - 9.05) / 1.5, M02 = Math.pow(10, 1.5 * ctl.M2 + 9.05);
      const r = M0 / M02;
      st({ A: f(A, 3), M0: M0.toExponential(2), Mw: f(Mw, 2), Mw2: f((Math.log10(M0) - 9.1) / 1.5, 2), E: (M0 / 2e4).toExponential(2),
        r: r >= 1 ? f(r, r < 10 ? 2 : 0) + ' 倍' : '1/' + f(1 / r, 1 / r < 10 ? 2 : 0) });
      p.clear();
      p.fn(M => M, { c: 'ink', w: 1.4 });
      p.fn(M => sat(M, 6.3, 4), { c: 'blue', w: 2 });
      p.fn(M => sat(M, 7.4, 6), { c: 'red', w: 2 });
      p.fn(M => sat(M, 8.3, 4), { c: 'gray', w: 2 });
      p.markers([7.6, 7.4], [7.3, 7.2], { c: 'ink', over: true });
      if (Mw >= 4 && Mw <= 9.6) p.vline(Mw);
      p.text(Math.min(9.4, Math.max(4.1, Mw)), 4.25, '目前 M_w = ' + Mw.toFixed(2), { anchor: Mw > 8 ? 'end' : 'start', dx: Mw > 8 ? -6 : 6, cls: 'txt-muted' });
    }
    draw();
  })();

  /* ================= W6 紀錄檢視器 ================= */
  (() => {
    const body = SD.widget('w-record'); if (!body || !STA.length) return;
    const ctl = SD.controls(body, [
      { k: 'sta', label: '測站', type: 'select', value: STA[0], options: STA.map(k => [k, staLabel(k)]) },
      { k: 'c', label: '分量', type: 'select', value: '0', options: [['0', '水平 1'], ['1', '水平 2'], ['2', '垂直 U-D']] }
    ], draw);
    SD.div(body, 'w-note', 'PEER NGA-West2 的 .AT2 紀錄，已在 Python 以 3.3 的流程（去平均、1 s 漸縮、0.1–30 Hz 零相位帶通）處理並降為 100 Hz 嵌入；網頁再做梯形積分與積分後 0.1 Hz 高通，並即時計算所有參數。');
    const pa = new Plot(body, { h: 170, xlim: [0, 90], ylim: [-0.3, 0.3], ylabel: 'a (g)', m: { b: 28 } });
    const pv = new Plot(body, { h: 170, xlim: [0, 90], ylim: [-20, 20], ylabel: 'v (cm/s)', m: { b: 28 } });
    const pd = new Plot(body, { h: 190, xlim: [0, 90], ylim: [-5, 5], ylabel: 'd (cm)', xlabel: '時間 t (s)' });
    const ph = new Plot(body, { h: 200, xlim: [0, 90], ylim: [0, 1], ylabel: 'Husid H(t)', xlabel: '時間 t (s)' });
    pa.legend([{ label: '紀錄（實測）', c: 'ink', w: 1 }, { label: '峰值', c: 'red', line: false, marker: 'circle', fill: true }]);
    ph.legend([{ label: '正規化 Arias intensity', c: 'blue', w: 2 }, { label: '強震段 t_5–t_95', c: 'gray', w: 8, op: .4 }]);
    const st = SD.stats(body, [
      { k: 'pga', label: 'PGA (g)', hl: true }, { k: 'pgv', label: 'PGV (cm/s)', hl: true }, { k: 'pgd', label: 'PGD (cm)' }, { k: 'ratio', label: 'PGV/PGA (s)' },
      { k: 'ai', label: 'I_a (m/s)' }, { k: 'cav', label: 'CAV (m/s)' }, { k: 'ds', label: 'D_{5–95} (s)' }, { k: 'tm', label: 'T_m (s)' }
    ]);
    const cmp = SD.div(body, 'w-note');
    function draw() {
      const sta = ctl.sta, codes = Object.keys(REC[sta].comps), code = codes[+ctl.c] || codes[0];
      const R = getComp(sta, code), P = R.P, T = R.t[R.t.length - 1];
      const ag = R.a.map(x => x / G), vc = R.v.map(x => x * 100), dc = R.d.map(x => x * 100);
      const lim = (x, m) => { const v = m * 1.15; return [-v, v]; };
      [pa, pv, pd, ph].forEach(p => p.setLim([0, Math.ceil(T / 10) * 10], null));
      pa.setLim(null, lim(ag, P.PGA)); pv.setLim(null, lim(vc, P.PGV)); pd.setLim(null, lim(dc, P.PGD));
      [[pa, ag, P.iA], [pv, vc, P.iV], [pd, dc, P.iD]].forEach(([p, y, i]) => {
        p.clear(); p.hline(0, { dash: null, w: 0.6 });
        const [X, Y] = thin(R.t, y); p.line(X, Y, { c: 'ink', w: 0.9 });
        p.markers([R.t[i]], [y[i]], { c: 'red', fill: true, over: true });
      });
      ph.clear(); ph.band(P.t5, P.t95);
      const [X, Y] = thin(R.t, P.H, 600); ph.line(X, Y, { c: 'blue', w: 2 });
      ph.hline(0.05); ph.hline(0.95);
      ph.text(P.t95, 0.5, ' D_{5–95} = ' + P.Ds.toFixed(1) + ' s', { cls: 'txt' });
      st({ pga: f(P.PGA, 4), pgv: f(P.PGV, 2), pgd: f(P.PGD, 2), ratio: f(P.PGV / (P.PGA * 981), 3), ai: f(P.AI, 3), cav: f(P.CAV, 2), ds: f(P.Ds, 2), tm: f(P.Tm, 3) });
      const ref = (REC[sta].full || {})[code];
      cmp.innerHTML = `分量 ${code}（${compName(code)}）。以原始取樣率算出的值：PGA ${ref ? ref.PGA.toFixed(4) : '—'} g、PGV ${ref ? ref.PGV.toFixed(2) : '—'} cm/s、PGD ${ref ? ref.PGD.toFixed(2) : '—'} cm、I<sub>a</sub> ${ref ? ref.AI.toFixed(3) : '—'} m/s、D<sub>5–95</sub> ${ref ? ref.Ds.toFixed(2) : '—'} s，與上方網頁即時值的差異來自降取樣，都在 1 % 以內。注意 PGA、PGV、PGD 出現在不同時刻：越積分越偏長週期。`;
    }
    draw();
  })();

  /* ================= W7 延時（合成） ================= */
  (() => {
    const body = SD.widget('w-dur'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'Ts', label: '強震段長度 (s)', min: 2, max: 30, step: 0.5, value: 15, fmt: v => v.toFixed(1) },
      { k: 'fq', label: '主要頻率 (Hz)', min: 0.5, max: 8, step: 0.1, value: 2, fmt: v => v.toFixed(1) },
      { k: 'pga', label: 'PGA (g)', min: 0.05, max: 0.6, step: 0.01, value: 0.3, fmt: v => v.toFixed(2) }
    ], draw);
    SD.div(body, 'w-note', '合成地動（教學用，不是真實紀錄）：a(t) = 包絡線 × 正弦，包絡線 2 s 漸強、維持「強震段長度」、再指數衰減，並把峰值調成指定的 PGA。灰色是固定的參考地動（強震段 5 s、2 Hz、PGA 相同）。');
    const pa = new Plot(body, { h: 200, xlim: [0, 45], ylim: [-0.35, 0.35], ylabel: 'a (g)', m: { b: 28 } });
    const ph = new Plot(body, { h: 220, xlim: [0, 45], ylim: [0, 1], ylabel: 'Husid H(t)', xlabel: '時間 t (s)' });
    pa.legend([{ label: '合成地動', c: 'blue', w: 1 }, { label: '參考地動（強震段 5 s）', c: 'gray', w: 1 }]);
    ph.legend([{ label: '合成地動', c: 'blue', w: 2 }, { label: '參考地動', c: 'gray', w: 2 }]);
    const st = SD.stats(body, [
      { k: 'pga', label: 'PGA (g)' }, { k: 'ia', label: 'I_a (m/s)', hl: true }, { k: 'iar', label: 'I_a / 參考' }, { k: 'cav', label: 'CAV (m/s)', hl: true },
      { k: 'ds', label: 'D_{5–95} (s)' }, { k: 'nc', label: '強震段循環數 ≈ f·D_{5–95}' }, { k: 'pgv', label: 'PGV (cm/s)' }
    ]);
    const dt = 0.01, N = 4500, t = tvec(N, dt);
    function motion(Ts, fq, pga) {
      const a = new Float64Array(N);
      for (let i = 0; i < N; i++) {
        const tt = t[i], env = tt < 2 ? (tt / 2) ** 2 : tt < 2 + Ts ? 1 : Math.exp(-0.6 * (tt - 2 - Ts));
        a[i] = env * Math.sin(TAU * fq * tt);
      }
      const m = DSP.absmax(a).v;
      for (let i = 0; i < N; i++) a[i] *= pga * G / m;
      const { v, d } = DSP.integrate(a, dt, true);
      return { a, P: DSP.params(a, v, d, dt) };
    }
    function draw() {
      const ref = motion(5, 2, ctl.pga), M = motion(ctl.Ts, ctl.fq, ctl.pga);
      pa.setLim(null, [-ctl.pga * 1.15, ctl.pga * 1.15]);
      pa.clear(); ph.clear();
      let [X, Y] = thin(t, ref.a.map(x => x / G)); pa.line(X, Y, { c: 'gray', w: 0.8 });
      [X, Y] = thin(t, M.a.map(x => x / G)); pa.line(X, Y, { c: 'blue', w: 0.9 });
      ph.band(M.P.t5, M.P.t95);
      [X, Y] = thin(t, ref.P.H, 500); ph.line(X, Y, { c: 'gray', w: 2 });
      [X, Y] = thin(t, M.P.H, 500); ph.line(X, Y, { c: 'blue', w: 2 });
      ph.hline(0.05); ph.hline(0.95);
      st({ pga: f(M.P.PGA, 2), ia: f(M.P.AI, 3), iar: f(M.P.AI / ref.P.AI, 2) + ' 倍', cav: f(M.P.CAV, 2), ds: f(M.P.Ds, 1), nc: f(ctl.fq * M.P.Ds, 0), pgv: f(M.P.PGV, 1) });
    }
    draw();
  })();

  /* ================= W8 Fourier 與取樣 ================= */
  (() => {
    const body = SD.widget('w-fourier'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'f1', label: 'f_1 (Hz)', min: 0.5, max: 120, step: 0.5, value: 5, fmt: v => v.toFixed(1) },
      { k: 'a1', label: '振幅 A_1', min: 0, max: 10, step: 0.5, value: 10, fmt: v => v.toFixed(1) },
      { k: 'f2', label: 'f_2 (Hz)', min: 0.5, max: 120, step: 0.5, value: 14, fmt: v => v.toFixed(1) },
      { k: 'a2', label: '振幅 A_2', min: 0, max: 10, step: 0.5, value: 10, fmt: v => v.toFixed(1) },
      { k: 'f3', label: 'f_3 (Hz)', min: 0.5, max: 120, step: 0.5, value: 60, fmt: v => v.toFixed(1) },
      { k: 'a3', label: '振幅 A_3', min: 0, max: 10, step: 0.5, value: 0, fmt: v => v.toFixed(1) },
      { k: 'fs', label: '取樣率 f_s (Hz)', min: 10, max: 250, step: 1, value: 100, fmt: v => v.toFixed(0) },
      { k: 'T', label: '紀錄長度 T (s)', min: 0.5, max: 20, step: 0.5, value: 4, fmt: v => v.toFixed(1) }
    ].map(d => Object.assign(d, { label: SD.rich(d.label) })), draw);
    SD.buttons(body, '情境', [
      { label: '講義：幾個等振幅正弦', on: () => ctl.set({ f1: 2, a1: 10, f2: 5, a2: 10, f3: 9, a3: 10, fs: 100, T: 4 }) },
      { label: '手機量 60 Hz', on: () => ctl.set({ f1: 14, a1: 10, f2: 28, a2: 0, f3: 60, a3: 10, fs: 100, T: 4 }) },
      { label: '紀錄太短：5 與 5.5 Hz', on: () => ctl.set({ f1: 5, a1: 10, f2: 5.5, a2: 10, f3: 60, a3: 0, fs: 100, T: 1 }) }
    ]);
    const pt = new Plot(body, { h: 230, xlim: [0, 0.5], ylim: [-30, 30], xlabel: '時間 t (s)（只顯示前 0.5 s）', ylabel: 'a(t)' });
    pt.legend([{ label: '真實的連續訊號', c: 'gray', w: 1 }, { label: '取樣點（量到的資料）', c: 'ink', dash: true, marker: 'circle' }]);
    const pf = new Plot(body, { h: 260, xlim: [0, 50], ylim: [0, 12], xlabel: '頻率 f (Hz)', ylabel: '振幅 2|X_k|/N' });
    pf.legend([{ label: '由取樣點算的振幅譜', c: 'blue', w: 1.6 }, { label: '真實頻率', c: 'gray', dash: true }, { label: '混疊後出現的位置', c: 'red', line: false, marker: 'circle', fill: true }]);
    const st = SD.stats(body, [
      { k: 'fn', label: 'Nyquist f_N = f_s/2 (Hz)', hl: true }, { k: 'df', label: '解析度 Δf = 1/T (Hz)' }, { k: 'N', label: '點數 N' },
      { k: 'al', label: '各成分出現在 (Hz)', hl: true }
    ]);
    SD.div(body, 'w-note', '縱軸用 2|X_k|/N，可以直接讀成正弦波的振幅；FAS = Δt|X_k| 只差一個常數。把 f_3 拉過 Nyquist：真實頻率（灰虛線）跑出圖外，譜上卻在 |f − f_s·round(f/f_s)| 冒出一根假峰。把紀錄縮短到 1 s：5 與 5.5 Hz 擠成一個峰，因為 Δf = 1 Hz 分不開它們。');
    function draw() {
      const comps = [[ctl.f1, ctl.a1], [ctl.f2, ctl.a2], [ctl.f3, ctl.a3]].filter(c => c[1] > 0);
      const sig = tt => comps.reduce((s, [fq, A]) => s + A * Math.sin(TAU * fq * tt), 0);
      const fs = ctl.fs, dt = 1 / fs, N = Math.max(4, Math.round(ctl.T * fs));
      const xs = new Float64Array(N), ts = tvec(N, dt);
      for (let i = 0; i < N; i++) xs[i] = sig(ts[i]);
      const amp = comps.reduce((s, c) => s + c[1], 0) || 1;
      pt.setLim(null, [-amp * 1.1, amp * 1.1]); pt.clear();
      pt.fn(sig, { c: 'gray', w: 1, n: 1500 });
      const k = ts.findIndex(v => v > 0.5); const n2 = k < 0 ? N : k;
      pt.line(Array.from(ts.slice(0, n2)), Array.from(xs.slice(0, n2)), { c: 'ink', dash: '5 4', w: 1 });
      if (n2 <= 130) pt.markers(Array.from(ts.slice(0, n2)), Array.from(xs.slice(0, n2)), { c: 'ink', r: 3.4 });
      const F = DSP.fft(xs, Math.max(4096, N));
      const fr = [], A = [];
      for (let j = 0; j <= F.N / 2; j++) { fr.push(j * fs / F.N); A.push(2 * Math.hypot(F.re[j], F.im[j]) / N); }
      pf.setLim([0, fs / 2], [0, Math.max(...A, 1) * 1.15]); pf.clear();
      const al = comps.map(([fq]) => Math.abs(fq - fs * Math.round(fq / fs)));
      comps.forEach(([fq]) => { if (fq <= fs / 2) pf.vline(fq); });
      pf.line(fr, A, { c: 'blue', w: 1.6 });
      comps.forEach(([fq, a], i) => { if (Math.abs(al[i] - fq) > 1e-9) pf.markers([al[i]], [a], { c: 'red', fill: true, over: true }); });
      st({ fn: f(fs / 2, 1), df: f(1 / (N * dt), 3), N: String(N),
        al: comps.map(([fq], i) => fq.toFixed(1) + '→' + al[i].toFixed(1)).join('<br>') || '—' },
        { al: comps.some(([fq]) => fq > fs / 2) });
    }
    draw();
  })();

  /* ================= W9 岩盤 vs 軟土 ================= */
  (() => {
    const body = SD.widget('w-site'); if (!body || !STA.length) return;
    const ctl = SD.controls(body, [
      { k: 'A', label: '測站 A', type: 'select', value: 'HWA035', options: STA.map(k => [k, staLabel(k)]) },
      { k: 'B', label: '測站 B', type: 'select', value: 'CHY078', options: STA.map(k => [k, staLabel(k)]) },
      { k: 'c', label: '分量', type: 'select', value: '0', options: [['0', '水平 1'], ['1', '水平 2']] },
      { k: 'H', label: '土層厚度 H (m)', min: 5, max: 120, step: 1, value: 40, fmt: v => v.toFixed(0) },
      { k: 'Vs', label: '土層 V_s (m/s)', min: 80, max: 800, step: 5, value: 250, fmt: v => v.toFixed(0) }
    ].map(d => Object.assign(d, { label: SD.rich(d.label) })), draw);
    const row = SD.div(body, 'plot-row');
    const pf = new Plot(row, { w: 420, h: 300, xlim: [0.04, 10], ylim: [0.005, 2], logx: true, logy: true, xlabel: '週期 T = 1/f (s)', ylabel: 'FAS (m/s)，Konno–Ohmachi 平滑' });
    const ps = new Plot(row, { w: 420, h: 300, xlim: [0.04, 5], ylim: [0, 1], logx: true, xlabel: '自然週期 T (s)', ylabel: 'S_a (g)，ξ = 5 %' });
    pf.legend([{ label: '測站 A', c: 'ink', dash: true }, { label: '測站 B', c: 'red', dash: true }, { label: '峰值 T_p', c: 'ink', line: false, marker: 'circle' }]);
    ps.legend([{ label: '測站 A', c: 'ink', w: 2 }, { label: '測站 B', c: 'red', w: 2 }, { label: '土層 T = 4H/V_s', c: 'gray', dash: true }]);
    const st = SD.stats(body, [
      { k: 'tpA', label: 'A：FAS 峰值週期 (s)' }, { k: 'tpB', label: 'B：FAS 峰值週期 (s)' }, { k: 'saA', label: 'A：S_a 峰值週期 (s)' }, { k: 'saB', label: 'B：S_a 峰值週期 (s)' },
      { k: 'tmA', label: 'A：T_m (s)' }, { k: 'tmB', label: 'B：T_m (s)' }, { k: 'Ts', label: '土層 4H/V_s (s)', hl: true }
    ]);
    SD.div(body, 'w-note', '兩種「卓越週期」：FAS 的峰（地動本身的頻率內容）與 5 % 阻尼反應譜的峰（講義 3.3 表的 T_p 定義）。HWA035（V_s30 677 m/s）的能量集中在 0.1–0.5 s；CHY078（161 m/s）與 CHY065（250 m/s）明顯往長週期移。拖動 H、V_s 看看哪種土層週期會和哪一站的峰對上。距離不同也會改變振幅，比較時以「形狀」為主。');
    const fc = logspace(0.1, 25, 140), Ts = logspace(0.04, 5, 48);
    const cache = {};
    function spec(sta, ci) {
      const key = sta + ci; if (cache[key]) return cache[key];
      const code = horiz(sta)[ci] || horiz(sta)[0], R = getComp(sta, code);
      const F = DSP.fas(R.a, R.dt), sm = DSP.koSmooth(F.f, F.A, fc);
      let ip = 0; sm.forEach((v, i) => { if (v > sm[ip]) ip = i; });
      const sa = Ts.map(T => (TAU / T) ** 2 * S.peakDisp(R.a, R.dt, T, 0.05) / G);
      let is = 0; sa.forEach((v, i) => { if (v > sa[is]) is = i; });
      return (cache[key] = { sm, tp: 1 / fc[ip], fp: sm[ip], sa, ts: Ts[is], sap: sa[is], Tm: R.P.Tm });
    }
    function draw() {
      const ci = +ctl.c, A = spec(ctl.A, ci), B = spec(ctl.B, ci), per = fc.map(x => 1 / x);
      pf.clear(); ps.clear();
      const Tsoil = 4 * ctl.H / ctl.Vs;
      pf.vline(Tsoil); ps.vline(Tsoil);
      pf.line(per, A.sm, { c: 'ink', dash: '6 4', w: 1.6 }); pf.line(per, B.sm, { c: 'red', dash: '6 4', w: 1.6 });
      pf.markers([A.tp], [A.fp], { c: 'ink', over: true }); pf.markers([B.tp], [B.fp], { c: 'red', over: true });
      ps.setLim(null, [0, Math.max(...A.sa, ...B.sa) * 1.12]);
      ps.line(Ts, A.sa, { c: 'ink', w: 2 }); ps.line(Ts, B.sa, { c: 'red', w: 2 });
      st({ tpA: f(A.tp, 2), tpB: f(B.tp, 2), saA: f(A.ts, 2), saB: f(B.ts, 2), tmA: f(A.Tm, 2), tmB: f(B.Tm, 2), Ts: f(Tsoil, 2) });
    }
    draw();
  })();

  /* ================= W10 手機積分漂移 ================= */
  const PHONE = (() => {
    const d = window.CH02_PHONE; if (!d) return null;
    const out = { dt: d.dt, n: d.n, t: tvec(d.n, d.dt), ax: {} };
    out.get = k => {
      if (out.ax[k]) return out.ax[k];
      const raw = Float64Array.from(d[k], v => v / 1000), dt = d.dt, eps = DSP.mean(raw);
      const L0 = DSP.integrate(raw, dt, false);
      const dm = raw.map(v => v - eps), L1 = DSP.integrate(dm, dt, false);
      const bp = DSP.bandpass(dm, dt), L2 = DSP.integrate(bp, dt, false), L3 = DSP.integrate(bp, dt, true);
      return (out.ax[k] = { raw, eps, dm, bp, L: [{ a: raw, ...L0 }, { a: dm, ...L1 }, { a: bp, ...L2 }, { a: bp, ...L3 }] });
    };
    return out;
  })();
  (() => {
    const body = SD.widget('w-drift'); if (!body || !PHONE) return;
    const ctl = SD.controls(body, [
      { k: 'ax', label: '方向', type: 'select', value: 'z', options: [['x', 'X 左右'], ['y', 'Y 前後'], ['z', 'Z 上下']] },
      { k: 'lv', label: '處理到哪一步', type: 'select', value: '0', options: [['0', '① 原始資料直接積分'], ['1', '② 去平均後積分'], ['2', '③ 去平均＋0.1–30 Hz 帶通後積分'], ['3', '④ 完整流程（再加積分後高通）']] }
    ], draw);
    SD.buttons(body, '逐步', ['①', '②', '③', '④'].map((s, i) => ({ label: s, on: () => ctl.set({ lv: String(i) }) })));
    const pa = new Plot(body, { h: 180, xlim: [0, 162], ylim: [-1, 1], ylabel: 'a (m/s²)', m: { b: 28 } });
    const pv = new Plot(body, { h: 180, xlim: [0, 162], ylim: [-1, 1], ylabel: 'v (m/s)', m: { b: 28 } });
    const pd = new Plot(body, { h: 200, xlim: [0, 162], ylim: [-1, 1], ylabel: 'u', xlabel: '時間 t (s)' });
    pa.legend([{ label: '此步驟使用的加速度', c: 'ink', w: 1 }, { label: '平均值 ε', c: 'red', dash: true }]);
    pv.legend([{ label: '積分結果', c: 'blue', w: 1.4 }, { label: '偏差預測 εt、½εt²（①）', c: 'gray', dash: true }]);
    const st = SD.stats(body, [
      { k: 'eps', label: '平均值 ε (m/s²)' }, { k: 'vend', label: 'v 終值 (m/s)' }, { k: 'uend', label: 'u 終值', hl: true },
      { k: 'pgv', label: '|v|_{max} (m/s)', hl: true }, { k: 'pgd', label: '|u|_{max}' }, { k: 'pred', label: '預測 εt / ½εt²' }
    ]);
    const note = SD.div(body, 'w-note');
    const notes = [
      '① 偏差 ε 被積分成直線（速度）與拋物線（位移）。Z 方向 160 s 後「移動」了 −142.6 m。',
      '② 拿掉平均後終值小很多，但偏差在運轉中會變（X 方向尤其明顯），殘餘的低頻仍讓位移漂移數公尺。',
      '③ 帶通拿掉 0.1 Hz 以下的成分，速度已經合理；但每次積分都把 0.1 Hz 附近放大，位移仍有緩慢起伏。',
      '④ 每次積分後再高通一次：速度峰值約 0.004–0.013 m/s、位移約 1–2 mm。剩下週期 7–8 s 的 ±1 mm 起伏是雜訊被二次積分放大，洗衣機真實位移約 0.1 mm。'
    ];
    function draw() {
      const A = PHONE.get(ctl.ax), lv = +ctl.lv, L = A.L[lv], t = PHONE.t, T = t[t.length - 1];
      const big = Math.abs(L.d[L.d.length - 1]) > 0.5 || DSP.absmax(L.d).v > 0.05;
      const us = big ? 1 : 1000, uu = big ? 'm' : 'mm';
      const dd = L.d.map(x => x * us);
      const ma = DSP.absmax(L.a).v, mv = DSP.absmax(L.v).v, md = DSP.absmax(dd).v;
      pa.setLim(null, [-ma * 1.1, ma * 1.1]); pv.setLim(null, [-mv * 1.15, mv * 1.15]); pd.setLim(null, [-md * 1.15, md * 1.15]);
      pd.o.ylabel = 'u (' + uu + ')'; pd.drawAxes();
      pa.clear(); pv.clear(); pd.clear();
      let [X, Y] = thin(t, L.a); pa.line(X, Y, { c: 'ink', w: 0.8 });
      pa.hline(lv === 0 ? A.eps : 0, { c: 'red', dash: '6 4', w: 1.4 });
      [X, Y] = thin(t, L.v); pv.line(X, Y, { c: 'blue', w: 1.3 }); pv.hline(0, { dash: null, w: 0.6 });
      [X, Y] = thin(t, dd); pd.line(X, Y, { c: 'blue', w: 1.3 }); pd.hline(0, { dash: null, w: 0.6 });
      if (lv === 0) { pv.fn(tt => A.eps * tt, { c: 'gray', dash: '6 4', w: 1.2 }); pd.fn(tt => 0.5 * A.eps * tt * tt * us, { c: 'gray', dash: '6 4', w: 1.2 }); }
      const ue = L.d[L.d.length - 1] * us;
      st({ eps: f(A.eps, 4), vend: f(L.v[L.v.length - 1], 4), uend: f(ue, big ? 2 : 3) + ' ' + uu, pgv: f(mv, 4), pgd: f(md, big ? 2 : 3) + ' ' + uu,
        pred: f(A.eps * T, 3) + ' m/s / ' + f(0.5 * A.eps * T * T, 1) + ' m' });
      note.textContent = notes[lv];
    }
    draw();
  })();

  /* ================= W11 CSR ================= */
  (() => {
    const body = SD.widget('w-csr'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'amax', label: '地表 a_max (g)', min: 0.05, max: 0.6, step: 0.01, value: 0.25, fmt: v => v.toFixed(2) },
      { k: 'gwt', label: '地下水位深度 (m)', min: 0, max: 15, step: 0.5, value: 2, fmt: v => v.toFixed(1) },
      { k: 'g1', label: '水位以上 γ (kN/m³)', min: 15, max: 21, step: 0.1, value: 18, fmt: v => v.toFixed(1) },
      { k: 'g2', label: '水位以下 γ_sat (kN/m³)', min: 17, max: 22, step: 0.1, value: 19.5, fmt: v => v.toFixed(1) },
      { k: 'z', label: '深度 z (m)', min: 0.5, max: 20, step: 0.5, value: 6, fmt: v => v.toFixed(1) }
    ].map(d => Object.assign(d, { label: SD.rich(d.label) })), draw);
    SD.buttons(body, '範例', [
      { label: '講義例：z = 6 m、水位 2 m', on: () => ctl.set({ amax: 0.25, gwt: 2, g1: 18, g2: 19.5, z: 6 }) },
      { label: '測驗：水位很深', on: () => ctl.set({ gwt: 15, z: 6 }) }
    ]);
    const row = SD.div(body, 'plot-row');
    const p = new Plot(row, { w: 420, h: 340, xlim: [0, 0.4], ylim: [-20, 0], xlabel: 'CSR', ylabel: '深度 z (m)', yfmt: v => String(Math.abs(v)) });
    const p2 = new Plot(row, { w: 420, h: 340, xlim: [0, 400], ylim: [-20, 0], xlabel: '應力 (kPa)', ylabel: '深度 z (m)', yfmt: v => String(Math.abs(v)) });
    p.legend([{ label: 'CSR(z)', c: 'blue', w: 2 }, { label: '地下水位', c: 'gray', dash: true }, { label: '目前深度', c: 'blue', line: false, marker: 'circle', fill: true }]);
    p2.legend([{ label: 'σ_v', c: 'ink', w: 1.6 }, { label: 'σ′_v', c: 'blue', w: 2 }, { label: 'u', c: 'red', dash: true }, { label: 'r_d × 100', c: 'gray', w: 1.4 }]);
    const st = SD.stats(body, [
      { k: 'sv', label: 'σ_v (kPa)' }, { k: 'u', label: 'u (kPa)' }, { k: 'se', label: 'σ′_v (kPa)' }, { k: 'rd', label: 'r_d' },
      { k: 'tau', label: 'τ_av (kPa)' }, { k: 'csr', label: 'CSR', hl: true }
    ]);
    const rd = z => z <= 9.15 ? 1 - 0.00765 * z : 1.174 - 0.0267 * z;
    function at(z) {
      const w = ctl.gwt, sv = Math.min(z, w) * ctl.g1 + Math.max(0, z - w) * ctl.g2, u = Math.max(0, z - w) * 9.81, se = sv - u;
      return { sv, u, se, rd: rd(z), tau: 0.65 * ctl.amax * sv * rd(z), csr: 0.65 * ctl.amax * sv / se * rd(z) };
    }
    function draw() {
      const zs = linspace(0.2, 20, 200), R = zs.map(at);
      p.setLim([0, Math.max(0.3, ...R.map(r => r.csr)) * 1.1], null);
      p.clear(); p2.clear();
      p.hline(-ctl.gwt, { dash: '6 4' }); p2.hline(-ctl.gwt, { dash: '6 4' });
      p.line(R.map(r => r.csr), zs.map(z => -z), { c: 'blue', w: 2 });
      p2.line(R.map(r => r.sv), zs.map(z => -z), { c: 'ink', w: 1.6 });
      p2.line(R.map(r => r.se), zs.map(z => -z), { c: 'blue', w: 2 });
      p2.line(R.map(r => r.u), zs.map(z => -z), { c: 'red', dash: '6 4', w: 1.4 });
      p2.line(R.map(r => r.rd * 100), zs.map(z => -z), { c: 'gray', w: 1.4 });
      const c = at(ctl.z);
      p.markers([c.csr], [-ctl.z], { c: 'blue', fill: true, over: true });
      st({ sv: f(c.sv, 1), u: f(c.u, 1), se: f(c.se, 1), rd: f(c.rd, 3), tau: f(c.tau, 1), csr: f(c.csr, 3) });
    }
    draw();
  })();

  /* ================= W12 應變範圍 ================= */
  (() => {
    const body = SD.widget('w-strain'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'lg', label: '剪應變 log₁₀γ', min: -6, max: -1, step: 0.05, value: -4, fmt: v => '10^' + v.toFixed(2) }
    ], draw);
    SD.buttons(body, '典型情況', [
      { label: '震測 10⁻⁶', on: () => ctl.set({ lg: -5.7 }) },
      { label: '機器基礎 10⁻⁵', on: () => ctl.set({ lg: -5 }) },
      { label: '中度地震 10⁻⁴', on: () => ctl.set({ lg: -3.7 }) },
      { label: '液化 10⁻²', on: () => ctl.set({ lg: -1.7 }) }
    ]);
    const W = 760, H = 270, x0 = 120, x1 = 740;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot-svg' }); body.appendChild(svg);
    const g = svgEl('g', null, svg);
    const X = lg => x0 + (lg + 6) / 5 * (x1 - x0);
    const rows = [
      ['現象', [[-6, -4, '波傳、振動'], [-4, -3, '裂縫、差異沉陷'], [-2, -1, '滑動、壓密、液化']]],
      ['力學行為', [[-6, -5, '彈性'], [-5, -2, '彈塑性'], [-2, -1, '破壞']]],
      ['載重效應', [[-4, -1, '反覆次數、載重速率的影響']]],
      ['參數', [[-6, -3, 'G、ν、阻尼比 ξ'], [-2, -1, 'φ、c']]],
      ['現地試驗', [[-6, -5, '震測'], [-5, -4, '現地振動'], [-4, -2, '反覆載重試驗']]],
      ['室內試驗', [[-6, -5, '波傳'], [-5.5, -4, '共振柱'], [-4, -1, '反覆載重（三軸、單剪）']]]
    ];
    const p = new Plot(body, { w: 640, h: 270, xlim: [1e-6, 1e-1], ylim: [0, 1.05], logx: true, xlabel: '剪應變 γ', ylabel: 'G/G_{max}、ξ/ξ_{max}' });
    p.legend([{ label: 'G/G_{max}（示意）', c: 'blue', w: 2 }, { label: '阻尼比 ξ/ξ_{max}（示意）', c: 'red', w: 2 }]);
    const st = SD.stats(body, [{ k: 'g', label: 'γ' }, { k: 'beh', label: '力學行為', hl: true }, { k: 'par', label: '該用的參數' }, { k: 'test', label: '可用的室內試驗' }, { k: 'G', label: 'G/G_{max}（示意）' }]);
    SD.div(body, 'w-note', '上表依原講義 Ch0 p.5（Ishihara）的分區重繪，界線是量級而非精確值。下圖是 Hardin–Drnevich 雙曲線 G/G_max = 1/(1 + γ/γ_r)、γ_r = 0.05 % 的示意，單元 8 會用真實試驗資料取代它。');
    function draw() {
      clearG(g);
      for (let k = -6; k <= -1; k++) {
        svgEl('line', { x1: X(k), x2: X(k), y1: 16, y2: H - 22, class: 'guide' }, g);
        svgText(g, X(k), H - 6, '10^{' + k + '}', 'txt-muted', 'middle');
      }
      const lg = ctl.lg;
      rows.forEach(([name, segs], r) => {
        const y = 24 + r * 38;
        svgText(g, 8, y + 17, name, 'txt');
        segs.forEach(([a, b, s]) => {
          const on = lg >= a && lg <= b;
          svgEl('rect', { x: X(a) + 2, y, width: X(b) - X(a) - 4, height: 26, rx: 6, class: on ? 'f-blue' : 'f-band', opacity: on ? 0.9 : 1 }, g);
          const t = svgText(g, (X(a) + X(b)) / 2, y + 17, s, on ? 'txt' : 'txt', 'middle');
          if (on) t.setAttribute('style', 'fill:#fff;font-weight:600');
          t.setAttribute('font-size', '12');
        });
      });
      svgEl('line', { x1: X(lg), x2: X(lg), y1: 12, y2: H - 22, class: 's-red', 'stroke-width': 2 }, g);
      const gam = Math.pow(10, lg), Gr = 1 / (1 + gam / 5e-4);
      p.clear();
      p.fn(x => 1 / (1 + x / 5e-4), { c: 'blue', w: 2 });
      p.fn(x => 1 - 1 / (1 + x / 5e-4), { c: 'red', w: 2 });
      p.vline(gam, { c: 'red', dash: '4 3' });
      const pick = i => (rows[i][1].find(([a, b]) => lg >= a && lg <= b) || [0, 0, '—'])[2];
      st({ g: gam.toExponential(1), beh: pick(1), par: pick(3), test: pick(5), G: f(Gr, 3) });
    }
    draw();
  })();

  /* ================= W13 載重時間與循環數 ================= */
  (() => {
    const body = SD.widget('w-rate'); if (!body) return;
    const cats = [
      { k: 'blast', n: '衝擊、爆破', t: [1e-3, 1e-2], N: [1, 3], T: 0.005 },
      { k: 'eq', n: '地震', t: [0.05, 5], N: [10, 100], T: 0.5 },
      { k: 'pile', n: '打樁', t: [0.005, 0.05], N: [300, 3000], T: 0.02 },
      { k: 'mach', n: '機器基礎', t: [0.01, 0.1], N: [1e4, 1e8], T: 0.03 },
      { k: 'traf', n: '交通', t: [0.1, 1], N: [1e4, 1e7], T: 0.4 },
      { k: 'wave', n: '波浪', t: [5, 20], N: [1e3, 1e6], T: 10 },
      { k: 'park', n: '停車', t: [1e3, 1e5], N: [100, 1e4], T: 1e4 },
      { k: 'cons', n: '壓密、潛變', t: [1e6, 1e8], N: [1, 2], T: 1e7 }
    ];
    const ctl = SD.controls(body, [
      { k: 'cat', label: '載重', type: 'select', value: 'traf', options: cats.map(c => [c.k, c.n]) },
      { k: 'Tn', label: '系統自然週期 T_n (s)', min: 0.02, max: 3, log: true, snap: 2, value: 0.5, fmt: v => v.toFixed(2) }
    ].map(d => Object.assign(d, { label: SD.rich(d.label) })), draw);
    const p = new Plot(body, { w: 680, h: 360, xlim: [1e-3, 1e8], ylim: [1, 1e8], logx: true, logy: true, xlabel: '一次載重的時間（週期）(s)', ylabel: '循環次數 N' });
    p.legend([{ label: '動態區：載重時間 < 4T_n（β > 0.25）', c: 'gray', w: 8, op: .4 }, { label: '系統自然週期 T_n', c: 'red', dash: true }]);
    const st = SD.stats(body, [{ k: 'T', label: '代表載重週期 (s)' }, { k: 'b', label: 'β = T_n / T_{load}' }, { k: 'D', label: '放大 D（ξ = 5 %）', hl: true }, { k: 'j', label: '判斷', hl: true }]);
    SD.div(body, 'w-note', '橫軸是時間、縱軸是次數：左半部（載重快）才需要動力分析，右半部即使反覆很多次（停車、波浪）仍是靜態問題，只可能有疲勞效應。D 用單元 4 的動態放大因子，β 很小時 D ≈ 1，就是「準靜態」。');
    function draw() {
      p.clear();
      p.band(1e-3, 4 * ctl.Tn);
      p.vline(ctl.Tn, { c: 'red', dash: '6 4', w: 1.4 });
      cats.forEach(c => {
        const on = c.k === ctl.cat;
        const xa = p.sx(c.t[0]), xb = p.sx(c.t[1]), ya = p.sy(c.N[1]), yb = p.sy(c.N[0]);
        svgEl('rect', { x: xa, y: ya, width: Math.max(6, xb - xa), height: Math.max(8, yb - ya), rx: 3, class: (on ? 's-blue f-blue' : 's-ink f-surface'), 'stroke-width': 1.2, opacity: on ? 0.85 : 1 }, p.overG);
        // 方框夠大就把標籤放在框內，否則放在框上方（頂端太高才放框下方）
        const bw = Math.max(6, xb - xa), bh = Math.max(8, yb - ya), tw = 12.5 * c.n.length;
        const inside = bh >= 26 && bw >= tw + 6;
        let tx = (xa + xb) / 2, anchor = 'middle';
        if (tx - tw / 2 < p.x0 + 4) { tx = p.x0 + 4; anchor = 'start'; }   // 貼近左框線時改靠左對齊
        const ty = inside ? (ya + yb) / 2 + 4 : (ya - 6 < 22 ? yb + 15 : ya - 6);
        const t = svgEl('text', { x: tx, y: ty, 'text-anchor': anchor, class: on ? (inside ? 'txt' : 'txt-blue') : 'txt', 'font-size': 12, style: on && inside ? 'fill:#fff' : null }, p.overG);
        t.textContent = c.n;
      });
      const c = cats.find(x => x.k === ctl.cat), b = ctl.Tn / c.T, D = S.D(b, 0.05);
      st({ T: c.T >= 100 ? c.T.toExponential(0) : String(c.T), b: b < 1e-3 ? b.toExponential(1) : f(b, 3), D: f(D, 3),
        j: b < 0.1 ? '準靜態（慣性可忽略）' : b < 0.25 ? '接近靜態' : b <= 3 ? '動態：需考慮慣性與共振' : '衝擊型：慣性主導' }, { j: b >= 0.25 });
    }
    draw();
  })();

  /* ================= W14 HW1 Part A FAS ================= */
  (() => {
    const body = SD.widget('w-hw1a'); if (!body || !PHONE) return;
    const ctl = SD.controls(body, [
      { k: 'ax', label: '方向', type: 'select', value: 'y', options: [['x', 'X 左右'], ['y', 'Y 前後'], ['z', 'Z 上下']] },
      { k: 'seg', label: '時段', type: 'select', value: 'all', options: [['all', '全段 0–161.7 s'], ['0,13', '0–13 s 正反轉'], ['13,91', '13–91 s 注水'], ['91,131', '91–131 s 振動群'], ['136,162', '136–162 s 洗滌']] },
      { k: 'dec', label: '取樣率', type: 'select', value: '1', options: [['1', '原始 100.34 Hz'], ['2', '每 2 點取 1（50.17 Hz，無抗混疊）'], ['3', '每 3 點取 1（33.45 Hz，無抗混疊）']] }
    ], draw);
    const p = new Plot(body, { h: 320, xlim: [0, 50.2], ylim: [0, 1.3], xlabel: '頻率 f (Hz)', ylabel: 'FAS (m/s)' });
    p.legend([{ label: '帶通後加速度的 FAS', c: 'ink', w: 1 }, { label: '卓越頻率 f_p', c: 'red', line: false, marker: 'circle', fill: true }, { label: 'Nyquist f_N', c: 'gray', dash: true }]);
    const st = SD.stats(body, [
      { k: 'fs', label: 'f_s (Hz)' }, { k: 'fn', label: 'f_N (Hz)' }, { k: 'df', label: 'Δf = 1/(NΔt) (Hz)' },
      { k: 'fp', label: 'f_p (Hz)', hl: true }, { k: 'tp', label: 'T_p (s)' }, { k: 'A', label: 'FAS(f_p) (m/s)' }
    ]);
    const note = SD.div(body, 'w-note');
    function draw() {
      const A = PHONE.get(ctl.ax), dt0 = PHONE.dt, k = +ctl.dec;
      let i0 = 0, i1 = PHONE.n;
      if (ctl.seg !== 'all') { const [a, b] = ctl.seg.split(',').map(Number); i0 = Math.round(a / dt0); i1 = Math.min(PHONE.n, Math.round(b / dt0)); }
      const x = []; for (let i = i0; i < i1; i += k) x.push(A.bp[i]);
      const dt = dt0 * k, fs = 1 / dt, F = DSP.fas(x, dt, 4 * x.length);
      let ip = 1;
      for (let j = 1; j < F.f.length; j++) if (F.f[j] >= 0.1 && F.f[j] <= Math.min(30, fs / 2) && F.A[j] > F.A[ip]) ip = j;
      // 在原始 DFT 頻率格 k/(NΔt) 上精修
      const N = x.length, df = 1 / (N * dt);
      let best = { fq: F.f[ip], A: 0 };
      for (let m = Math.floor(F.f[ip] / df) - 3; m <= Math.ceil(F.f[ip] / df) + 3; m++) { const a = DSP.dftAt(x, dt, m * df); if (a > best.A) best = { fq: m * df, A: a }; }
      p.setLim([0, fs / 2], [0, best.A * 1.15]);
      p.clear();
      const [X, Y] = thin(F.f, F.A, 1500);
      p.line(X, Y, { c: 'ink', w: 0.9 });
      p.vline(fs / 2, { c: 'gray', dash: '6 4' });
      p.markers([best.fq], [best.A], { c: 'red', fill: true, over: true });
      st({ fs: f(fs, 2), fn: f(fs / 2, 2), df: f(df, 4), fp: f(best.fq, 2), tp: f(1 / best.fq, 3), A: f(best.A, 3) });
      note.textContent = k === 1
        ? '主峰在約 14 Hz（洗衣槽轉動），約 28 Hz 有倍頻。紀錄很長、Δf 很小，所以峰附近有許多細峰。改選時段：注水段幾乎只剩雜訊。'
        : `取樣率降為 ${fs.toFixed(2)} Hz、Nyquist ${(fs / 2).toFixed(2)} Hz：28 Hz 的倍頻超過 Nyquist，混疊到 ${Math.abs(28 - fs * Math.round(28 / fs)).toFixed(1)} Hz 附近出現一根假峰，單看頻譜無法分辨真假。`;
    }
    draw();
  })();

  /* ================= W15 RotD100 ================= */
  (() => {
    const body = SD.widget('w-rotd'); if (!body || !STA.length) return;
    const ctl = SD.controls(body, [
      { k: 'sta', label: '測站', type: 'select', value: STA[0], options: STA.map(k => [k, staLabel(k)]) },
      { k: 'th', label: '旋轉角 θ (°)', min: 0, max: 179, step: 1, value: 0, fmt: v => v.toFixed(0) }
    ], draw);
    SD.buttons(body, null, [{ label: '轉到 RotD100 的角度', primary: true, on: () => ctl.set({ th: cur.imax }) }]);
    const pt = new Plot(body, { h: 190, xlim: [0, 90], ylim: [-0.3, 0.3], xlabel: '時間 t (s)', ylabel: 'a_θ (g)' });
    pt.legend([{ label: '旋轉後的分量 a_θ = a_1cosθ + a_2sinθ', c: 'ink', w: 1 }, { label: '此方向的峰值', c: 'red', line: false, marker: 'circle', fill: true }]);
    const row = SD.div(body, 'plot-row');
    const pp = new Plot(row, { w: 420, h: 300, xlim: [0, 180], ylim: [0, 0.3], xticks: [0, 30, 60, 90, 120, 150, 180], xlabel: '旋轉角 θ (°)', ylabel: 'PGA(θ) (g)' });
    const po = new Plot(row, { w: 360, h: 330, xlim: [-0.3, 0.3], ylim: [-0.3, 0.3], xlabel: 'a_1 (g)', ylabel: 'a_2 (g)' });
    pp.legend([{ label: 'PGA(θ)', c: 'blue', w: 2 }, { label: 'RotD100', c: 'ink', dash: true }, { label: 'RotD50（中位數）', c: 'gray', dash: true }]);
    po.legend([{ label: '質點運動軌跡（實測）', c: 'ink', w: 0.8 }, { label: '目前方向 θ', c: 'red', w: 2 }]);
    const st = SD.stats(body, [
      { k: 'p', label: 'PGA(θ) (g)', hl: true }, { k: 'r100', label: 'PGA_{RotD100} (g)', hl: true }, { k: 'ang', label: 'RotD100 角度 (°)' },
      { k: 'r50', label: 'PGA_{RotD50} (g)' }, { k: 'h', label: '|a_1|_{max}、|a_2|_{max} (g)' }, { k: 'res', label: '合成向量最大值 (g)' }
    ]);
    let cur = { imax: 0 };
    const cache = {};
    function prep(sta) {
      if (cache[sta]) return cache[sta];
      const [c1, c2] = horiz(sta), A1 = getComp(sta, c1).a, A2 = getComp(sta, c2).a, n = Math.min(A1.length, A2.length);
      const pg = [];
      for (let d = 0; d < 180; d++) { const c = Math.cos(rad(d)), s = Math.sin(rad(d)); let m = 0; for (let i = 0; i < n; i++) { const v = Math.abs(A1[i] * c + A2[i] * s); if (v > m) m = v; } pg.push(m / G); }
      let res = 0; for (let i = 0; i < n; i++) res = Math.max(res, Math.hypot(A1[i], A2[i]));
      let imax = 0; pg.forEach((v, i) => { if (v > pg[imax]) imax = i; });
      const srt = [...pg].sort((a, b) => a - b), r50 = (srt[89] + srt[90]) / 2;
      return (cache[sta] = { c1, c2, A1, A2, n, pg, res: res / G, imax, r50, t: tvec(n, REC[sta].dt) });
    }
    function draw() {
      const D = cur = prep(ctl.sta), th = rad(ctl.th), c = Math.cos(th), s = Math.sin(th);
      const at = new Float64Array(D.n); for (let i = 0; i < D.n; i++) at[i] = (D.A1[i] * c + D.A2[i] * s) / G;
      const mx = D.pg[D.imax] * 1.15, pk = DSP.absmax(at);
      pt.setLim([0, Math.ceil(D.t[D.n - 1] / 10) * 10], [-mx, mx]); pt.clear();
      let [X, Y] = thin(D.t, at); pt.line(X, Y, { c: 'ink', w: 0.8 });
      pt.markers([D.t[pk.i]], [at[pk.i]], { c: 'red', fill: true, over: true });
      pp.setLim(null, [0, mx]); pp.clear();
      pp.line(D.pg.map((_, i) => i), D.pg, { c: 'blue', w: 2 });
      pp.hline(D.pg[D.imax], { c: 'ink', dash: '6 4', w: 1.2 }); pp.hline(D.r50, { c: 'gray', dash: '6 4', w: 1.2 });
      pp.markers([ctl.th], [D.pg[ctl.th]], { c: 'red', fill: true, over: true });
      po.setLim([-mx, mx], [-mx, mx]); po.clear();
      const step = Math.max(1, Math.floor(D.n / 5000)), ox = [], oy = [];
      for (let i = 0; i < D.n; i += step) { ox.push(D.A1[i] / G); oy.push(D.A2[i] / G); }
      po.line(ox, oy, { c: 'ink', w: 0.6, op: 0.8 });
      po.line([-mx * c, mx * c], [-mx * s, mx * s], { c: 'red', w: 2 });
      const m1 = DSP.absmax(D.A1).v / G, m2 = DSP.absmax(D.A2).v / G;
      st({ p: f(D.pg[ctl.th], 4), r100: f(D.pg[D.imax], 4), ang: String(D.imax), r50: f(D.r50, 4), h: f(m1, 3) + '、' + f(m2, 3), res: f(D.res, 4) });
    }
    draw();
  })();

  /* ================= W16 HW1 Part B 散佈圖 ================= */
  (() => {
    const body = SD.widget('w-hw1b'); const H = window.CH02_HW1; if (!body || !H) return;
    const K = H.keys, col = k => H.rows.map(r => r[K.indexOf(k)]);
    const lab = { rrup: 'R_rup (km)', vs30: 'V_s30 (m/s)', pga: 'PGA (g)', pgv: 'PGV (cm/s)', pgd: 'PGD (cm)', cav: 'CAV (m/s)', ai: 'AI (m/s)', tm: 'T_m (s)', ds: 'D_s (s)' };
    const ctl = SD.controls(body, [
      { k: 'cp', label: '分量', type: 'select', value: 'h', options: [['h', '水平 RotD100'], ['v', '垂直 V']] },
      { k: 'x', label: 'x 軸', type: 'select', value: 'pga', options: [['pga', 'PGA'], ['rrup', 'R_rup'], ['vs30', 'V_s30']] },
      { k: 'y', label: 'y 軸', type: 'select', value: 'ai', options: ['pga', 'pgv', 'pgd', 'cav', 'ai', 'tm', 'ds'].map(k => [k, k.toUpperCase().replace('TM', 'T_m').replace('DS', 'D_s').replace('AI', 'AI')]) }
    ], draw);
    SD.buttons(body, '題目', [
      { label: 'Q5(d) PGA–AI', on: () => ctl.set({ cp: 'h', x: 'pga', y: 'ai' }) },
      { label: 'Q5(e) PGA–T_m', on: () => ctl.set({ cp: 'h', x: 'pga', y: 'tm' }) },
      { label: 'Q6(a) R_rup–PGA', on: () => ctl.set({ cp: 'h', x: 'rrup', y: 'pga' }) },
      { label: 'Q7(a) 垂直 R_rup–PGA', on: () => ctl.set({ cp: 'v', x: 'rrup', y: 'pga' }) }
    ].map(b => Object.assign(b, { label: b.label.replace(/_/g, '') })));
    const p = new Plot(body, { w: 640, h: 360, xlim: [0.01, 1], ylim: [0.01, 1], logx: true, logy: true, xlabel: '', ylabel: '' });
    p.legend([{ label: '30 筆紀錄（實測）', c: 'ink', line: false, marker: 'circle' }, { label: '迴歸線', c: 'blue', w: 2 }, { label: '本頁嵌入的 4 站', c: 'red', line: false, marker: 'circle', fill: true }]);
    p.readout(() => '');
    const st = SD.stats(body, [{ k: 'n', label: 'N' }, { k: 'r', label: '相關係數 r', hl: true }, { k: 'b', label: '斜率 b', hl: true }, { k: 'a', label: '截距 a' }, { k: 'm', label: '模型' }]);
    const note = SD.div(body, 'w-note');
    const nice = (lo, hi) => [Math.pow(10, Math.floor(Math.log10(lo))), Math.pow(10, Math.ceil(Math.log10(hi)))];
    function draw() {
      const v = ctl.cp === 'v', yk = ctl.y + (v ? 'V' : ''), xk = ctl.x === 'pga' ? 'pga' + (v ? 'V' : '') : ctl.x;
      const xs = col(xk), ys = col(yk), sta = col('sta');
      const linY = ctl.y === 'tm' || ctl.y === 'ds';
      const lx = xs.map(Math.log10), ly = linY ? ys : ys.map(Math.log10);
      const n = xs.length, mx = lx.reduce((a, b) => a + b) / n, my = ly.reduce((a, b) => a + b) / n;
      let sxy = 0, sxx = 0, syy = 0; lx.forEach((x, i) => { sxy += (x - mx) * (ly[i] - my); sxx += (x - mx) ** 2; syy += (ly[i] - my) ** 2; });
      const b = sxy / sxx, a = my - b * mx, r = sxy / Math.sqrt(sxx * syy);
      const xl = nice(Math.min(...xs), Math.max(...xs));
      if (linY) { const lo = Math.min(...ys), hi = Math.max(...ys), pad = (hi - lo) * 0.15; p.o.logy = false; p.setLim(xl, [Math.max(0, lo - pad), hi + pad]); }
      else { p.o.logy = true; p.setLim(xl, nice(Math.min(...ys), Math.max(...ys))); }
      p.o.xlabel = (ctl.x === 'pga' ? (v ? 'PGA_V (g)' : 'PGA_{RotD100} (g)') : lab[ctl.x]);
      p.o.ylabel = (lab[ctl.y] || ctl.y).replace(' (', (v ? '_V' : '_{RotD100}') + ' (').replace('_s_', '_s').replace('_m_', '_m');
      p.drawAxes();
      p.clear();
      p.fn(x => linY ? a + b * Math.log10(x) : Math.pow(10, a + b * Math.log10(x)), { c: 'blue', w: 2 });
      p.markers(xs, ys, { c: 'ink' });
      const mine = sta.map((s, i) => STA.includes(s) ? i : -1).filter(i => i >= 0);
      p.markers(mine.map(i => xs[i]), mine.map(i => ys[i]), { c: 'red', fill: true, over: true, r: 4 });
      mine.forEach(i => p.text(xs[i], ys[i], sta[i], { dx: 7, dy: -6, cls: 'txt-red' }));
      p.readout(x => {
        let j = 0; xs.forEach((xx, i) => { if (Math.abs(Math.log10(xx / x)) < Math.abs(Math.log10(xs[j] / x))) j = i; });
        return `${sta[j]}：x = ${xs[j]}　y = ${ys[j]}`;
      });
      st({ n: String(n), r: f(r, 3), b: f(b, 3), a: f(a, 3), m: linY ? 'y = a + b log x' : 'log y = a + b log x' });
      note.innerHTML = SD.rich('表中數值直接讀自 B1_process.xlsx（Q9_Summary），迴歸由網頁重算，與報告表 B-9～B-11 的 r、b 相同（到小數第 2 位）。也可以選 x = V_s30 看場址的影響：這不是作業要求，但能看出 PGV 隨場址變軟而變大的趨勢。');
    }
    draw();
  })();
});
