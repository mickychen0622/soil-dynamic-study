/* =========================================================
   Ch8 土壤液化：互動元件（單元 12–14）
   - 三套 SPT 簡化法：Youd et al. (2001)、HBF (2012)、Boulanger & Idriss (2014)
   - 液化後體積應變與最大剪應變：Idriss & Boulanger (2008) 對 Ishihara & Yoshimine (1992) 的近似
   ========================================================= */
(function () {
  const PA = 100;                    // kPa（講義取 pa ≈ 100 kPa）
  const GW = 9.81;                   // 水單位重 kN/m³

  /* ---------- 應力 ---------- */
  function stresses(z, gwt, g1, g2) {
    const sv = z <= gwt ? g1 * z : g1 * gwt + g2 * (z - gwt);
    const u = Math.max(0, z - gwt) * GW;
    return { sv, u, se: sv - u };
  }
  /* 桿長修正（桿長 = 深度 + 地面以上約 1 m） */
  function cR(rod) { return rod < 3 ? 0.75 : rod < 4 ? 0.80 : rod < 6 ? 0.85 : rod < 10 ? 0.95 : 1.0; }

  /* ---------- 三套方法 ---------- */
  const M = {
    youd: {
      name: 'Youd et al. (2001)',
      rd: z => (1 - 0.4113 * Math.sqrt(z) + 0.04052 * z + 0.001753 * Math.pow(z, 1.5)) /
        (1 - 0.4177 * Math.sqrt(z) + 0.05729 * z - 0.006205 * Math.pow(z, 1.5) + 0.001210 * z * z),
      run(p) {
        const CN = Math.min(1.7, Math.sqrt(PA / p.se));
        const N1 = CN * p.N60;
        const FC = p.FC;
        const a = FC <= 5 ? 0 : FC >= 35 ? 5 : Math.exp(1.76 - 190 / (FC * FC));
        const b = FC <= 5 ? 1 : FC >= 35 ? 1.2 : 0.99 + Math.pow(FC, 1.5) / 1000;
        const Ncs = a + b * N1;
        const crr75 = Ncs >= 30 ? Infinity : 1 / (34 - Ncs) + Ncs / 135 + 50 / Math.pow(10 * Ncs + 45, 2) - 1 / 200;
        const msf = Math.pow(10, 2.24) / Math.pow(p.Mw, 2.56);
        const Ks = p.se > PA ? Math.pow(p.se / PA, 0.7 - 1) : 1;
        return { rd: this.rd(p.z), CN, N1, Ncs, crr75, msf, Ks, note: `α = ${a.toFixed(2)}、β = ${b.toFixed(3)}` };
      }
    },
    hbf: {
      name: 'HBF (2012)',
      rd: z => z <= 10 ? 1 - 0.01 * z : 1.2 - 0.03 * z,
      run(p) {
        const CN = 1 / Math.sqrt(p.se / 98.1);           // σ′ 以 kg/cm² 計
        const N1 = CN * p.N60;
        const Ks = p.FC <= 10 ? 1 : 1 + 0.07 * Math.sqrt(p.FC - 10);
        const Ncs = Ks * N1;
        const crr75 = Ncs >= 39 ? Infinity : 0.08 + 0.0035 * Ncs / (1 - Ncs / 39);
        const msf = Math.pow(p.Mw / 7.5, -3.3);
        return { rd: this.rd(p.z), CN, N1, Ncs, crr75, msf, Ks: 1, note: `細粒料 K_s = ${Ks.toFixed(3)}` };
      }
    },
    bi: {
      name: 'Boulanger & Idriss (2014)',
      rd(z, Mw) {
        const a = -1.012 - 1.126 * Math.sin(z / 11.73 + 5.133);
        const b = 0.106 + 0.118 * Math.sin(z / 11.28 + 5.142);
        return Math.exp(a + b * Mw);
      },
      run(p) {
        const FC = p.FC;
        const dN = Math.exp(1.63 + 9.7 / (FC + 0.01) - Math.pow(15.7 / (FC + 0.01), 2));
        let Ncs = p.N60 + dN, CN = 1, N1 = p.N60;
        for (let i = 0; i < 50; i++) {                  // C_N 與 (N1)60cs 疊代
          const m = 0.784 - 0.0768 * Math.sqrt(Math.min(Math.max(Ncs, 0), 46));
          CN = Math.min(1.7, Math.pow(PA / p.se, m));
          N1 = CN * p.N60;
          const nx = N1 + dN;
          if (Math.abs(nx - Ncs) < 1e-6) { Ncs = nx; break; }
          Ncs = nx;
        }
        const n = Math.min(Ncs, 46);
        const crr75 = Math.exp(n / 14.1 + Math.pow(n / 126, 2) - Math.pow(n / 23.6, 3) + Math.pow(n / 25.4, 4) - 2.8);
        const mmax = Math.min(2.2, 1.09 + Math.pow(n / 31.5, 2));
        const msf = 1 + (mmax - 1) * (8.64 * Math.exp(-p.Mw / 4) - 1.325);
        const Cs = Math.min(0.3, 1 / (18.9 - 2.55 * Math.sqrt(n)));
        const Ks = Math.min(1.1, 1 - Cs * Math.log(p.se / PA));
        return { rd: this.rd(p.z, p.Mw), CN, N1, Ncs, crr75, msf, Ks, note: `Δ(N1)60 = ${dN.toFixed(2)}（疊代）` };
      }
    }
  };

  /** 單一深度的完整評估 */
  function evalPoint(p, method) {
    const s = stresses(p.z, p.gwt, p.g1, p.g2);
    const N60 = p.Nm * p.ER / 60 * cR(p.z + 1);
    const q = Object.assign({}, p, { se: s.se, N60 });
    const r = M[method].run(q);
    const csr = 0.65 * p.amax * s.sv / s.se * r.rd;
    const crr = r.crr75 * r.msf * r.Ks;
    return Object.assign(r, s, { N60, CR: cR(p.z + 1), csr, crr, fs: crr / csr });
  }

  /** 液化後體積應變與最大剪應變（I&B 2008，Ishihara & Yoshimine 圖的近似） */
  function postLiq(fs, Ncs) {
    const Dr = Math.min(1, Math.sqrt(Math.max(Ncs, 0) / 46));
    const glim = Math.max(0, 1.859 * Math.pow(1.1 - Dr, 3));
    const Fa = 0.032 + 4.7 * Dr - 6.0 * Dr * Dr;
    let g;
    if (fs >= 2) g = 0;
    else if (fs <= Fa) g = glim;
    else g = Math.min(glim, Math.max(0, 0.035 * (1 - Fa) * (2 - fs) / (fs - Fa)));
    const ev = 1.5 * Math.exp(-2.5 * Dr) * Math.min(0.08, g);
    return { Dr, gmax: g, ev };
  }

  window.LIQ = { M, stresses, evalPoint, postLiq, cR };
})();

document.addEventListener('sd:ready', () => {
  const f = SD.f;
  const L = window.LIQ;
  const fsTxt = v => !isFinite(v) || v > 2 ? '> 2' : v.toFixed(2);

  /* =========================================================
     1. #w-ru（12.3）：反覆三軸試驗 — CRR 曲線、N_L 與 r_u 累積
     ========================================================= */
  (() => {
    const body = SD.widget('w-ru'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'dr', label: '相對密度 D_r (%)', min: 30, max: 80, step: 1, value: 47, fmt: v => v.toFixed(0) + ' %' },
      { k: 'csr', label: '試驗的反覆應力比 CSR', min: 0.06, max: 0.45, step: 0.005, value: 0.13, fmt: v => v.toFixed(3) },
      { k: 'b', label: 'CRR 曲線指數 b', min: 0.15, max: 0.40, step: 0.01, value: 0.25, fmt: v => v.toFixed(2) }
    ], draw);
    SD.buttons(body, '情境', [
      { label: '鬆砂 D_r = 47 %', on: () => ctl.set({ dr: 47, csr: 0.13, b: 0.25 }) },
      { label: '密砂 D_r = 75 %', on: () => ctl.set({ dr: 75, csr: 0.30, b: 0.25 }) },
      { label: '同一 CSR，密砂', on: () => ctl.set({ dr: 75, csr: 0.13, b: 0.25 }) }
    ]);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pC = new Plot(grid, { h: 280, logx: true, xlim: [1, 1000], ylim: [0, 0.6], xlabel: '液化圈數 N_L（對數）', ylabel: 'CSR' });
    pC.legend([{ label: 'CRR 曲線 CRR = aN^{−b}', c: 'blue', w: 2 }, { label: '試驗的 CSR', c: 'red', dash: true }, { label: '15 圈（M = 7.5 基準）', c: 'gray', dash: '3 3' }]);
    const pR = new Plot(grid, { h: 280, xlim: [0, 30], ylim: [0, 1.05], xlabel: '圈數 N', ylabel: 'r_u = Δu/σ′_0（藍）、應變示意（紅）' });
    pR.legend([{ label: 'r_u 包絡（Seed et al. 1975）', c: 'blue', w: 2 }, { label: '單圈應變（示意，相對值）', c: 'red' }]);
    const st = SD.stats(body, [
      { k: 'n1', label: '(N_1)_{60} ≈ 46 D_r²' },
      { k: 'crr15', label: 'CRR（15 圈）' },
      { k: 'nl', label: '液化圈數 N_L', hl: true },
      { k: 'type', label: '破壞型態', hl: true }
    ]);

    function draw() {
      const Dr = ctl.dr / 100, b = ctl.b;
      const n1 = 46 * Dr * Dr;
      const yo = L.M.youd.run({ N60: n1, se: 100, FC: 0, Mw: 7.5, z: 5 });
      const crr15 = Math.min(isFinite(yo.crr75) ? yo.crr75 : 0.6, 0.6);
      const a = crr15 * Math.pow(15, b);
      const NL = Math.pow(ctl.csr / a, -1 / b);
      pC.clear();
      pC.fn(N => a * Math.pow(N, -b), { c: 'blue', w: 2 });
      pC.hline(ctl.csr, { c: 'red', dash: '5 4', w: 1.5 });
      pC.vline(15, { c: 'gray', dash: '3 3' });
      if (NL >= 1 && NL <= 1000) pC.markers([NL], [ctl.csr], { c: 'red', fill: true, over: true });
      const xmax = Math.max(10, Math.min(200, Math.ceil(NL * 1.3)));
      pR.setLim([0, xmax], [0, 1.05]);
      pR.clear();
      pR.hline(1, { c: 'gray', dash: '3 3' });
      const ru = N => N >= NL ? 1 : 0.5 + Math.asin(2 * Math.pow(N / NL, 1 / 0.7) - 1) / Math.PI;
      pR.fn(ru, { c: 'blue', w: 2, n: 400 });
      const dense = Dr >= 0.6;
      // 應變：r_u < 0.8 時很小；之後鬆砂失控、密砂每圈有限增加
      pR.fn(N => {
        const r = ru(N);
        if (r < 0.8) return 0.03 * r;
        if (N < NL) return 0.03 + 0.25 * Math.pow((r - 0.8) / 0.2, 3);
        return dense ? Math.min(0.8, 0.28 + 0.06 * (N - NL)) : Math.min(1.0, 0.28 + 1.2 * (N - NL));
      }, { c: 'red', w: 1.6, n: 400 });
      if (NL <= xmax) pR.vline(NL, { c: 'red', dash: '4 3' });
      st({
        n1: f(n1, 1), crr15: f(crr15, 3),
        nl: NL > 1e4 ? '> 10⁴（不液化）' : NL < 1 ? '< 1' : f(NL, 1),
        type: NL > 1e4 ? '—' : dense ? '反覆流動性（應變有限）' : '流動型（應變失控）'
      });
    }
    draw();
  })();

  /* =========================================================
     2. #w-neq（13.2）：不規則歷時 → N_eq
     ========================================================= */
  (() => {
    const body = SD.widget('w-neq'); if (!body || !window.CH6_MOTION) return;
    const R = window.CH6_MOTION;
    const tau = R.a.map(v => v * R.unit);
    const amax = Math.max(...tau.map(Math.abs));
    const tn = tau.map(v => v / amax);                        // τ/τmax（剛性土柱：τ ∝ a）
    const t = tn.map((_, i) => i * R.dt);
    // 以零交越切成半圈，每半圈取最大 |τ|
    const halves = [];
    let s = 0;
    for (let i = 1; i <= tn.length; i++) {
      if (i === tn.length || Math.sign(tn[i]) !== Math.sign(tn[i - 1]) && tn[i] !== 0) {
        let pk = 0, ip = s;
        for (let j = s; j < i; j++) if (Math.abs(tn[j]) > pk) { pk = Math.abs(tn[j]); ip = j; }
        if (pk > 0.02) halves.push({ amp: pk, t: ip * R.dt });
        s = i;
      }
    }
    const ctl = SD.controls(body, [
      { k: 'b', label: 'CRR 曲線指數 b', min: 0.15, max: 0.40, step: 0.01, value: 0.25, fmt: v => v.toFixed(2) + `（1/b = ${(1 / v).toFixed(1)}）` }
    ], draw);
    const grid = SD.div(body, 'plot-row plot-grid');
    const pT = new Plot(grid, { h: 260, xlim: [0, t[t.length - 1]], ylim: [-1.05, 1.05], xlabel: '時間 t (s)', ylabel: 'τ / τ_{max}' });
    pT.legend([{ label: `${R.name}（剛性土柱：τ ∝ a）`, c: 'ink', w: 1 }, { label: '±0.65τ_{max}', c: 'red', dash: true }]);
    const pN = new Plot(grid, { h: 260, xlim: [0, t[t.length - 1]], ylim: [0, 20], xlabel: '時間 t (s)', ylabel: '累積等值圈數 N_{eq}' });
    pN.legend([{ label: '累積 N_{eq}', c: 'blue', w: 2 }, { label: '只數圈數（不加權）/ 2', c: 'gray', dash: true }]);
    const st = SD.stats(body, [
      { k: 'h', label: '半圈數（振幅 > 2 %）' },
      { k: 'big', label: '振幅 ≥ 0.65τ_{max} 的半圈' },
      { k: 'neq', label: 'N_{eq}（0.65τ_{max}）', hl: true },
      { k: 'm', label: '對應規模（講義 p.33 表內插）', hl: true }
    ]);
    pT.line(t, tn, { c: 'ink', w: 0.8 });
    pT.hline(0.65, { c: 'red', dash: '5 4' }); pT.hline(-0.65, { c: 'red', dash: '5 4' });
    const tbl = [[2.5, 5.25], [5.5, 6.0], [10, 6.75], [15, 7.5], [26, 8.5]];
    function mOf(n) {
      if (n <= tbl[0][0]) return '≤ 5.25';
      for (let i = 1; i < tbl.length; i++) if (n <= tbl[i][0]) {
        const [n0, m0] = tbl[i - 1], [n1, m1] = tbl[i];
        return 'M ≈ ' + (m0 + (m1 - m0) * (n - n0) / (n1 - n0)).toFixed(2);
      }
      return '≥ 8.5';
    }
    function draw() {
      const b = ctl.b;
      let acc = 0, cnt = 0;
      const xs = [0], ys = [0], ys2 = [0];
      halves.forEach(h => {
        acc += 0.5 * Math.pow(h.amp / 0.65, 1 / b); cnt += 0.5;
        xs.push(h.t); ys.push(acc); ys2.push(cnt / 2);
      });
      xs.push(t[t.length - 1]); ys.push(acc); ys2.push(cnt / 2);
      pN.setLim([0, t[t.length - 1]], [0, Math.max(10, Math.ceil(Math.max(acc, cnt / 2) * 1.15))]);
      pN.clear();
      pN.line(xs, ys, { c: 'blue', w: 2 });
      pN.line(xs, ys2, { c: 'gray', dash: '5 4' });
      st({ h: halves.length, big: halves.filter(h => h.amp >= 0.65).length, neq: f(acc, 1), m: mOf(acc) });
    }
    draw();
  })();

  /* =========================================================
     3. #w-spt（13.4、14.1）：單一深度的 SPT 簡化法
     ========================================================= */
  (() => {
    const body = SD.widget('w-spt'); if (!body) return;
    const ctl = SD.controls(body, [
      { k: 'method', type: 'select', label: '方法', value: 'youd', options: [['youd', 'Youd et al. (2001)'], ['hbf', 'HBF (2012)'], ['bi', 'Boulanger & Idriss (2014)']] },
      { k: 'z', label: '深度 z (m)', min: 1, max: 20, step: 0.5, value: 6, fmt: v => v.toFixed(1) + ' m' },
      { k: 'gwt', label: '地下水位 (m)', min: 0, max: 6, step: 0.5, value: 1, fmt: v => v.toFixed(1) + ' m' },
      { k: 'amax', label: '地表 a_max (g)', min: 0.05, max: 0.6, step: 0.01, value: 0.24, fmt: v => v.toFixed(2) + ' g' },
      { k: 'Mw', label: '規模 M_w', min: 5.5, max: 8.5, step: 0.1, value: 7.0, fmt: v => v.toFixed(1) },
      { k: 'Nm', label: '實測 N_m', min: 1, max: 50, step: 1, value: 11, fmt: v => v.toFixed(0) },
      { k: 'FC', label: '細粒料含量 FC (%)', min: 0, max: 50, step: 1, value: 12, fmt: v => v.toFixed(0) + ' %' },
      { k: 'ER', label: '能量比 ER (%)', min: 40, max: 100, step: 5, value: 70, fmt: v => v.toFixed(0) + ' %' }
    ], draw);
    SD.buttons(body, '', [
      { label: '講義例題（p.50）', on: () => ctl.set({ method: 'youd', z: 6, gwt: 1, amax: 0.24, Mw: 7.0, Nm: 11, FC: 12, ER: 70 }) },
      { label: '同一點改用 HBF', on: () => ctl.set({ method: 'hbf', z: 6, gwt: 1, amax: 0.24, Mw: 7.0, Nm: 11, FC: 12, ER: 70 }) },
      { label: '同一點改用 B&I', on: () => ctl.set({ method: 'bi', z: 6, gwt: 1, amax: 0.24, Mw: 7.0, Nm: 11, FC: 12, ER: 70 }) }
    ]);
    const tw = SD.div(body, 'tablewrap');
    const st = SD.stats(body, [
      { k: 'csr', label: 'CSR（載重）' },
      { k: 'crr', label: 'CRR（抵抗）' },
      { k: 'fs', label: 'FS_L = CRR/CSR', hl: true },
      { k: 'need', label: 'FS_L = 1.2 所需 CRR_{7.5}' }
    ]);
    SD.div(body, 'mini-note', '單位重固定：地下水位以上 18.5、以下 19.0 kN/m³；C_B = C_S = 1.0；桿長取深度 + 1 m 查 C_R。HBF 的 C_N 依講義不設上限；Youd 的 K_σ 取 f = 0.7。');
    function draw() {
      const p = { z: ctl.z, gwt: ctl.gwt, g1: 18.5, g2: 19.0, amax: ctl.amax, Mw: ctl.Mw, Nm: ctl.Nm, FC: ctl.FC, ER: ctl.ER };
      const r = L.evalPoint(p, ctl.method);
      const rows = [
        ['σ_{v0}、u、σ′_{v0}', `${f(r.sv, 1)}、${f(r.u, 1)}、${f(r.se, 1)} kPa`],
        ['r_d', f(r.rd, 3)],
        ['CSR = 0.65(a_{max}/g)(σ_{v0}/σ′_{v0})r_d', f(r.csr, 3)],
        ['N_{60} = N_m(ER/60)C_R', `${f(r.N60, 2)}（C_R = ${r.CR}）`],
        ['C_N → (N_1)_{60}', `${f(r.CN, 3)} → ${f(r.N1, 2)}`],
        ['(N_1)_{60cs}', `${f(r.Ncs, 2)}（${r.note}）`],
        ['CRR_{7.5}', isFinite(r.crr75) ? f(r.crr75, 3) : '∞（太密，不液化）'],
        ['MSF、K_σ', `${f(r.msf, 3)}、${f(r.Ks, 3)}`],
        ['CRR = CRR_{7.5}·MSF·K_σ', isFinite(r.crr) ? f(r.crr, 3) : '∞']
      ];
      tw.innerHTML = `<table class="data compact">${rows.map(x => `<tr><th>${SD.rich(x[0])}</th><td class="num">${SD.rich(x[1])}</td></tr>`).join('')}</table>`;
      st({ csr: f(r.csr, 3), crr: isFinite(r.crr) ? f(r.crr, 3) : '∞', fs: fsTxt(r.fs), need: f(1.2 * r.csr / (r.msf * r.Ks), 3) },
        { fs: r.fs < 1 });
    }
    draw();
  })();

  /* =========================================================
     4. #w-boring（14.2）：整孔 FS_L、P_L、沉陷與側潰
     ========================================================= */
  (() => {
    const body = SD.widget('w-boring'); if (!body) return;
    const DEF = [  // [底深 m, N_m, FC %, 類別]
      [2, 9, 15, 's'], [3, 7, 15, 's'], [5, 9, 12, 's'], [7, 11, 12, 's'], [9, 13, 12, 's'],
      [11, 18, 8, 's'], [13, 22, 5, 's'], [16, 27, 5, 's'], [20, 32, 5, 's']
    ];
    let layers = DEF.map(r => r.slice());
    const ctl = SD.controls(body, [
      { k: 'method', type: 'select', label: '方法', value: 'youd', options: [['youd', 'Youd et al. (2001)'], ['hbf', 'HBF (2012)'], ['bi', 'Boulanger & Idriss (2014)']] },
      { k: 'gwt', label: '地下水位 (m)', min: 0, max: 6, step: 0.5, value: 1, fmt: v => v.toFixed(1) + ' m' },
      { k: 'amax', label: '地表 a_max (g)', min: 0.05, max: 0.6, step: 0.01, value: 0.24, fmt: v => v.toFixed(2) + ' g' },
      { k: 'Mw', label: '規模 M_w', min: 5.5, max: 8.5, step: 0.1, value: 7.0, fmt: v => v.toFixed(1) },
      { k: 'ER', label: '能量比 ER (%)', min: 40, max: 100, step: 5, value: 70, fmt: v => v.toFixed(0) + ' %' },
      { k: 'S', label: '地表坡度 S (%)（側潰）', min: 0.2, max: 3.5, step: 0.1, value: 1.0, fmt: v => v.toFixed(1) + ' %' }
    ], draw);
    SD.buttons(body, '', [
      { label: '預設鑽孔', on: () => { layers = DEF.map(r => r.slice()); buildTable(); ctl.set({ gwt: 1, amax: 0.24, Mw: 7.0, ER: 70 }); } },
      { label: '把 3–9 m 改良到 N_m = 25', on: () => { layers.forEach(r => { if (r[0] > 3 && r[0] <= 9) r[1] = 25; }); buildTable(); draw(); } },
      { label: '地下水位降到 4 m', on: () => ctl.set({ gwt: 4 }) }
    ]);
    const tw = SD.div(body, 'tablewrap');
    function buildTable() {
      let top = 0;
      tw.innerHTML = `<table class="data edit"><tr><th>深度 (m)</th><th>N<sub>m</sub></th><th>FC (%)</th><th>土類</th></tr>${layers.map((r, i) => {
        const h = `<tr><td class="lab">${top}–${r[0]}</td><td><input data-i="${i}" data-c="1" type="number" value="${r[1]}"></td><td><input data-i="${i}" data-c="2" type="number" value="${r[2]}"></td><td><select data-i="${i}" data-c="3"><option value="s"${r[3] === 's' ? ' selected' : ''}>砂、非塑性粉土</option><option value="c"${r[3] === 'c' ? ' selected' : ''}>黏土（PI ≥ 7，不評估）</option></select></td></tr>`;
        top = r[0]; return h;
      }).join('')}</table>`;
      tw.querySelectorAll('input,select').forEach(e => e.addEventListener('change', () => {
        const i = +e.dataset.i, c = +e.dataset.c;
        layers[i][c] = c === 3 ? e.value : (parseFloat(e.value) || 0);
        draw();
      }));
    }
    buildTable();
    const grid = SD.div(body, 'plot-row plot-grid');
    const yf = v => Math.abs(v).toFixed(0);
    const pF = new Plot(grid, { h: 320, xlim: [0, 2], ylim: [-20, 0], yfmt: yf, xlabel: 'FS_L', ylabel: '深度 z (m)' });
    pF.legend([{ label: 'FS_L', c: 'ink', w: 1.8, marker: 'circle' }, { label: 'FS_L = 1', c: 'red', dash: true }]);
    const pW = new Plot(grid, { h: 320, xlim: [0, 10], ylim: [-20, 0], yfmt: yf, xlabel: 'F·W（藍）、W = 10 − 0.5z（灰）', ylabel: '深度 z (m)' });
    pW.legend([{ label: 'F·W，面積 = P_L', c: 'blue', w: 2 }, { label: 'W(z)', c: 'gray', dash: true }]);
    const st = SD.stats(body, [
      { k: 'pl', label: 'P_L', hl: true },
      { k: 'sev', label: '分級（Iwasaki）', hl: true },
      { k: 'liq', label: 'FS_L < 1 的總厚度' },
      { k: 'set', label: '液化後沉陷（I&B 2008）' },
      { k: 'ld', label: '側潰 LD = (S + 0.2)·LDI' }
    ]);

    function draw() {
      const dz = 0.25;
      const zs = [], fs = [], fw = [];
      let pl = 0, thick = 0, set = 0, ldi = 0;
      for (let z = dz / 2; z < 20; z += dz) {
        const lay = layers.find(r => z < r[0]);
        let F = 0, FS = NaN;
        if (lay && lay[3] === 's' && z > ctl.gwt) {
          const r = L.evalPoint({ z, gwt: ctl.gwt, g1: 18.5, g2: 19.0, amax: ctl.amax, Mw: ctl.Mw, Nm: lay[1], FC: lay[2], ER: ctl.ER }, ctl.method);
          FS = r.fs;
          F = Math.max(0, 1 - (isFinite(FS) ? FS : 9));
          if (FS < 1) thick += dz;
          const pq = L.postLiq(isFinite(FS) ? FS : 9, r.Ncs);
          set += pq.ev * dz;
          if (z <= 15) ldi += pq.gmax * dz;          // 只計 15 m 內（常用上限）
        }
        const W = 10 - 0.5 * z;
        pl += F * W * dz;
        zs.push(z); fs.push(FS); fw.push(F * W);
      }
      pF.clear();
      pF.vline(1, { c: 'red', dash: '5 4', w: 1.4 });
      pF.hline(-ctl.gwt, { c: 'blue', dash: '3 3' });
      pF.line(fs.map(v => Math.min(v, 2)), zs.map(z => -z), { c: 'ink', w: 1.8 });
      pW.clear();
      pW.line([10, 0], [0, -20], { c: 'gray', dash: '5 4' });
      const xs = [], ys = [];
      zs.forEach((z, i) => { xs.push(0, fw[i], fw[i], 0); ys.push(-(z - dz / 2), -(z - dz / 2), -(z + dz / 2), -(z + dz / 2)); });
      pW.line(xs, ys, { c: 'blue', w: 1.6 });
      const sev = pl > 15 ? '嚴重' : pl > 5 ? '中度' : '輕微';
      st({ pl: f(pl, 1), sev, liq: f(thick, 2) + ' m', set: f(set * 100, 1) + ' cm', ld: f((ctl.S + 0.2) * ldi, 2) + ' m' },
        { pl: pl > 15 });
    }
    draw();
  })();
});
