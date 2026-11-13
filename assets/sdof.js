/* =========================================================
   SDOF 計算核心（Ch4；Ch3 共振柱、Ch6 也會用到）
   符號：ω 自然頻率、ξ 阻尼比、β = ω̄/ω、ρs = p0/k
   ========================================================= */
(function () {
  const TAU = 2 * Math.PI;

  const S = {
    /** 動態放大因子 D(β, ξ) */
    D(b, x) { return 1 / Math.sqrt((1 - b * b) ** 2 + (2 * x * b) ** 2); },
    /** 相位落後 φ（rad，0 ~ π） */
    phi(b, x) { return Math.atan2(2 * x * b, 1 - b * b); },
    /** 傳遞率 TR = D √(1 + (2ξβ)²) */
    TR(b, x) { return S.D(b, x) * Math.sqrt(1 + (2 * x * b) ** 2); },
    /** 共振峰位置與高度 */
    peak(x) {
      if (2 * x * x >= 1) return null;
      return { b: Math.sqrt(1 - 2 * x * x), D: 1 / (2 * x * Math.sqrt(1 - x * x)) };
    },
    /** 由共振放大倍數反解 ξ：1/(2ξ√(1−ξ²)) = Dmax */
    xiFromDmax(Dm) {
      if (Dm <= 1) return NaN;
      let a = 1e-6, b = 0.7071;
      const f = x => 1 / (2 * x * Math.sqrt(1 - x * x)) - Dm;
      for (let i = 0; i < 80; i++) { const m = (a + b) / 2; (f(m) > 0 ? (a = m) : (b = m)); }
      return (a + b) / 2;
    },

    /** 自由振動 u(t)，含三種阻尼情形；w 為 rad/s */
    free(t, w, x, u0, v0) {
      if (x < 1 - 1e-9) {
        const wd = w * Math.sqrt(1 - x * x);
        const A = (v0 + x * w * u0) / wd, B = u0;
        return Math.exp(-x * w * t) * (A * Math.sin(wd * t) + B * Math.cos(wd * t));
      }
      if (x <= 1 + 1e-9) {
        return (u0 + (v0 + w * u0) * t) * Math.exp(-w * t);
      }
      const wh = w * Math.sqrt(x * x - 1);
      const A = (v0 + x * w * u0) / wh, B = u0;
      return Math.exp(-x * w * t) * (A * Math.sinh(wh * t) + B * Math.cosh(wh * t));
    },

    /** 由靜止開始受 p0 sin ω̄t 的反應比 R(t) = u/ρs（ξ < 1） */
    fromRest(t, w, x, b) {
      const wb = b * w;
      if (x < 1e-9 && Math.abs(b - 1) < 1e-9) {
        return 0.5 * (Math.sin(w * t) - w * t * Math.cos(w * t));
      }
      const D = S.D(b, x), ph = S.phi(b, x);
      const wd = w * Math.sqrt(1 - x * x);
      const B = D * Math.sin(ph);
      const A = (x * w * B - D * wb * Math.cos(ph)) / wd;
      return Math.exp(-x * w * t) * (A * Math.sin(wd * t) + B * Math.cos(wd * t)) + D * Math.sin(wb * t - ph);
    },

    /**
     * Newmark 平均加速度法：SDOF 受地表加速度 ag（m/s²）
     * 回傳最大相對位移
     */
    peakDisp(ag, dt, T, x) {
      const w = TAU / T, m = 1, k = w * w, c = 2 * x * w;
      const g = 0.5, bt = 0.25;
      const kh = k + g / (bt * dt) * c + m / (bt * dt * dt);
      const a1 = m / (bt * dt) + g / bt * c, a2 = m / (2 * bt) + dt * (g / (2 * bt) - 1) * c;
      let u = 0, v = 0, a = -ag[0], umax = 0;
      for (let i = 0; i < ag.length - 1; i++) {
        const dp = -(ag[i + 1] - ag[i]) * m + a1 * v + a2 * a;
        const du = dp / kh;
        const dv = g / (bt * dt) * du - g / bt * v + dt * (1 - g / (2 * bt)) * a;
        const da = du / (bt * dt * dt) - v / (bt * dt) - a / (2 * bt);
        u += du; v += dv; a += da;
        if (Math.abs(u) > umax) umax = Math.abs(u);
      }
      return umax;
    },
    /** 同上，回傳整段相對位移歷時 */
    history(ag, dt, T, x) {
      const w = TAU / T, m = 1, k = w * w, c = 2 * x * w;
      const g = 0.5, bt = 0.25;
      const kh = k + g / (bt * dt) * c + m / (bt * dt * dt);
      const a1 = m / (bt * dt) + g / bt * c, a2 = m / (2 * bt) + dt * (g / (2 * bt) - 1) * c;
      let u = 0, v = 0, a = -ag[0];
      const out = [0];
      for (let i = 0; i < ag.length - 1; i++) {
        const dp = -(ag[i + 1] - ag[i]) * m + a1 * v + a2 * a;
        const du = dp / kh;
        const dv = g / (bt * dt) * du - g / bt * v + dt * (1 - g / (2 * bt)) * a;
        const da = du / (bt * dt * dt) - v / (bt * dt) - a / (2 * bt);
        u += du; v += dv; a += da;
        out.push(u);
      }
      return out;
    },

    /** 線性內插：在 xs 遞增的資料中，找 y 跨越 level 的 x */
    crossing(xs, ys, level, fromIdx, dir) {
      for (let i = fromIdx; dir > 0 ? i < xs.length - 1 : i > 0; i += dir) {
        const j = i + dir;
        const ya = ys[i], yb = ys[j];
        if ((ya - level) * (yb - level) <= 0 && ya !== yb) {
          return xs[i] + (level - ya) / (yb - ya) * (xs[j] - xs[i]);
        }
      }
      return NaN;
    },
    /** 三點拋物線內插求峰 */
    parabolicPeak(xs, ys) {
      let i = 0;
      for (let j = 1; j < ys.length; j++) if (ys[j] > ys[i]) i = j;
      if (i === 0 || i === ys.length - 1) return { x: xs[i], y: ys[i], i };
      const [x0, x1, x2] = [xs[i - 1], xs[i], xs[i + 1]], [y0, y1, y2] = [ys[i - 1], ys[i], ys[i + 1]];
      const den = (x0 - x1) * (x0 - x2) * (x1 - x2);
      const A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den;
      const B = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den;
      const C = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / den;
      const xv = -B / (2 * A);
      return { x: xv, y: C - B * B / (4 * A), i };
    },
    TAU
  };
  window.SDOF = S;
})();
