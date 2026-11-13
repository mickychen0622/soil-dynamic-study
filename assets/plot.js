/* =========================================================
   Plot：輕量 SVG 繪圖（依作圖規範）
   - 四邊框線、刻度朝內、上右也有刻度、不畫格線
   - 圖例置於圖上方（不壓曲線）
   - 顏色 class：ink（黑）、blue、red、gray
   - 標籤支援 a_{b}、a^{b} 上下標
   ========================================================= */
(function () {
  const NS = 'http://www.w3.org/2000/svg';
  let uid = 0;

  function el(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    if (attrs) for (const k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function richText(t, s) {
    s = String(s);
    const re = /([_^])\{([^}]*)\}|([_^])([^\s{])/g;
    let last = 0, m;
    const add = (str, kind) => {
      if (!str) return;
      const sp = el('tspan', null, t);
      sp.textContent = str;
      if (kind) {
        sp.setAttribute('baseline-shift', kind === '_' ? '-25%' : '38%');
        sp.setAttribute('font-size', '72%');
      }
    };
    while ((m = re.exec(s))) {
      add(s.slice(last, m.index));
      add(m[2] !== undefined ? m[2] : m[4], m[1] || m[3]);
      last = re.lastIndex;
    }
    add(s.slice(last));
  }

  function niceStep(range, n) {
    const raw = range / n;
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    const f = raw / p;
    return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
  }
  function linTicks(a, b, n) {
    const s = niceStep(b - a, n || 5);
    const out = [];
    for (let v = Math.ceil(a / s - 1e-9) * s; v <= b + s * 1e-9; v += s) out.push(+v.toPrecision(12));
    return { major: out, step: s };
  }
  function fmt(v, step) {
    if (Math.abs(v) < 1e-12) return '0';
    const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    return v.toFixed(Math.min(d, 6));
  }
  function logLabel(k) {
    const v = Math.pow(10, k);
    if (k >= -3 && k <= 4) return v >= 1 ? String(Math.round(v)) : v.toFixed(-k);
    return '10^{' + k + '}';
  }

  class Plot {
    constructor(host, o) {
      this.o = Object.assign({
        w: 640, h: 340, m: { l: 64, r: 18, t: 14, b: 48 },
        xlim: [0, 1], ylim: [0, 1], logx: false, logy: false, xlabel: '', ylabel: '',
        xticks: null, yticks: null, nx: 6, ny: 5, xfmt: null, yfmt: null
      }, o || {});
      this.o.m = Object.assign({ l: 64, r: 18, t: 14, b: 48 }, (o && o.m) || {});
      this.host = typeof host === 'string' ? document.querySelector(host) : host;
      this.wrap = document.createElement('div');
      this.wrap.className = 'plot';
      this.legendEl = document.createElement('div');
      this.legendEl.className = 'plot-legend';
      this.wrap.appendChild(this.legendEl);
      const { w, h } = this.o;
      this.svg = el('svg', { viewBox: `0 0 ${w} ${h}`, class: 'plot-svg', role: 'img' }, this.wrap);
      this.id = 'p' + (++uid);
      const defs = el('defs', null, this.svg);
      const cp = el('clipPath', { id: this.id + 'c' }, defs);
      this.clipRect = el('rect', null, cp);
      this.bandG = el('g', { 'clip-path': `url(#${this.id}c)` }, this.svg);
      this.dataG = el('g', { 'clip-path': `url(#${this.id}c)` }, this.svg);
      this.axesG = el('g', null, this.svg);
      this.overG = el('g', null, this.svg);
      this.readEl = null;
      this.host.appendChild(this.wrap);
      this.drawAxes();
    }
    get x0() { return this.o.m.l; }
    get x1() { return this.o.w - this.o.m.r; }
    get y0() { return this.o.h - this.o.m.b; }
    get y1() { return this.o.m.t; }
    tx(v) { return this.o.logx ? Math.log10(v) : v; }
    ty(v) { return this.o.logy ? Math.log10(v) : v; }
    sx(x) {
      const [a, b] = this.o.xlim.map(v => this.tx(v));
      return this.x0 + (this.tx(x) - a) / (b - a) * (this.x1 - this.x0);
    }
    sy(y) {
      const [a, b] = this.o.ylim.map(v => this.ty(v));
      return this.y0 - (this.ty(y) - a) / (b - a) * (this.y0 - this.y1);
    }
    ix(px) {
      const [a, b] = this.o.xlim.map(v => this.tx(v));
      const t = a + (px - this.x0) / (this.x1 - this.x0) * (b - a);
      return this.o.logx ? Math.pow(10, t) : t;
    }
    setLim(xlim, ylim) {
      if (xlim) this.o.xlim = xlim;
      if (ylim) this.o.ylim = ylim;
      this.drawAxes();
    }
    _ticks(axis) {
      const log = axis === 'x' ? this.o.logx : this.o.logy;
      const lim = axis === 'x' ? this.o.xlim : this.o.ylim;
      const given = axis === 'x' ? this.o.xticks : this.o.yticks;
      const userFmt = axis === 'x' ? this.o.xfmt : this.o.yfmt;
      if (log) {
        const ka = Math.floor(Math.log10(lim[0]) + 1e-9), kb = Math.ceil(Math.log10(lim[1]) - 1e-9);
        const major = [], minor = [];
        for (let k = ka; k <= kb; k++) {
          const v = Math.pow(10, k);
          if (v >= lim[0] * (1 - 1e-9) && v <= lim[1] * (1 + 1e-9)) major.push({ v, s: userFmt ? userFmt(v) : logLabel(k) });
          for (let j = 2; j <= 9; j++) { const u = j * v; if (u > lim[0] && u < lim[1]) minor.push(u); }
        }
        return { major, minor };
      }
      const nn = axis === 'x' ? this.o.nx : this.o.ny;
      const t = given ? { major: given, step: given.length > 1 ? Math.abs(given[1] - given[0]) : 1 } : linTicks(lim[0], lim[1], nn);
      return { major: t.major.map(v => ({ v, s: userFmt ? userFmt(v) : fmt(v, t.step) })), minor: [] };
    }
    drawAxes() {
      const g = this.axesG;
      while (g.firstChild) g.removeChild(g.firstChild);
      const { x0, x1, y0, y1 } = this;
      this.clipRect.setAttribute('x', x0); this.clipRect.setAttribute('y', y1);
      this.clipRect.setAttribute('width', x1 - x0); this.clipRect.setAttribute('height', y0 - y1);
      el('rect', { x: x0, y: y1, width: x1 - x0, height: y0 - y1, class: 'axis-frame' }, g);
      const L = 6, Lm = 3.5;
      const xt = this._ticks('x'), yt = this._ticks('y');
      xt.major.forEach(({ v, s }) => {
        const X = this.sx(v);
        if (X < x0 - .5 || X > x1 + .5) return;
        el('line', { x1: X, x2: X, y1: y0, y2: y0 - L, class: 'tick' }, g);
        el('line', { x1: X, x2: X, y1: y1, y2: y1 + L, class: 'tick' }, g);
        const t = el('text', { x: X, y: y0 + 17, 'text-anchor': 'middle', class: 'tick-label' }, g);
        richText(t, s);
      });
      xt.minor.forEach(v => {
        const X = this.sx(v);
        el('line', { x1: X, x2: X, y1: y0, y2: y0 - Lm, class: 'tick' }, g);
        el('line', { x1: X, x2: X, y1: y1, y2: y1 + Lm, class: 'tick' }, g);
      });
      yt.major.forEach(({ v, s }) => {
        const Y = this.sy(v);
        if (Y < y1 - .5 || Y > y0 + .5) return;
        el('line', { x1: x0, x2: x0 + L, y1: Y, y2: Y, class: 'tick' }, g);
        el('line', { x1: x1, x2: x1 - L, y1: Y, y2: Y, class: 'tick' }, g);
        const t = el('text', { x: x0 - 7, y: Y + 4.5, 'text-anchor': 'end', class: 'tick-label' }, g);
        richText(t, s);
      });
      yt.minor.forEach(v => {
        const Y = this.sy(v);
        el('line', { x1: x0, x2: x0 + Lm, y1: Y, y2: Y, class: 'tick' }, g);
        el('line', { x1: x1, x2: x1 - Lm, y1: Y, y2: Y, class: 'tick' }, g);
      });
      if (this.o.xlabel) {
        const t = el('text', { x: (x0 + x1) / 2, y: this.o.h - 8, 'text-anchor': 'middle', class: 'axis-label' }, g);
        richText(t, this.o.xlabel);
      }
      if (this.o.ylabel) {
        const cy = (y0 + y1) / 2;
        const t = el('text', { x: 16, y: cy, 'text-anchor': 'middle', class: 'axis-label', transform: `rotate(-90 16 ${cy})` }, g);
        richText(t, this.o.ylabel);
      }
    }
    clear() {
      [this.dataG, this.overG, this.bandG].forEach(g => { while (g.firstChild) g.removeChild(g.firstChild); });
      return this;
    }
    path(xs, ys) {
      let d = '', pen = false;
      for (let i = 0; i < xs.length; i++) {
        const x = xs[i], y = ys[i];
        const ok = isFinite(x) && isFinite(y) && !(this.o.logx && x <= 0) && !(this.o.logy && y <= 0);
        if (!ok) { pen = false; continue; }
        let Y = this.sy(y);
        Y = Math.max(-5000, Math.min(5000, Y));
        d += (pen ? 'L' : 'M') + this.sx(x).toFixed(2) + ' ' + Y.toFixed(2);
        pen = true;
      }
      return d;
    }
    line(xs, ys, opt) {
      opt = Object.assign({ c: 'ink', w: 1.33, dash: null, op: 1 }, opt || {});
      const dash = opt.dash === true ? '6 4' : opt.dash;
      return el('path', {
        d: this.path(xs, ys), fill: 'none', class: 's-' + opt.c, 'stroke-width': opt.w,
        'stroke-dasharray': dash, 'stroke-linejoin': 'round', 'stroke-linecap': 'round', opacity: opt.op
      }, opt.over ? this.overG : this.dataG);
    }
    fn(f, opt) {
      opt = opt || {};
      const n = opt.n || 400;
      const [a, b] = opt.range || this.o.xlim;
      const xs = [], ys = [];
      for (let i = 0; i <= n; i++) {
        const x = this.o.logx ? a * Math.pow(b / a, i / n) : a + (b - a) * i / n;
        xs.push(x); ys.push(f(x));
      }
      return this.line(xs, ys, opt);
    }
    markers(xs, ys, opt) {
      opt = Object.assign({ c: 'ink', shape: 'circle', fill: false, r: 4.2, over: false }, opt || {});
      const g = el('g', null, opt.over ? this.overG : this.dataG);
      for (let i = 0; i < xs.length; i++) {
        if (!isFinite(xs[i]) || !isFinite(ys[i])) continue;
        const X = this.sx(xs[i]), Y = this.sy(ys[i]);
        const cls = 's-' + opt.c + ' ' + (opt.fill ? 'f-' + opt.c : 'f-surface');
        const r = opt.r;
        if (opt.shape === 'circle') el('circle', { cx: X, cy: Y, r, class: cls, 'stroke-width': 1.2 }, g);
        else if (opt.shape === 'square') el('rect', { x: X - r * .9, y: Y - r * .9, width: r * 1.8, height: r * 1.8, class: cls, 'stroke-width': 1.2 }, g);
        else if (opt.shape === 'tri') el('path', { d: `M${X} ${Y - r * 1.1}L${X + r} ${Y + r * .75}L${X - r} ${Y + r * .75}Z`, class: cls, 'stroke-width': 1.2 }, g);
      }
      return g;
    }
    vline(x, opt) {
      opt = Object.assign({ c: 'gray', w: 1, dash: '4 3' }, opt || {});
      return el('line', { x1: this.sx(x), x2: this.sx(x), y1: this.y0, y2: this.y1, class: 's-' + opt.c, 'stroke-width': opt.w, 'stroke-dasharray': opt.dash }, this.dataG);
    }
    hline(y, opt) {
      opt = Object.assign({ c: 'gray', w: 1, dash: '4 3' }, opt || {});
      return el('line', { x1: this.x0, x2: this.x1, y1: this.sy(y), y2: this.sy(y), class: 's-' + opt.c, 'stroke-width': opt.w, 'stroke-dasharray': opt.dash }, this.dataG);
    }
    band(xa, xb) {
      const X0 = this.sx(xa), X1 = this.sx(xb);
      return el('rect', { x: Math.min(X0, X1), y: this.y1, width: Math.abs(X1 - X0), height: this.y0 - this.y1, class: 'f-band' }, this.bandG);
    }
    text(x, y, s, opt) {
      opt = Object.assign({ anchor: 'start', cls: 'txt', dx: 0, dy: 0 }, opt || {});
      const t = el('text', { x: this.sx(x) + opt.dx, y: this.sy(y) + opt.dy, 'text-anchor': opt.anchor, class: opt.cls }, this.overG);
      richText(t, s);
      return t;
    }
    legend(items) {
      this.legendEl.innerHTML = '';
      items.forEach(it => {
        const span = document.createElement('span');
        span.className = 'li';
        const s = el('svg', { viewBox: '0 0 30 12' });
        if (it.line !== false) {
          el('line', { x1: 1, x2: 29, y1: 6, y2: 6, class: 's-' + (it.c || 'ink'), 'stroke-width': it.w || 1.5, 'stroke-dasharray': it.dash === true ? '5 3' : it.dash }, s);
        }
        if (it.marker) {
          const cls = 's-' + (it.c || 'ink') + ' ' + (it.fill ? 'f-' + (it.c || 'ink') : 'f-surface');
          if (it.marker === 'square') el('rect', { x: 11.5, y: 2.5, width: 7, height: 7, class: cls, 'stroke-width': 1.2 }, s);
          else if (it.marker === 'tri') el('path', { d: 'M15 1.8L19 9.6L11 9.6Z', class: cls, 'stroke-width': 1.2 }, s);
          else el('circle', { cx: 15, cy: 6, r: 3.6, class: cls, 'stroke-width': 1.2 }, s);
        }
        span.appendChild(s);
        const lab = document.createElement('span');
        lab.innerHTML = window.SD && SD.rich ? SD.rich(it.label) : it.label;
        span.appendChild(lab);
        this.legendEl.appendChild(span);
      });
      return this;
    }
    readout(fn) {
      if (!this.readEl) {
        this.readEl = document.createElement('div');
        this.readEl.className = 'plot-readout';
        this.wrap.appendChild(this.readEl);
        const hit = el('rect', { x: this.x0, y: this.y1, width: this.x1 - this.x0, height: this.y0 - this.y1, fill: 'transparent', style: 'cursor:crosshair' }, this.svg);
        const guide = el('line', { class: 'guide', y1: this.y1, y2: this.y0, visibility: 'hidden' }, this.svg);
        const move = (ev) => {
          const pt = this.svg.createSVGPoint();
          const e = ev.touches ? ev.touches[0] : ev;
          pt.x = e.clientX; pt.y = e.clientY;
          const p = pt.matrixTransform(this.svg.getScreenCTM().inverse());
          if (p.x < this.x0 || p.x > this.x1) return;
          guide.setAttribute('x1', p.x); guide.setAttribute('x2', p.x); guide.setAttribute('visibility', 'visible');
          this.readEl.innerHTML = this._readFn(this.ix(p.x));
        };
        hit.addEventListener('pointermove', move);
        hit.addEventListener('pointerdown', move);
        hit.addEventListener('pointerleave', () => { guide.setAttribute('visibility', 'hidden'); this.readEl.innerHTML = this._readHint || ''; });
      }
      this._readFn = fn;
      this._readHint = '將游標移到圖上可讀值';
      this.readEl.innerHTML = this._readHint;
      return this;
    }
  }

  window.Plot = Plot;
  window.svgEl = el;
  window.svgRich = richText;
})();
