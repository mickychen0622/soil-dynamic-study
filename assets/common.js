/* =========================================================
   共用：導覽、目錄、進度、測驗、控制項、KaTeX
   每頁 <body data-page="ch4" data-root="../">
   ========================================================= */
(function () {
  /* 學習路線：以觀念分單元，不照講義章節；作業獨立成頁 */
  const SITE = {
    title: '土壤動力學 學習網站',
    stages: ['地震輸入', '振動與頻率', '波與土壤性質', '地盤反應', '土壤液化'],
    units: [
      { n: 1, stage: 0, short: '動態問題', name: '什麼是動態問題' },
      { n: 2, stage: 0, short: '地震來源', name: '地震從哪裡來' },
      { n: 3, stage: 0, short: '地動參數', name: '描述地動：時間域參數' },
      { n: 4, stage: 1, short: 'SDOF 振動', name: '單自由度振動' },
      { n: 5, stage: 1, short: '頻率觀點', name: '頻率的觀點：頻譜、反應譜與卓越週期' },
      { n: 6, stage: 2, short: '波傳', name: '波在土中傳播' },
      { n: 7, stage: 2, short: '量 Vs', name: '量測剪力波速與 Gmax' },
      { n: 8, stage: 2, short: '非線性', name: '土壤的非線性：模數與阻尼' },
      { n: 9, stage: 2, short: '室內試驗', name: '室內動態試驗' },
      { n: 10, stage: 3, short: '線性地盤反應', name: '一維線性地盤反應' },
      { n: 11, stage: 3, short: '非線性地盤反應', name: '等值線性與非線性分析' },
      { n: 12, stage: 4, short: '液化易致性', name: '液化易致性與反覆試驗' },
      { n: 13, stage: 4, short: '觸發評估', name: '液化觸發評估：簡化法' },
      { n: 14, stage: 4, short: '液化後果', name: '液化潛能與液化後果' }
    ].map(u => Object.assign(u, { id: 'u' + String(u.n).padStart(2, '0'), href: 'units/u' + String(u.n).padStart(2, '0') + '.html' })),
    hws: [
      { id: 'hw1', short: 'HW1', name: '強地動紀錄', href: 'hw/hw1.html' },
      { id: 'hw2', short: 'HW2', name: 'MASW 與場址分類', href: 'hw/hw2.html' },
      { id: 'hw3', short: 'HW3', name: '動態土壤性質與振動', href: 'hw/hw3.html' }
    ]
  };
  /* ---------- 主題 ---------- */
  function applyTheme(t) {
    if (t === 'light' || t === 'dark') document.documentElement.setAttribute('data-theme', t);
    else document.documentElement.removeAttribute('data-theme');
  }
  try { applyTheme(localStorage.getItem('sd-theme')); } catch (e) { }

  /* ---------- 左側導覽（寬螢幕固定在左；窄螢幕由 ☰ 叫出）：
       各階段單元、三份作業，目前這一頁展開它的各節並隨捲動標示 ---------- */
  function buildSidebar(page, root) {
    const secs = [...document.querySelectorAll('main section.sec[id]')].filter(s => s.dataset.num || s.dataset.short);
    const secList = () => `<div class="sn-secs">${secs.map(s => {
      const num = s.dataset.num || '', short = s.dataset.short || s.dataset.title || '';
      const n = num && !short.startsWith(num) ? `<span class="n">${num}</span>` : '';
      return `<a href="#${s.id}">${n}<span>${short}</span></a>`;
    }).join('')}</div>`;
    const item = (o, inner) => `<a class="sn-item${o.id === page ? ' cur' : ''}" href="${root}${o.href}">${inner}</a>${o.id === page && secs.length ? secList() : ''}`;
    const units = SITE.stages.map((s, i) => `<div class="sn-stage">階段 ${i + 1} · ${s}</div>${
      SITE.units.filter(u => u.stage === i).map(u => item(u, `<b>${u.n}</b><span>${u.name}</span>`)).join('')}`).join('');
    const hws = SITE.hws.map(h => item(h, `<b>${h.short.replace('HW', 'H')}</b><span>${h.short} ${h.name}</span>`)).join('');
    const side = document.createElement('aside');
    side.className = 'sidenav';
    side.innerHTML = `<div class="sn-top"><a class="brand" href="${root}index.html">土壤動力學</a>
        <button class="icon-btn" id="themeBtn" title="切換深淺色">◐</button></div>
      <nav class="sn-scroll">
        <a class="sn-item sn-home${page === 'home' ? ' cur' : ''}" href="${root}index.html"><b>⌂</b><span>首頁 · 學習路線</span></a>
        ${units}<div class="sn-stage">作業</div>${hws}</nav>`;
    const mbar = document.createElement('header');
    mbar.className = 'mbar';
    mbar.innerHTML = `<button class="icon-btn sn-toggle" title="導覽">☰</button><a class="ln-title" href="#top">${document.body.dataset.title || '土壤動力學'}</a>`;
    const scrim = document.createElement('div');
    scrim.className = 'sn-scrim';
    document.body.prepend(scrim);
    document.body.prepend(mbar);
    document.body.prepend(side);
    document.body.classList.add('has-side');

    const setOpen = v => document.body.classList.toggle('sn-open', v);
    mbar.querySelector('.sn-toggle').addEventListener('click', () => setOpen(!document.body.classList.contains('sn-open')));
    scrim.addEventListener('click', () => setOpen(false));
    side.addEventListener('click', e => { if (e.target.closest('a')) setOpen(false); });
    document.addEventListener('keydown', e => { if (e.key === 'Escape') setOpen(false); });
    side.querySelector('#themeBtn').addEventListener('click', () => {
      let cur = null;
      try { cur = localStorage.getItem('sd-theme'); } catch (e) { }
      const dark = cur ? cur === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
      const next = dark ? 'light' : 'dark';
      try { localStorage.setItem('sd-theme', next); } catch (e) { }
      applyTheme(next);
    });

    // 目前所在節：側欄高亮，並把它捲進側欄可視範圍
    const box = side.querySelector('.sn-scroll');
    const curItem = box.querySelector('.sn-item.cur');
    if (curItem) box.scrollTop = Math.max(0, curItem.offsetTop - 80);
    const links = new Map([...box.querySelectorAll('.sn-secs a')].map(a => [a.getAttribute('href').slice(1), a]));
    if (!links.size) return;
    const io = new IntersectionObserver(entries => {
      entries.forEach(e => {
        if (!e.isIntersecting) return;
        const a = links.get(e.target.id);
        if (!a) return;
        links.forEach(x => x.classList.remove('active'));
        a.classList.add('active');
        const top = a.offsetTop, bot = top + a.offsetHeight;
        if (top < box.scrollTop + 60 || bot > box.scrollTop + box.clientHeight - 60)
          box.scrollTo({ top: top - box.clientHeight / 2, behavior: 'smooth' });
      });
    }, { rootMargin: '-30% 0px -65% 0px' });
    secs.forEach(s => io.observe(s));
  }

  /* ---------- 寬螢幕並排：每段說明文字與它後面的互動圖／示意圖組成左右兩欄
       規則：從最近的 h3（或節首）起、最多 6 個文字區塊，接上後面的圖；
       兩個圖緊鄰時也並排。窄螢幕由 CSS 退回單欄，閱讀順序不變。 ---------- */
  function pairBlocks() {
    const isVis = n => n.nodeType === 1 && n.matches('[id^="w-"], .widget, figure.diagram');
    const isText = n => n.nodeType === 3 && n.textContent.trim();
    const MAXTEXT = 6;
    document.querySelectorAll('main section.sec').forEach(sec => {
      let buf = [], lastVis = null;
      [...sec.childNodes].forEach(n => {
        if (n.nodeType === 3 && !n.textContent.trim()) { buf.push(n); return; }
        if (n.nodeType === 8) return;
        if (isVis(n)) {
          const content = buf.filter(x => x.nodeType === 1 || isText(x));
          if (content.length) {
            // 太長就只取最後幾個區塊，前面的維持全寬
            let keep = buf;
            if (content.length > MAXTEXT) keep = buf.slice(buf.indexOf(content[content.length - MAXTEXT]));
            const pair = document.createElement('div'); pair.className = 'pair';
            const pt = document.createElement('div'); pt.className = 'pt';
            const pv = document.createElement('div'); pv.className = 'pv';
            sec.insertBefore(pair, keep[0]);
            keep.forEach(x => pt.appendChild(x));
            pv.appendChild(n);
            pair.append(pt, pv);
            lastVis = null;
          } else if (lastVis && lastVis.parentNode === sec) {
            const pair = document.createElement('div'); pair.className = 'pair pair2';
            const a = document.createElement('div'); a.className = 'pv';
            const b = document.createElement('div'); b.className = 'pv';
            sec.insertBefore(pair, lastVis);
            a.appendChild(lastVis); b.appendChild(n);
            pair.append(a, b);
            lastVis = null;
          } else lastVis = n;
          buf = [];
          return;
        }
        lastVis = null;
        if (n.nodeType === 1 && n.matches('.sec-head, .src, .quiz, h2, .fcards, ul.pitfalls, .grid2, .grid3')) { buf = []; return; }
        if (n.nodeType === 1 && n.matches('h3')) { buf = [n]; return; }
        buf.push(n);
      });
    });
  }

  /* ---------- 測驗 ---------- */
  function initQuizzes() {
    document.querySelectorAll('.quiz').forEach(q => {
      const fb = q.querySelector('.fb');
      const show = (ok) => {
        fb.classList.add('show');
        fb.classList.toggle('ok', ok);
        fb.classList.toggle('no', !ok);
        let v = fb.querySelector('.verdict');
        if (!v) { v = document.createElement('span'); v.className = 'verdict'; fb.prepend(v); }
        v.textContent = ok ? '答對了。' : '再想想，解說如下：';
      };
      if (q.dataset.answer !== undefined) {
        const ans = +q.dataset.answer;
        const opts = [...q.querySelectorAll('.opt')];
        opts.forEach((b, i) => b.addEventListener('click', () => {
          opts.forEach(o => o.classList.remove('right', 'wrong'));
          b.classList.add(i === ans ? 'right' : 'wrong');
          opts[ans].classList.add('right');
          show(i === ans);
        }));
      } else if (q.dataset.num !== undefined) {
        const target = +q.dataset.num, tol = +(q.dataset.tol || 0.02);
        const inp = q.querySelector('input.ans');
        const chk = () => {
          const v = parseFloat(inp.value);
          if (!isFinite(v)) return;
          show(Math.abs(v - target) <= Math.abs(target) * tol);
        };
        q.querySelector('.btn').addEventListener('click', chk);
        inp.addEventListener('keydown', e => { if (e.key === 'Enter') chk(); });
      }
    });
  }

  /* ---------- 互動元件工具 ---------- */
  const SD = {
    widget(id) {
      const host = document.getElementById(id);
      if (!host) return null;
      host.classList.add('widget');
      host.innerHTML = `<div class="w-head"><span class="tag">互動</span><span class="wt">${host.dataset.title || ''}</span></div><div class="w-body"></div>`;
      return host.querySelector('.w-body');
    },
    div(parent, cls, html) {
      const d = document.createElement('div');
      if (cls) d.className = cls;
      if (html) d.innerHTML = html;
      parent.appendChild(d);
      return d;
    },
    controls(parent, defs, onChange) {
      const box = SD.div(parent, 'controls');
      const vals = {}, setters = {};
      defs.forEach(d => {
        const row = SD.div(box, 'ctrl' + (d.type === 'check' ? ' check' : ''));
        const id = 'c' + Math.random().toString(36).slice(2, 8);
        if (d.type === 'select') {
          // <option> 不能放上下標：把 V_s30、G_{max} 寫成 Vs30、Gmax
          const plain = t => SD.rich(t).replace(/<\/?(sub|sup)>/g, '');
          row.innerHTML = `<label for="${id}">${d.label}</label><select id="${id}">${d.options.map(o => `<option value="${o[0]}">${plain(o[1])}</option>`).join('')}</select>`;
          const s = row.querySelector('select');
          s.value = d.value; vals[d.k] = d.value;
          s.addEventListener('change', () => { vals[d.k] = s.value; onChange && onChange(d.k); });
          setters[d.k] = v => { s.value = v; vals[d.k] = v; };
        } else if (d.type === 'check') {
          row.innerHTML = `<input type="checkbox" id="${id}"><label for="${id}">${d.label}</label>`;
          const c = row.querySelector('input');
          c.checked = !!d.value; vals[d.k] = !!d.value;
          c.addEventListener('change', () => { vals[d.k] = c.checked; onChange && onChange(d.k); });
          setters[d.k] = v => { c.checked = !!v; vals[d.k] = !!v; };
        } else if (d.type === 'number') {
          row.innerHTML = `<label for="${id}">${d.label}</label><input type="number" id="${id}" step="${d.step || 'any'}">`;
          const n = row.querySelector('input');
          n.value = d.value; vals[d.k] = +d.value;
          n.addEventListener('input', () => { if (isFinite(parseFloat(n.value))) { vals[d.k] = parseFloat(n.value); onChange && onChange(d.k); } });
          setters[d.k] = v => { n.value = v; vals[d.k] = +v; };
        } else {
          const fmt = d.fmt || (v => v);
          const toSlider = d.log ? v => 1000 * Math.log(v / d.min) / Math.log(d.max / d.min) : v => v;
          const fromSlider = d.log ? s => d.min * Math.pow(d.max / d.min, s / 1000) : s => +s;
          row.innerHTML = `<label for="${id}">${d.label}</label><input type="range" id="${id}"
            min="${d.log ? 0 : d.min}" max="${d.log ? 1000 : d.max}" step="${d.log ? 1 : d.step}"><output></output>`;
          const r = row.querySelector('input'), o = row.querySelector('output');
          const put = v => {
            vals[d.k] = v; o.textContent = fmt(v);
            const lo = +r.min, hi = +r.max;
            r.style.setProperty('--pct', (100 * (r.value - lo) / (hi - lo)) + '%');
          };
          r.value = toSlider(d.value); put(d.value);
          r.addEventListener('input', () => {
            let v = fromSlider(r.value);
            if (d.log && d.snap) v = +v.toPrecision(d.snap);
            put(v); onChange && onChange(d.k);
          });
          setters[d.k] = v => { r.value = toSlider(v); put(v); };
        }
      });
      vals.set = (obj) => { for (const k in obj) setters[k] && setters[k](obj[k]); onChange && onChange('*'); };
      return vals;
    },
    buttons(parent, label, items) {
      const row = SD.div(parent, 'btnrow');
      if (label) { const s = document.createElement('span'); s.className = 'lbl'; s.textContent = label; row.appendChild(s); }
      const btns = items.map(it => {
        const b = document.createElement('button');
        b.className = 'btn' + (it.primary ? ' primary' : '');
        b.innerHTML = SD.rich(it.label);
        b.addEventListener('click', () => it.on(b));
        row.appendChild(b);
        return b;
      });
      return { row, btns };
    },
    stats(parent, defs) {
      const box = SD.div(parent, 'stats');
      const map = {};
      defs.forEach(d => {
        const s = SD.div(box, 'stat' + (d.hl ? ' hl' : ''), `<div class="k">${SD.rich(d.label)}</div><div class="v">–</div>`);
        map[d.k] = s;
      });
      return (obj, warn) => {
        for (const k in obj) if (map[k]) map[k].querySelector('.v').innerHTML = SD.rich(obj[k]);
        if (warn) for (const k in map) map[k].classList.toggle('warn', !!warn[k]);
      };
    },
    tabs(parent, names, onSwitch) {
      const bar = SD.div(parent, 'tabs');
      const panes = names.map(() => SD.div(parent, 'tabpane'));
      const btns = names.map((n, i) => {
        const b = document.createElement('button');
        b.textContent = n;
        b.addEventListener('click', () => sel(i));
        bar.appendChild(b);
        return b;
      });
      function sel(i) {
        btns.forEach((b, j) => b.classList.toggle('on', i === j));
        panes.forEach((p, j) => p.classList.toggle('on', i === j));
        onSwitch && onSwitch(i);
      }
      sel(0);
      return panes;
    },
    /**
     * 符號上下標：V_s、γ_c、R_rup、PGA_in、W_{D}、e^{−ξωt}、x^2 → HTML
     * 不動檔名與識別字（A1_resonant_column.py、B1_process.xlsx、masing_model_tool.html）：
     * 底線前的符號必須是獨立的 1–3 個大寫字母、1 個小寫字母或希臘字母，且前面不是英數字、底線、點、斜線。
     */
    rich(s) {
      return String(s)
        .replace(/((?<![A-Za-z_\\])[A-Z]{1,3}′?|(?<![A-Za-z0-9_\\])[a-z]′?|(?<=[0-9])[a-z](?=_(?:\{[^}]*\}|[A-Za-z0-9]{1,2}(?![A-Za-z0-9_])))|[α-ωΑ-Ωϕϑ]′?|\))_(\{[^}]*\}|[A-Za-z0-9]{1,6}(?![A-Za-z0-9_]|\.[A-Za-z]))/g,
          (m, base, sub) => `${base}<sub>${sub[0] === '{' ? sub.slice(1, -1) : sub}</sub>`)
        .replace(/\^\{([^}]*)\}/g, '<sup>$1</sup>')
        .replace(/([A-Za-z0-9α-ω)])\^([0-9A-Za-z−-])/g, '$1<sup>$2</sup>');
    },
    /** 對已渲染的頁面文字做一次 rich（略過 KaTeX、程式碼、SVG、表單選項） */
    richify(root) {
      const w = document.createTreeWalker(root || document.body, NodeFilter.SHOW_TEXT);
      const todo = [];
      let n;
      while ((n = w.nextNode())) {
        const p = n.parentElement;
        if (!p || p.closest('.katex, script, style, code, pre, svg, select, option, textarea, sub, sup')) continue;
        if (/[_^]/.test(n.textContent)) todo.push(n);
      }
      todo.forEach(t => {
        const esc = t.textContent.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
        const html = SD.rich(esc);
        if (html === esc) return;
        const span = document.createElement('span');
        span.innerHTML = html;
        t.replaceWith(...span.childNodes);
      });
    },
    f(v, d) { return (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0)) ? v.toExponential(2) : v.toFixed(d === undefined ? 3 : d); },
    rng(seed) { // 可重現的亂數
      let s = seed >>> 0 || 1;
      const u = () => { s ^= s << 13; s >>>= 0; s ^= s >> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
      u.normal = () => { let a = 0; while (a === 0) a = u(); return Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * u()); };
      return u;
    }
  };
  window.SD = SD;
  window.SITE = SITE;
  // 在各頁元件程式執行前先排好並排區塊（common.js 位於 </main> 之後，此時內容已解析）
  if (document.querySelector('main')) pairBlocks();

  /* ---------- 啟動 ---------- */
  document.addEventListener('DOMContentLoaded', () => {
    const page = document.body.dataset.page || '';
    const root = document.body.dataset.root || '';
    buildSidebar(page, root);
    initQuizzes();
    if (window.renderMathInElement) {
      renderMathInElement(document.body, {
        delimiters: [{ left: '\\[', right: '\\]', display: true }, { left: '\\(', right: '\\)', display: false }],
        throwOnError: false, ignoredClasses: ['widget']
      });
    }
    document.dispatchEvent(new Event('sd:ready'));
    // 各章元件建好之後，把頁面上殘留的 a_b 寫法轉成上下標
    SD.richify(document.querySelector('main') || document.body);
  });
})();
