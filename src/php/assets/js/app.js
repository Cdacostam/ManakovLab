/* Manakov Lab — front-end controller
 * - form state, presets, live LaTeX previews (server-side parser via api/latex.php)
 * - runs simulations through api/run.php
 * - renders results with Plotly (3D surface, space–time maps, final profiles,
 *   conservation diagnostics, errors and convergence) and exports data.
 */
(() => {
  'use strict';

  // ------------------------------------------------------------------ utils
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => [...r.querySelectorAll(s)];
  const PRESETS = window.MANAKOV_PRESETS || [];
  const API = { latex: 'api/latex.php', run: 'api/run.php' };
  const COEFFS = ['beta', 'gamma', 'g', 'V'];
  const SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };

  const state = {
    model: 'classic',
    coeffLatex: { beta: '1', gamma: '1', g: '0', V: '0' },
    coeffValid: { beta: true, gamma: true, g: true, V: true },
    result: null,
    input: null,
    setup: null,
    drawn: new Set(),
    tab: 'evo',
    running: false,
  };

  const debounce = (fn, ms) => { let h; return (...a) => { clearTimeout(h); h = setTimeout(() => fn(...a), ms); }; };
  const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

  function sci(v, digits = 3) {
    if (v === null || v === undefined || !Number.isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a === 0) return '0';
    if (a >= 1e-3 && a < 1e5) return Number(v.toPrecision(digits)).toString();
    const [m, e] = v.toExponential(digits - 1).split('e');
    return `${m} × 10${String(parseInt(e, 10)).split('').map((c) => SUP[c] || c).join('')}`;
  }
  const fixed = (v, d = 3) => (v === null || v === undefined || !Number.isFinite(v) ? '—' : v.toFixed(d));

  function tex(el, s, display = false) {
    if (!window.katex) { el.textContent = s; return; }
    try { window.katex.render(s, el, { displayMode: display, throwOnError: false, strict: 'ignore' }); }
    catch (e) { el.textContent = s; }
  }
  function renderMathIn(el) {
    if (window.renderMathInElement) {
      window.renderMathInElement(el, {
        delimiters: [{ left: '\\[', right: '\\]', display: true }, { left: '\\(', right: '\\)', display: false }],
        throwOnError: false,
      });
    }
    $$('.tex[data-tex]', el).forEach((n) => tex(n, n.dataset.tex));
  }

  // ------------------------------------------------------------ form access
  const num = (id) => parseFloat($(id).value);
  function setVal(id, v) { const el = $(id); if (el && v !== undefined && v !== null) el.value = v; }

  function setModel(model) {
    state.model = model;
    $$('.segmented button').forEach((b) => {
      b.setAttribute('aria-checked', String(b.dataset.model === model));
      b.tabIndex = b.dataset.model === model ? 0 : -1;
    });
    $$('[data-models]').forEach((el) => { el.hidden = !el.dataset.models.split(' ').includes(model); });
    updateExactUI();
    updateEquation();
    updateCost();
  }

  // ----------------------------------------------------------------- presets
  function applyPreset(id) {
    const p = PRESETS.find((x) => x.id === id) || PRESETS[0];
    if (!p) return;
    $('#preset').value = p.id;
    $('#preset-desc').textContent = p.description;
    const P = p.params;
    setVal('#p-L', P.L); setVal('#p-Nx', P.Nx); setVal('#p-T', P.T);
    setVal('#p-dt', P.dt ?? 0.005); setVal('#p-alpha', P.alpha ?? 1.5);
    setVal('#p-nu', P.nu ?? 0.8); setVal('#p-N', P.N ?? 1000); setVal('#p-r', P.r ?? 'auto');
    COEFFS.forEach((c) => { $(`#c-${c}`).value = p.coeffs[c]; });
    $('#exact').checked = !!p.exact;
    $('#convergence').checked = !!p.convergence;
    $('#solitons').innerHTML = '';
    p.solitons.forEach((s) => addSoliton(s));
    setModel(p.model);
    COEFFS.forEach((c) => previewCoeff(c));
  }

  // ---------------------------------------------------------------- solitons
  function addSoliton(s = { eta: 0.5, xi: 0, x0: 0, p1re: 1, p1im: 0, p2re: 0, p2im: 0 }) {
    const host = $('#solitons');
    if (host.children.length >= 3) return;
    const node = $('#soliton-row').content.firstElementChild.cloneNode(true);
    $$('input', node).forEach((inp) => { inp.value = s[inp.dataset.k] ?? 0; });
    $('.remove', node).addEventListener('click', () => { node.remove(); renumberSolitons(); });
    host.appendChild(node);
    renderMathIn(node);
    renumberSolitons();
  }
  function renumberSolitons() {
    const rows = $$('#solitons .soliton');
    rows.forEach((r, i) => {
      $('.soliton-title', r).textContent = `Soliton ${i + 1}`;
      $('.remove', r).disabled = rows.length === 1;
    });
    $('#add-soliton').disabled = rows.length >= 3;
  }
  function readSolitons() {
    return $$('#solitons .soliton').map((r) => {
      const o = {};
      $$('input', r).forEach((inp) => { o[inp.dataset.k] = parseFloat(inp.value); });
      return o;
    });
  }

  // ------------------------------------------------------- LaTeX previews
  const previewTimers = {};
  function previewCoeff(name) {
    clearTimeout(previewTimers[name]);
    previewTimers[name] = setTimeout(async () => {
      const inp = $(`#c-${name}`);
      const box = $(`[data-coeff="${name}"] .preview`);
      const qs = new URLSearchParams({ expr: inp.value, vars: inp.dataset.vars });
      try {
        const r = await fetch(API.latex, { method: 'POST', body: qs, headers: { Accept: 'application/json' } });
        if (!r.ok) {
          // e.g. a web-application firewall answering 403/414 with an HTML page: the server is reachable
          box.classList.add('error'); inp.classList.add('invalid');
          box.textContent = `The server rejected this expression (HTTP ${r.status}). Simplify it or check its characters.`;
          state.coeffValid[name] = false;
          updateEquation();
          return;
        }
        const d = await r.json();
        if (d.ok) {
          box.classList.remove('error'); inp.classList.remove('invalid');
          tex(box, d.latex);
          state.coeffLatex[name] = d.latex; state.coeffValid[name] = true;
          if (COEFFS.every((c) => state.coeffValid[c]) && $('#alert').textContent.startsWith('Fix the coefficient')) showAlert('');
        } else {
          box.classList.add('error'); inp.classList.add('invalid');
          box.textContent = d.error + (d.position !== undefined ? `  (at character ${d.position + 1})` : '');
          state.coeffValid[name] = false;
        }
      } catch (e) {
        box.classList.add('error');
        box.textContent = 'Preview unavailable (server not reachable).';
      }
      updateEquation();
    }, 220);
  }

  // ------------------------------------------------------ equation card
  const SCHEME = {
    classic: 'Strang split-step Fourier · 2nd order',
    space: 'Split-step Fourier, Riesz symbol |k|^α · 2nd order',
  };
  // Scheme and expected order actually used: nu = 1 runs Strang; the server caps the automatic grading at r = 10.
  function schemeText(m, nu, r, order) {
    if (m !== 'caputo') return SCHEME[m];
    if (nu >= 1) return 'ν = 1: Strang split-step Fourier, uniform Δt = T/N · 2nd order';
    const p = order ?? Math.min(r * nu, 2 - nu);
    return `L1 scheme on graded mesh, r = ${fixed(r, 2)} · order min(rν, 2−ν) = ${fixed(p, 2)}`;
  }
  function gradingInput(nu) {
    const rIn = $('#p-r').value.trim();
    return rIn === '' || rIn === 'auto' ? Math.min(10, (2 - nu) / nu) : parseFloat(rIn);
  }
  function updateEquation() {
    const m = state.model;
    const exact = $('#exact').checked;
    const lhs = {
      classic: 'i\\,\\partial_t q_j+\\beta(t)\\,\\partial_x^2 q_j',
      space: 'i\\,\\partial_t q_j-\\beta(t)\\,(-\\partial_x^2)^{\\alpha/2} q_j',
      caputo: 'i\\,{}^{C}\\!D_t^{\\nu} q_j+\\beta(t)\\,\\partial_x^2 q_j',
    }[m];
    const forcing = m === 'caputo' && exact ? '+f_j(x,t)' : '';
    tex($('#eq-main'), `${lhs}+2\\gamma(t)\\big(|q_1|^2+|q_2|^2\\big)q_j+V(x,t)\\,q_j=i\\,g(t)\\,q_j${forcing},\\qquad j=1,2`, true);

    const L = state.coeffLatex;
    const derived = exact && m !== 'caputo';
    const gam = derived ? '\\beta(t)\\,e^{-2\\int_0^t g(s)\\,ds}\\;\\;(\\text{derived})' : L.gamma;
    const V = derived ? '0\\;\\;(\\text{exact mode})' : L.V;
    let extra = '';
    if (m === 'space') extra = `\\\\ \\alpha&=${fixed(num('#p-alpha'), 2)}`;
    let badge = SCHEME[m];
    if (m === 'caputo') {
      const nu = num('#p-nu');
      const r = gradingInput(nu);
      extra = nu >= 1 ? `\\\\ \\nu&=${fixed(nu, 2)}` : `\\\\ \\nu&=${fixed(nu, 2)} &\\qquad r&=${fixed(r, 3)}`;
      badge = schemeText(m, nu, r);
    }
    tex($('#eq-coeffs'),
      `\\begin{aligned}\\beta(t)&=${L.beta} &\\qquad g(t)&=${L.g}\\\\ \\gamma(t)&=${gam} &\\qquad V(x,t)&=${V}${extra}\\end{aligned}`, true);
    $('#scheme-badge').textContent = badge;
  }

  function updateExactUI() {
    const exact = $('#exact').checked;
    const m = state.model;
    $('#convergence').disabled = !exact;
    const derived = exact && m !== 'caputo';
    $('.derived-note').hidden = !derived;
    $$('[data-coeff="gamma"], [data-coeff="V"]').forEach((el) => el.classList.toggle('disabled', derived));
    const hints = {
      classic: 'Exact mode uses the first soliton. γ(t) is derived from β and g and V is switched off, so the exact solution is the gauge-transformed Manakov soliton.',
      space: 'Exact mode uses the first soliton as a stationary fractional soliton (ξ ignored for α < 2), carried through β(t), g(t) by the gauge reduction. γ is derived and V switched off.',
      caputo: 'Exact mode uses the manufactured solution q* = c(1 + t^ν + i t²)ψ(x) with ψ from the first soliton; the forcing is computed automatically. All coefficients, including V, are used.',
    };
    $('#exact-hint').textContent = exact ? hints[m] : 'Without exact mode the run reports conservation diagnostics only.';
  }

  function updateCost() {
    const Nx = num('#p-Nx'); const T = num('#p-T');
    let txt = '';
    if (state.model === 'caputo') {
      const N = num('#p-N');
      if (num('#p-nu') >= 1) {                     // nu = 1 runs the split-step solver with N uniform steps
        const sec = (N * Nx * Math.log2(Nx)) / 3.5e6 + 0.3;
        txt = `ν = 1: ${N} uniform Strang steps, Δt = T/N × ${Nx} grid points · rough estimate ${sec < 1 ? '< 1' : sec.toFixed(0)} s (PHP, no JIT)`;
      } else {
        const sec = (N * N * Nx) / 1.0e7 + 0.3;
        txt = `${N} graded time steps · work ∝ N²·Nx, memory ∝ N·Nx · rough estimate ${sec < 1 ? '< 1' : sec.toFixed(0)} s (PHP, no JIT)`;
      }
    } else {
      const steps = Math.round(T / num('#p-dt'));
      const sec = (steps * Nx * Math.log2(Nx)) / 3.5e6 + 0.3;
      txt = `${steps.toLocaleString()} time steps × ${Nx} grid points · rough estimate ${sec < 1 ? '< 1' : sec.toFixed(0)} s (PHP, no JIT)`;
    }
    if ($('#convergence').checked && $('#exact').checked) txt += ' (+ convergence runs)';
    $('#cost-hint').textContent = Number.isFinite(Nx) && Number.isFinite(T) ? txt : '';
  }

  function collectSetup() {
    const params = { L: num('#p-L'), Nx: num('#p-Nx'), T: num('#p-T') };
    if (state.model === 'caputo') {
      params.nu = num('#p-nu'); params.N = num('#p-N');
      const r = $('#p-r').value.trim();
      params.r = r === '' || r === 'auto' ? 'auto' : parseFloat(r);
    } else {
      params.dt = num('#p-dt');
      if (state.model === 'space') params.alpha = num('#p-alpha');
    }
    const coeffs = {};
    COEFFS.forEach((c) => { coeffs[c] = $(`#c-${c}`).value; });
    return {
      model: state.model, params, coeffs, solitons: readSolitons(),
      exact: $('#exact').checked, convergence: $('#exact').checked && $('#convergence').checked,
    };
  }

  // -------------------------------------------------------------------- run
  let busyTimer = null;
  function setBusy(on, detail = '') {
    state.running = on;
    $('#busy').hidden = !on;
    $('#run').disabled = on;
    $('#busy-detail').textContent = detail;
    clearInterval(busyTimer);
    if (on) {
      const t0 = performance.now();
      busyTimer = setInterval(() => { $('#busy-time').textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`; }, 100);
    }
  }
  function showAlert(msg) {
    const a = $('#alert');
    a.hidden = !msg;
    a.textContent = msg || '';
    $('#status').textContent = msg ? 'Run failed.' : '';
    $('#status').classList.toggle('error', !!msg);
  }

  async function run(ev) {
    ev?.preventDefault();
    if (state.running) return;
    const bad = COEFFS.filter((c) => !state.coeffValid[c]);
    if (bad.length) { showAlert(`Fix the coefficient expression(s): ${bad.join(', ')}.`); $(`#c-${bad[0]}`).focus(); return; }
    $$('input.invalid').forEach((i) => { if (!i.closest('.coeff')) i.classList.remove('invalid'); });
    showAlert('');
    const setup = collectSetup();
    const detail = setup.model === 'caputo'
      ? `${setup.params.N} steps, Nx = ${setup.params.Nx}` : `${Math.round(setup.params.T / setup.params.dt)} steps, Nx = ${setup.params.Nx}`;
    setBusy(true, detail + (setup.convergence ? ' · with convergence study' : ''));
    $('#status').textContent = 'Running…';
    const ctrl = new AbortController();
    const kill = setTimeout(() => ctrl.abort(), 300000);
    try {
      const r = await fetch(API.run, {
        method: 'POST', headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
        body: JSON.stringify(setup), signal: ctrl.signal,
      });
      let d;
      try { d = await r.json(); } catch (e) { throw new Error(`Server returned HTTP ${r.status} without a JSON body.`); }
      if (!d.ok) {
        if (d.field) {
          const el = $(`#c-${d.field}`) || $(`#p-${d.field}`);
          if (el) { el.classList.add('invalid'); el.focus(); }
        }
        throw new Error(d.error || 'Unknown error');
      }
      state.result = d.result; state.input = d.input; state.setup = setup;
      state.drawn.clear();
      renderAll();
      const RM = d.result.meta;      // the run itself: scheme, grading and order used by the server
      $('#scheme-badge').textContent = schemeText(d.result.model, RM.nu, RM.r, RM.expectedOrder);
      $('#status').textContent = `Done in ${d.result.meta.runtime_total.toFixed(2)} s (server).`;
      const rt = $('#res-title'); if (rt) { rt.setAttribute('tabindex', '-1'); rt.focus(); }
    } catch (e) {
      showAlert(e.name === 'AbortError' ? 'The request timed out.' : e.message);
    } finally {
      clearTimeout(kill);
      setBusy(false);
    }
  }

  // --------------------------------------------------------------- plotting
  const hasPlotly = () => !!window.Plotly;
  const PLOT_CONFIG = {
    responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
    toImageButtonOptions: { format: 'png', scale: 2, filename: 'manakov-lab' },
  };
  function isObj(v) { return v && typeof v === 'object' && !Array.isArray(v); }
  function merge(a, b) {
    const out = { ...a };
    Object.entries(b || {}).forEach(([k, v]) => { out[k] = isObj(v) && isObj(a[k]) ? merge(a[k], v) : v; });
    return out;
  }
  function theme() {
    return {
      ink: cssVar('--ink'), ink2: cssVar('--ink-2'), muted: cssVar('--muted'), grid: cssVar('--grid'),
      axis: cssVar('--axis'), surface: cssVar('--surface'), c1: cssVar('--series-1'), c2: cssVar('--series-2'),
      dark: document.documentElement.dataset.theme === 'dark',
    };
  }
  function axis(th, extra = {}) {
    return merge({
      gridcolor: th.grid, linecolor: th.axis, tickcolor: th.axis, zeroline: false, showline: true,
      ticks: 'outside', ticklen: 4, automargin: true, exponentformat: 'power', title: { font: { size: 12, color: th.ink2 } },
    }, extra);
  }
  function layout(extra = {}) {
    const th = theme();
    return merge({
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: 'system-ui, -apple-system, "Segoe UI", sans-serif', size: 12, color: th.ink2 },
      margin: { l: 64, r: 18, t: 30, b: 48 },
      xaxis: axis(th), yaxis: axis(th),
      legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', bgcolor: 'rgba(0,0,0,0)', font: { color: th.ink } },
      hovermode: 'x unified',
      hoverlabel: { bgcolor: th.surface, bordercolor: th.grid, font: { color: th.ink } },
    }, extra);
  }
  function colorscale(th) {
    return th.dark
      ? [[0, '#1a1a19'], [0.12, '#0d366b'], [0.28, '#184f95'], [0.45, '#256abf'], [0.62, '#3987e5'], [0.78, '#6da7ec'], [0.9, '#9ec5f4'], [1, '#cde2fb']]
      : [[0, '#f4f8fd'], [0.12, '#cde2fb'], [0.28, '#9ec5f4'], [0.45, '#6da7ec'], [0.62, '#3987e5'], [0.78, '#256abf'], [0.9, '#184f95'], [1, '#0d366b']];
  }
  const line = (color, dash = 'solid', width = 2) => ({ color, dash, width });
  function plot(id, traces, lay) {
    if (!hasPlotly()) return;
    window.Plotly.react(id, traces, lay, PLOT_CONFIG);
  }
  const posOrNull = (a) => a.map((v) => (v !== null && v > 0 ? v : null));

  function fieldOf(S, which) {
    if (which === 'I1') return S.I1;
    if (which === 'I2') return S.I2;
    return S.I1.map((row, i) => row.map((v, j) => v + S.I2[i][j]));
  }
  const FIELD_LABEL = { sum: '|q₁|²+|q₂|²', I1: '|q₁|²', I2: '|q₂|²' };

  function drawEvolution() {
    const R = state.result; const S = R.surface; const th = theme();
    const which = $('#evo-field').value; const view = $('#evo-view').value;
    const z = fieldOf(S, which); const lab = FIELD_LABEL[which];
    const cbar = { thickness: 10, len: 0.7, outlinewidth: 0, tickfont: { color: th.ink2, size: 11 }, title: { text: lab, side: 'right', font: { color: th.ink2, size: 12 } } };
    if (view === 'surface') {
      const sceneAxis = (title) => ({
        title: { text: title, font: { color: th.ink2, size: 13 } }, color: th.ink2, gridcolor: th.grid,
        zerolinecolor: th.axis, showbackground: true, backgroundcolor: th.dark ? '#151514' : '#f6f5f2', showspikes: false,
      });
      plot('plot-evo', [{
        type: 'surface', x: S.x, y: S.t, z, colorscale: colorscale(th), colorbar: cbar,
        lighting: { ambient: 0.72, diffuse: 0.62, specular: 0.06, roughness: 0.85, fresnel: 0.1 },
        contours: { z: { show: true, usecolormap: true, project: { z: true }, width: 1 } },
        hovertemplate: `x = %{x:.3f}<br>t = %{y:.3f}<br>${lab} = %{z:.4g}<extra></extra>`,
      }], layout({
        margin: { l: 0, r: 0, t: 0, b: 0 }, hovermode: 'closest',
        scene: {
          xaxis: sceneAxis('x'), yaxis: sceneAxis('t'), zaxis: sceneAxis(lab),
          aspectmode: 'manual', aspectratio: { x: 1.7, y: 1.25, z: 0.62 },
          camera: { eye: { x: 1.45, y: -1.75, z: 0.95 }, center: { x: 0, y: 0, z: -0.12 } },
        },
      }));
    } else {
      plot('plot-evo', [{
        type: 'heatmap', x: S.x, y: S.t, z, colorscale: colorscale(th), colorbar: cbar, zsmooth: 'best',
        hovertemplate: `x = %{x:.3f}<br>t = %{y:.3f}<br>${lab} = %{z:.4g}<extra></extra>`,
      }], layout({ hovermode: 'closest', xaxis: { title: { text: 'x' } }, yaxis: { title: { text: 't' } } }));
    }
    ['I1', 'I2'].forEach((k, i) => {
      plot(`plot-map${i + 1}`, [{
        type: 'heatmap', x: S.x, y: S.t, z: S[k], colorscale: colorscale(th), zsmooth: 'best',
        colorbar: { ...cbar, title: { text: FIELD_LABEL[k], side: 'right', font: { color: th.ink2, size: 12 } } },
        hovertemplate: `x = %{x:.3f}<br>t = %{y:.3f}<br>${FIELD_LABEL[k]} = %{z:.4g}<extra></extra>`,
      }], layout({ hovermode: 'closest', margin: { t: 10 }, xaxis: { title: { text: 'x' } }, yaxis: { title: { text: 't' } } }));
    });
  }

  function drawFinal() {
    const R = state.result; const th = theme(); const x = R.x; const F = R.final; const E = R.exact;
    const T = fixed(F.T, 3);
    $('#cap-final').textContent = `Moduli |q₁(x,T)|, |q₂(x,T)| at T = ${T}` + (E ? ' — solid: numerical, dashed: exact' : '');
    const tr = [
      { x, y: F.q1.abs, name: '|q₁| numerical', line: line(th.c1), hovertemplate: '%{y:.5g}' },
      { x, y: F.q2.abs, name: '|q₂| numerical', line: line(th.c2), hovertemplate: '%{y:.5g}' },
    ];
    if (E) {
      tr.push({ x, y: E.q1abs, name: '|q₁| exact', line: line(th.c1, 'dash', 2), opacity: 0.9, hovertemplate: '%{y:.5g}' });
      tr.push({ x, y: E.q2abs, name: '|q₂| exact', line: line(th.c2, 'dash', 2), opacity: 0.9, hovertemplate: '%{y:.5g}' });
    }
    tr.push({ x, y: R.initial.abs1, name: '|q₁| at t = 0', line: line(th.muted, 'dot', 1.5), visible: 'legendonly', hovertemplate: '%{y:.5g}' });
    tr.push({ x, y: R.initial.abs2, name: '|q₂| at t = 0', line: line(th.muted, 'dashdot', 1.5), visible: 'legendonly', hovertemplate: '%{y:.5g}' });
    plot('plot-final', tr, layout({ xaxis: { title: { text: 'x' } }, yaxis: { title: { text: `|qⱼ(x, ${T})|` }, rangemode: 'tozero' } }));

    const I = F.q1.abs.map((a, i) => a * a + F.q2.abs[i] * F.q2.abs[i]);
    const lt = [{ x, y: posOrNull(I), name: 'numerical', line: line(th.c1), hovertemplate: '%{y:.3e}' }];
    if (E) lt.push({ x, y: posOrNull(E.q1abs.map((a, i) => a * a + E.q2abs[i] * E.q2abs[i])), name: 'exact', line: line(th.ink, 'dash', 1.5), hovertemplate: '%{y:.3e}' });
    plot('plot-final-log', lt, layout({
      showlegend: !!E, xaxis: { title: { text: 'x' } },
      yaxis: { title: { text: '|q₁|²+|q₂|²' }, type: 'log', range: [-16, Math.log10(Math.max(...I)) + 0.5] },
    }));

    plot('plot-final-reim', [
      { x, y: F.q1.re, name: 'Re q₁', line: line(th.c1), hovertemplate: '%{y:.4g}' },
      { x, y: F.q1.im, name: 'Im q₁', line: line(th.c1, 'dot', 1.6), hovertemplate: '%{y:.4g}' },
      { x, y: F.q2.re, name: 'Re q₂', line: line(th.c2), hovertemplate: '%{y:.4g}' },
      { x, y: F.q2.im, name: 'Im q₂', line: line(th.c2, 'dot', 1.6), hovertemplate: '%{y:.4g}' },
    ], layout({ xaxis: { title: { text: 'x' } }, yaxis: { title: { text: `qⱼ(x, ${T})` } } }));

    // summary table
    const dx = R.meta.dx;
    const P1 = F.q1.abs.reduce((s, a) => s + a * a, 0) * dx;
    const P2 = F.q2.abs.reduce((s, a) => s + a * a, 0) * dx;
    let imax = 0; I.forEach((v, i) => { if (v > I[imax]) imax = i; });
    $('#final-table').innerHTML = `<table class="data"><thead><tr><th>Quantity at t = ${T}</th><th class="num">q₁</th><th class="num">q₂</th><th class="num">total</th></tr></thead><tbody>
      <tr><td>Power ∫|qⱼ|² dx</td><td class="num">${sci(P1, 6)}</td><td class="num">${sci(P2, 6)}</td><td class="num">${sci(P1 + P2, 6)}</td></tr>
      <tr><td>Power fraction</td><td class="num">${fixed(P1 / (P1 + P2), 4)}</td><td class="num">${fixed(P2 / (P1 + P2), 4)}</td><td class="num">1</td></tr>
      <tr><td>Peak intensity / location</td><td class="num" colspan="2">x = ${fixed(x[imax], 3)}</td><td class="num">${sci(I[imax], 5)}</td></tr>
      ${E ? `<tr><td>Error max|q − q<sub>exact</sub>|</td><td class="num" colspan="2"></td><td class="num">${sci(E.err_final, 4)}</td></tr>` : ''}
    </tbody></table>`;
  }

  function drawConservation() {
    const R = state.result; const th = theme(); const t = R.t; const D = R.diagnostics;
    const P = D.P1.map((v, i) => v + D.P2[i]);
    plot('plot-power', [
      { x: t, y: D.P1, name: 'P₁', line: line(th.c1), hovertemplate: '%{y:.6g}' },
      { x: t, y: D.P2, name: 'P₂', line: line(th.c2), hovertemplate: '%{y:.6g}' },
      { x: t, y: P, name: 'P₁ + P₂', line: line(th.ink, 'dash', 1.6), hovertemplate: '%{y:.6g}' },
    ], layout({ xaxis: { title: { text: 't' } }, yaxis: { title: { text: 'power' }, rangemode: 'tozero' } }));
    const bal = P.map((v, i) => v / (P[0] * D.gainFactor[i]) - 1);
    plot('plot-balance', [{ x: t, y: bal, name: 'deviation', line: line(th.c1), hovertemplate: '%{y:.3e}' }],
      layout({ showlegend: false, xaxis: { title: { text: 't' } }, yaxis: { title: { text: 'relative deviation' }, tickformat: '.2~e' } }));
    plot('plot-peak', [{ x: t, y: D.peak, name: 'peak', line: line(th.c1), hovertemplate: '%{y:.5g}' }],
      layout({ showlegend: false, xaxis: { title: { text: 't' } }, yaxis: { title: { text: 'max |q|²' }, rangemode: 'tozero' } }));
    const H0 = D.H[0];
    plot('plot-ham', [{ x: t, y: D.H.map((h) => (h - H0) / Math.abs(H0 || 1)), name: 'H', line: line(th.c1), hovertemplate: '%{y:.3e}' }],
      layout({ showlegend: false, xaxis: { title: { text: 't' } }, yaxis: { title: { text: '(H − H₀)/|H₀|' }, tickformat: '.2~e' } }));
  }

  function drawErrors() {
    const R = state.result; const th = theme(); const E = R.exact; const C = R.convergence;
    $('#err-empty').hidden = !!E;
    $('#err-content').hidden = !E;
    if (!E) return;
    plot('plot-err', [{ x: R.t, y: posOrNull(E.err_t), name: 'error', line: line(th.c1), hovertemplate: '%{y:.3e}' }],
      layout({ showlegend: false, xaxis: { title: { text: 't' } }, yaxis: { title: { text: 'max |q − q_exact|' }, type: 'log' } }));
    const host = $('#plot-conv');
    if (!C) {
      if (hasPlotly()) window.Plotly.purge(host);
      host.innerHTML = '<div class="empty small"><p>Enable <strong>Convergence study</strong> to measure the observed order.</p></div>';
      $('#conv-table').innerHTML = ''; $('#cap-conv').textContent = 'Convergence study';
      return;
    }
    host.innerHTML = '';
    const caputo = R.model === 'caputo';
    const pts = C.h.map((h, i) => [h, C.err[i]]).filter((p) => p[1] && p[1] > 0);
    const [hl, el] = pts[pts.length - 1];
    const p = C.expected;
    const ref = pts.map(([h]) => el * Math.pow(h / hl, caputo ? -p : p));
    $('#cap-conv').textContent = caputo ? `Convergence in the number of time steps N (expected order ${fixed(p, 2)})` : `Convergence in Δt (expected order ${fixed(p, 0)})`;
    plot('plot-conv', [
      { x: pts.map((q) => q[0]), y: pts.map((q) => q[1]), name: 'measured', mode: 'lines+markers', line: line(th.c1), marker: { size: 9, color: th.c1, line: { color: th.surface, width: 2 } }, hovertemplate: '%{y:.3e}' },
      { x: pts.map((q) => q[0]), y: ref, name: `slope ${caputo ? '−' : ''}${fixed(p, 2)}`, mode: 'lines', line: line(th.muted, 'dash', 1.6), hoverinfo: 'skip' },
    ], layout({
      hovermode: 'closest',
      xaxis: { title: { text: caputo ? 'time steps N' : 'Δt' }, type: 'log' },
      yaxis: { title: { text: 'max error' }, type: 'log' },
    }));
    const rows = C.h.map((h, i) => `<tr><td class="num">${caputo ? h : sci(h, 3)}</td><td class="num">${sci(C.err[i], 4)}</td><td class="num">${C.order[i] == null ? '—' : fixed(C.order[i], 3)}</td></tr>`).join('');
    $('#conv-table').innerHTML = `<table class="data"><thead><tr><th class="num">${caputo ? 'N' : 'Δt'}</th><th class="num">max error</th><th class="num">observed order</th></tr></thead><tbody>${rows}</tbody>
      <tfoot><tr><td colspan="3" class="muted">Expected order: ${fixed(p, 3)}${caputo ? ` (min(rν, 2 − ν) with r = ${fixed(R.meta.r, 3)})` : ' (Strang splitting)'}${floorNote(R, C)}</td></tr></tfoot></table>`;
  }

  // The time mesh actually used: T/N for the split-step solver (also for nu = 1), graded for the L1 scheme.
  function meshLabel(R) {
    const M = R.meta;
    if (M.dtEffective == null) return `graded mesh, r = ${fixed(M.r, 2)}`;
    const inv = 1 / M.dtEffective;
    const dt = Math.abs(inv - Math.round(inv)) < 1e-9 * inv ? `1/${Math.round(inv)}` : sci(M.dtEffective, 6);
    return `${R.model === 'caputo' ? 'uniform Strang steps, ' : ''}Δt = T/N = ${dt}`;
  }

  // Why some observed orders are omitted: reference floor (with its dominant part) and unresolved coefficients.
  function floorNote(R, C) {
    const parts = (R.exact && R.exact.refFloorParts) || {};
    const names = { images: 'periodic images', spatial: 'spatial resolution', profile: 'profile residual', roundoff: 'round-off' };
    const main = Object.keys(parts).sort((a, b) => parts[b] - parts[a])[0];
    let s = '';
    if (C.refFloor) s += ` · orders omitted where the error is below 10 × ${sci(C.refFloor, 2)}, the reference floor${main ? ` (mainly ${names[main] || main})` : ''}`;
    if (C.unresolved && C.unresolved.length) s += ` · orders omitted at Δt = ${C.unresolved.map((h) => sci(h, 3)).join(', ')}, which do not resolve the coefficients`;
    if (C.saturated && C.saturated.length) s += ` · orders omitted at ${R.model === 'caputo' ? 'N' : 'Δt'} = ${C.saturated.map((h) => (R.model === 'caputo' ? h : sci(h, 3))).join(', ')}, where the error stops decreasing${C.dtIndependent ? ` (step-independent part ≈ ${sci(C.dtIndependent, 2)})` : ''}`;
    if (C.preasymptotic && C.preasymptotic.length) s += ` · orders omitted where the error exceeds 10% of the peak (not asymptotic)`;
    if (C.untested && C.untested.length) s += ` · only two levels: the saturation test needs three, so no order is published (${R.model === 'caputo' ? 'use N ≥ 34' : 'use Δt ≤ T/40'})`;
    if (R.exact && R.exact.errorSampling === 'every step') s += ' · errors measured at every time step';
    return s;
  }

  function drawData() {
    const R = state.result; const M = R.meta; const I = state.input;
    const modelName = { classic: 'Variable-coefficient Manakov', space: 'Space-fractional (Riesz) Manakov', caputo: 'Time-fractional (Caputo) Manakov' }[R.model];
    const rows = [
      ['Model', modelName], ['Scheme', M.scheme],
      ['Domain', `[−${fixed(M.L / 2, 2)}, ${fixed(M.L / 2, 2)}), periodic`], ['Grid', `Nx = ${M.Nx}, Δx = ${sci(M.dx, 4)}`],
      ['Final time', `T = ${M.T}`],
      R.model === 'caputo' && M.dtEffective == null ? ['Time mesh', `N = ${M.N}, graded t_n = T(n/N)^r, r = ${fixed(M.r, 3)}`] : ['Time step', `${meshLabel(R)} (${M.steps} steps${M.dt != null && Math.abs(M.dtEffective - M.dt) > 1e-12 * M.T ? `; requested Δt = ${M.dt}` : ''})`],
    ];
    if (R.model === 'space') rows.push(['Fractional order', `α = ${M.alpha}`]);
    if (R.model === 'caputo') rows.push(['Caputo order', `ν = ${M.nu}`]);
    rows.push(['Exact mode', M.exact ? (M.gammaDerived ? 'on (γ derived from β, g)' : 'on (manufactured solution)') : 'off']);
    rows.push(['Run time', `${M.runtime_main.toFixed(2)} s solve, ${M.runtime_total.toFixed(2)} s total`]);
    rows.push(['Engine', `Manakov Lab v${M.version} · ${M.engine}`]);
    if (M.picardMax) rows.push(['Nonlinear solve', `Anderson-accelerated fixed point (depth ${M.andersonDepth}): mean ${M.picardMean}, max ${M.picardMax} iterations/step`]);
    $('#run-summary').innerHTML = `<table class="data"><tbody>${rows.map(([a, b]) => `<tr><th>${esc(a)}</th><td>${esc(b)}</td></tr>`).join('')}</tbody></table>`;
    const names = { beta: '\\beta(t)', gamma: '\\gamma(t)', g: 'g(t)', V: 'V(x,t)' };
    const used = { ...I.latex }; const py = { ...I.canonical };
    if (M.gammaDerived) {
      used.gamma = '\\beta(t)\\,e^{-2\\int_0^t g(s)\\,ds}'; py.gamma = 'derived (exact mode)';
      used.V = '0'; py.V = 'switched off (exact mode)';
    }
    $('#run-coeffs').innerHTML = `<table class="data"><tbody>${COEFFS.map((c) => `<tr><th><span class="tex" data-tex="${esc(names[c])}"></span></th><td><span class="tex" data-tex="${esc(used[c])}"></span></td><td><code>${esc(py[c])}</code></td></tr>`).join('')}</tbody></table>`;
    $$('#run-coeffs .tex').forEach((n) => tex(n, n.dataset.tex));
  }

  const DRAW = { evo: drawEvolution, final: drawFinal, cons: drawConservation, err: drawErrors, data: drawData };
  function drawTab(tab) {
    if (!state.result || state.drawn.has(tab)) return;
    DRAW[tab]();
    state.drawn.add(tab);
  }

  function renderKPIs() {
    const R = state.result; const M = R.meta; const D = R.diagnostics;
    const P = D.P1.map((v, i) => v + D.P2[i]);
    const tiles = [
      { l: 'Solve time', v: `${M.runtime_main.toFixed(2)} s`, s: `${M.runtime_total.toFixed(2)} s incl. diagnostics` },
      { l: 'Time steps', v: M.steps.toLocaleString(), s: meshLabel(R) },
      { l: 'Grid', v: `${M.Nx}`, s: `Δx = ${sci(M.dx, 3)}` },
    ];
    if (R.model === 'caputo') {
      tiles.push({ l: 'Power P(T)/P(0)', v: fixed(P[P.length - 1] / P[0], 4), s: M.nu < 1 ? 'not conserved for ν < 1' : '' });
    } else {
      const dev = Math.max(...P.map((v, i) => Math.abs(v / (P[0] * D.gainFactor[i]) - 1)));
      tiles.push({ l: 'Power-law deviation', v: sci(dev, 2), s: 'relative to P₀·exp(2∫g)', good: dev < 1e-8 });
    }
    if (R.exact) tiles.push({ l: 'Max error', v: sci(R.exact.err_max, 3), s: `at T: ${sci(R.exact.err_final, 3)}${R.exact.errorSampling === 'every step' ? ' · max over every step' : ''}` });
    if (R.convergence) {
      const o = R.convergence.order.filter((v) => v != null);
      const last = o[o.length - 1];
      const ok = last != null && Math.abs(last - R.convergence.expected) < 0.1;
      tiles.push({ l: 'Observed order', v: fixed(last, 2), s: `expected ${fixed(R.convergence.expected, 2)}${ok ? ' ✓' : ''}`, good: ok });
    }
    $('#kpis').innerHTML = tiles.map((k) => `<div class="kpi"><div class="kpi-label">${esc(k.l)}</div><div class="kpi-value">${esc(k.v)}</div><div class="kpi-sub${k.good ? ' good' : ''}">${esc(k.s)}</div></div>`).join('');
    $('#kpis').hidden = false;
  }

  function renderAll() {
    const R = state.result;
    $('#empty').hidden = true;
    $('#tabs-wrap').hidden = false;
    $('#print-report').disabled = false;
    renderKPIs();
    const w = $('#warnings');
    w.hidden = !R.warnings.length;
    w.innerHTML = R.warnings.length ? `<ul>${R.warnings.map((m) => `<li>${esc(m)}</li>`).join('')}</ul>` : '';
    drawTab(state.tab);
  }

  // ------------------------------------------------------------------- tabs
  function selectTab(tab) {
    state.tab = tab;
    $$('.tabs [role="tab"]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === tab)));
    $$('.tabpanel').forEach((p) => { p.hidden = p.dataset.panel !== tab; });
    drawTab(tab);
    if (hasPlotly()) $$(`[data-panel="${tab}"] .js-plotly-plot`).forEach((el) => window.Plotly.Plots.resize(el));
  }

  // ----------------------------------------------------------------- export
  function download(name, text, type = 'text/csv') {
    const url = URL.createObjectURL(new Blob([text], { type }));
    const a = Object.assign(document.createElement('a'), { href: url, download: name });
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function exportData(kind) {
    const R = state.result; if (!R) return;
    const tag = `${R.model}-T${R.meta.T}`;
    if (kind === 'final') {
      const F = R.final; const E = R.exact;
      const head = ['x', 're_q1', 'im_q1', 'abs_q1', 're_q2', 'im_q2', 'abs_q2'].concat(E ? ['abs_q1_exact', 'abs_q2_exact'] : []);
      const rows = R.x.map((x, i) => [x, F.q1.re[i], F.q1.im[i], F.q1.abs[i], F.q2.re[i], F.q2.im[i], F.q2.abs[i]]
        .concat(E ? [E.q1abs[i], E.q2abs[i]] : []).join(','));
      download(`manakov-final-${tag}.csv`, [head.join(','), ...rows].join('\n'));
    } else if (kind === 'diag') {
      const D = R.diagnostics; const E = R.exact;
      const head = ['t', 'P1', 'P2', 'P', 'peak_intensity', 'H', 'exp(2*int g)'].concat(E ? ['error'] : []);
      const rows = R.t.map((t, i) => [t, D.P1[i], D.P2[i], D.P1[i] + D.P2[i], D.peak[i], D.H[i], D.gainFactor[i]]
        .concat(E ? [E.err_t[i]] : []).join(','));
      download(`manakov-diagnostics-${tag}.csv`, [head.join(','), ...rows].join('\n'));
    } else if (kind === 'surface') {
      const S = R.surface; const out = ['t,x,abs_q1_sq,abs_q2_sq'];
      S.t.forEach((t, i) => S.x.forEach((x, j) => out.push(`${t},${x},${S.I1[i][j]},${S.I2[i][j]}`)));
      download(`manakov-spacetime-${tag}.csv`, out.join('\n'));
    } else if (kind === 'json') {
      download(`manakov-result-${tag}.json`, JSON.stringify({ input: state.input, result: R }), 'application/json');
    } else if (kind === 'config') {
      download(`manakov-setup-${tag}.json`, JSON.stringify(state.setup, null, 2), 'application/json');
    }
  }

  function printReport() {
    if (!state.result) return;
    const panels = $$('.tabpanel');
    panels.forEach((p) => { p.hidden = false; });
    Object.keys(DRAW).forEach((t) => drawTab(t));
    if (hasPlotly()) $$('.js-plotly-plot').forEach((el) => window.Plotly.Plots.resize(el));
    const restore = () => { selectTab(state.tab); window.removeEventListener('afterprint', restore); };
    window.addEventListener('afterprint', restore);
    setTimeout(() => window.print(), 400);
  }

  // ------------------------------------------------------------------ theme
  function toggleTheme() {
    const root = document.documentElement;
    root.dataset.theme = root.dataset.theme === 'dark' ? 'light' : 'dark';
    try { localStorage.setItem('manakov-theme', root.dataset.theme); } catch (e) { /* storage unavailable */ }
    state.drawn.clear();
    if (state.result) drawTab(state.tab);
  }

  // ------------------------------------------------------------------- init
  function init() {
    if (!hasPlotly() || window.__plotlyFailed) $('#plotly-banner').hidden = false;
    renderMathIn(document.body);

    $('#preset').addEventListener('change', (e) => applyPreset(e.target.value));
    $$('.segmented button').forEach((b) => b.addEventListener('click', () => setModel(b.dataset.model)));
    // ARIA radio group: arrow keys move the selection, only the checked radio is in the tab order
    $('.segmented').addEventListener('keydown', (e) => {
      const keys = { ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1 };
      if (!(e.key in keys)) return;
      e.preventDefault();
      const bs = $$('.segmented button'); const i = bs.findIndex((b) => b.dataset.model === state.model);
      const n = bs[(i + keys[e.key] + bs.length) % bs.length];
      setModel(n.dataset.model); n.focus();
    });
    COEFFS.forEach((c) => $(`#c-${c}`).addEventListener('input', () => previewCoeff(c)));
    $('#add-soliton').addEventListener('click', () => addSoliton());
    $('#exact').addEventListener('change', () => { updateExactUI(); updateEquation(); updateCost(); });
    $('#convergence').addEventListener('change', updateCost);
    ['#p-L', '#p-Nx', '#p-T', '#p-dt', '#p-N', '#p-nu', '#p-alpha', '#p-r'].forEach((id) => $(id).addEventListener('input', debounce(() => { updateCost(); updateEquation(); $(id).classList.remove('invalid'); }, 150)));
    $('#setup').addEventListener('submit', run);
    $$('.tabs [role="tab"]').forEach((b) => b.addEventListener('click', () => selectTab(b.dataset.tab)));
    $('.tabs').addEventListener('keydown', (e) => {
      if (e.key !== 'ArrowRight' && e.key !== 'ArrowLeft') return;
      const tabs = $$('.tabs [role="tab"]'); const i = tabs.findIndex((b) => b.dataset.tab === state.tab);
      const n = tabs[(i + (e.key === 'ArrowRight' ? 1 : tabs.length - 1)) % tabs.length];
      n.focus(); selectTab(n.dataset.tab);
    });
    $('#evo-field').addEventListener('change', () => { state.drawn.delete('evo'); drawTab('evo'); });
    $('#evo-view').addEventListener('change', () => { state.drawn.delete('evo'); drawTab('evo'); });
    $$('[data-export]').forEach((b) => b.addEventListener('click', () => exportData(b.dataset.export)));
    $('#print-report').addEventListener('click', printReport);
    $('#theme-toggle').addEventListener('click', toggleTheme);
    document.addEventListener('keydown', (e) => { if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') run(e); });

    applyPreset(PRESETS[0]?.id);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
