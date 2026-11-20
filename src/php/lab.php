<?php
declare(strict_types=1);

require __DIR__ . '/src/bootstrap.php';

use ManakovLab\Presets;

$cfg = manakov_config();
$presets = Presets::all();
$plotlySrc = is_file(__DIR__ . '/' . $cfg['PLOTLY_LOCAL']) ? $cfg['PLOTLY_LOCAL'] : $cfg['PLOTLY_CDN'];
$h = static fn($s): string => htmlspecialchars((string)$s, ENT_QUOTES | ENT_SUBSTITUTE, 'UTF-8');
$groups = [];
foreach ($presets as $p) {
    $groups[$p['group']][] = $p;
}
ini_set('serialize_precision', '-1');   // shortest round-trip floats in the embedded presets
$jsonFlags = JSON_HEX_TAG | JSON_HEX_AMP | JSON_HEX_APOS | JSON_HEX_QUOT | JSON_UNESCAPED_UNICODE;
header('X-Content-Type-Options: nosniff');
header('Referrer-Policy: same-origin');
?>
<!doctype html>
<html lang="en" data-theme="light">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Manakov Lab — laboratory</title>
  <meta name="description" content="Interactive numerical laboratory for the variable-coefficient Manakov system and its Riesz space-fractional and Caputo time-fractional variants.">
  <script>
    // apply the saved theme before first paint
    try { const t = localStorage.getItem('manakov-theme'); if (t) document.documentElement.dataset.theme = t;
          else if (matchMedia('(prefers-color-scheme: dark)').matches) document.documentElement.dataset.theme = 'dark'; } catch (e) {}
  </script>
  <link rel="stylesheet" href="assets/vendor/katex/katex.min.css">
  <link rel="stylesheet" href="assets/css/app.css">
  <script defer src="assets/vendor/katex/katex.min.js"></script>
  <script defer src="assets/vendor/katex/contrib/auto-render.min.js"></script>
  <script defer src="assets/vendor/katex/contrib/copy-tex.min.js"></script>
  <script defer src="<?= $h($plotlySrc) ?>" onerror="window.__plotlyFailed = true"></script>
  <script defer src="assets/js/app.js"></script>
</head>
<body>
<a class="skip" href="#results">Skip to results</a>

<header class="topbar">
  <a class="brand" href="./" title="Manakov Lab home">
    <svg class="brand-mark" viewBox="0 0 40 40" aria-hidden="true">
      <path d="M2 32 C 10 32, 12 6, 16 6 S 22 32, 30 32" fill="none" stroke="var(--series-1)" stroke-width="2.6" stroke-linecap="round"/>
      <path d="M10 32 C 18 32, 20 14, 24 14 S 30 32, 38 32" fill="none" stroke="var(--series-2)" stroke-width="2.6" stroke-linecap="round"/>
    </svg>
    <div>
      <div class="brand-name">Manakov Lab</div>
      <div class="brand-sub">Variable-coefficient · space-fractional · Caputo time-fractional</div>
    </div>
  </a>
  <nav class="topnav" aria-label="Sections">
    <a href="./">Home</a>
    <a href="#equation">Model</a>
    <a href="#results">Results</a>
    <a href="#theory">Theory &amp; methods</a>
    <button id="theme-toggle" class="icon-btn" type="button" aria-label="Toggle dark mode" title="Toggle dark mode">
      <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>
    </button>
  </nav>
</header>

<div id="plotly-banner" class="banner" hidden>
  The charting library could not be loaded, so charts are unavailable. Numerical results, the run report
  and the CSV/JSON exports remain available.
</div>

<main class="layout">
  <!-- ================================================================ CONTROLS -->
  <aside class="panel" aria-label="Simulation setup">
    <form id="setup" novalidate>
      <section class="group">
        <h2 class="group-title">Experiment</h2>
        <label class="field">
          <span class="label">Preset</span>
          <select id="preset">
            <?php foreach ($groups as $g => $items): ?>
              <optgroup label="<?= $h($g) ?>">
                <?php foreach ($items as $p): ?>
                  <option value="<?= $h($p['id']) ?>"><?= $h($p['title']) ?></option>
                <?php endforeach; ?>
              </optgroup>
            <?php endforeach; ?>
          </select>
        </label>
        <p id="preset-desc" class="hint"></p>
      </section>

      <section class="group">
        <h2 class="group-title">Model</h2>
        <div class="segmented" role="radiogroup" aria-label="Model">
          <button type="button" role="radio" data-model="classic" aria-checked="true">Classical</button>
          <button type="button" role="radio" data-model="space" aria-checked="false">Space-fractional</button>
          <button type="button" role="radio" data-model="caputo" aria-checked="false">Caputo (time)</button>
        </div>
      </section>

      <section class="group">
        <h2 class="group-title">Coefficients</h2>
        <p class="hint">Use <code>t</code> (and <code>x</code> in V), <code>pi</code>, <code>e</code>, + − * / ^ and
          sin, cos, tan, sinh, cosh, tanh, sech, exp, log, sqrt, abs, atan, sign, heaviside.
          Implicit products like <code>0.5cos(2t)</code> are allowed.</p>
        <?php foreach ([['beta', '\\beta(t)', 'dispersion', 't'], ['gamma', '\\gamma(t)', 'nonlinearity', 't'],
                        ['g', 'g(t)', 'gain (+) / loss (−)', 't'], ['V', 'V(x,t)', 'external potential', 'xt']] as [$id, $tex, $what, $vars]): ?>
          <div class="field coeff" data-coeff="<?= $id ?>">
            <label for="c-<?= $id ?>" class="label"><span class="tex" data-tex="<?= $h($tex) ?>"></span><span class="label-note"><?= $h($what) ?></span></label>
            <input id="c-<?= $id ?>" name="<?= $id ?>" type="text" spellcheck="false" autocomplete="off" data-vars="<?= $vars ?>">
            <div class="preview" aria-live="polite"></div>
          </div>
        <?php endforeach; ?>
        <p class="hint derived-note" hidden>In exact mode γ(t) is derived from the integrability condition
          \(\gamma=\beta\,e^{-2\int_0^t g}\) and the value above is ignored; V is switched off.</p>
      </section>

      <section class="group">
        <h2 class="group-title">Discretization</h2>
        <div class="grid2">
          <label class="field"><span class="label">Domain length \(L\)</span><input id="p-L" type="number" step="any" min="5" max="1000"></label>
          <label class="field"><span class="label">Grid points \(N_x\)</span>
            <select id="p-Nx"><option>128</option><option>256</option><option>512</option><option>1024</option><option>2048</option></select></label>
          <label class="field"><span class="label">Final time \(T\)</span><input id="p-T" type="number" step="any" min="0.01" max="200"></label>
          <label class="field" data-models="classic space"><span class="label">Time step \(\Delta t\)</span><input id="p-dt" type="number" step="any" min="0.00001" max="1"></label>
          <label class="field" data-models="space"><span class="label">Order \(\alpha\in[0.5,2]\), recommended \(1<\alpha\le2\)</span><input id="p-alpha" type="number" step="0.05" min="0.5" max="2"></label>
          <label class="field" data-models="caputo"><span class="label">Caputo order \(\nu\in(0,1]\)</span><input id="p-nu" type="number" step="0.05" min="0.05" max="1"></label>
          <label class="field" data-models="caputo"><span class="label">Time steps \(N\)</span><input id="p-N" type="number" step="1" min="16" max="3000"></label>
          <label class="field" data-models="caputo"><span class="label">Mesh grading \(r\)</span><input id="p-r" type="text" placeholder="auto = (2−ν)/ν"></label>
        </div>
        <p class="hint" id="cost-hint"></p>
      </section>

      <section class="group">
        <h2 class="group-title">Initial data <span class="muted">— bright vector solitons</span></h2>
        <p class="hint">Each soliton: amplitude \(2\eta\), velocity parameter \(\xi\) (phase \(e^{-2i\xi x}\)), centre \(x_0\),
          polarization \((c_1,c_2)\in\mathbb C^2\) (normalised).</p>
        <div id="solitons"></div>
        <button type="button" id="add-soliton" class="btn ghost small">+ Add soliton</button>
      </section>

      <section class="group">
        <h2 class="group-title">Validation</h2>
        <label class="check"><input type="checkbox" id="exact"> <span>Exact-solution mode <span class="muted">(error vs. exact solution)</span></span></label>
        <label class="check"><input type="checkbox" id="convergence"> <span>Convergence study <span class="muted">(4 resolutions)</span></span></label>
        <p class="hint" id="exact-hint"></p>
      </section>

      <div class="run-bar">
        <button type="submit" id="run" class="btn primary">
          <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M7 4.5v15l12-7.5z"/></svg>
          <span>Run simulation</span>
        </button>
        <div id="status" class="status" role="status" aria-live="polite"></div>
      </div>
    </form>
  </aside>

  <!-- ================================================================ CONTENT -->
  <div class="content">
    <section class="card" id="equation" aria-labelledby="eq-title">
      <div class="card-head">
        <h2 id="eq-title">Governing equation</h2>
        <span class="badge" id="scheme-badge"></span>
      </div>
      <div id="eq-main" class="eq-main"></div>
      <div id="eq-coeffs" class="eq-coeffs"></div>
    </section>

    <section class="card results" id="results" aria-labelledby="res-title">
      <div class="card-head">
        <h2 id="res-title">Results</h2>
        <div class="head-actions">
          <button type="button" class="btn ghost small" id="print-report" disabled>Print report</button>
        </div>
      </div>

      <div id="kpis" class="kpis" hidden></div>
      <div id="warnings" class="notes" hidden></div>
      <div id="alert" class="alert" role="alert" hidden></div>

      <div id="empty" class="empty">
        <svg viewBox="0 0 120 60" aria-hidden="true"><path d="M4 52 C 30 52, 36 8, 48 8 S 66 52, 92 52" fill="none" stroke="var(--series-1)" stroke-width="2"/><path d="M28 52 C 54 52, 60 22, 72 22 S 90 52, 116 52" fill="none" stroke="var(--series-2)" stroke-width="2"/></svg>
        <p>Choose a preset or set up your own problem, then press <strong>Run simulation</strong>.</p>
      </div>

      <div id="tabs-wrap" hidden>
        <div class="tabs" role="tablist" aria-label="Result views">
          <button role="tab" aria-selected="true" data-tab="evo">3D evolution</button>
          <button role="tab" aria-selected="false" data-tab="final">Final time</button>
          <button role="tab" aria-selected="false" data-tab="cons">Conservation</button>
          <button role="tab" aria-selected="false" data-tab="err">Errors &amp; convergence</button>
          <button role="tab" aria-selected="false" data-tab="data">Run data &amp; export</button>
        </div>

        <div class="tabpanel" data-panel="evo" role="tabpanel">
          <div class="toolbar">
            <label>Field
              <select id="evo-field">
                <option value="sum">|q₁|² + |q₂|²</option>
                <option value="I1">|q₁|²</option>
                <option value="I2">|q₂|²</option>
              </select></label>
            <label>View
              <select id="evo-view">
                <option value="surface">3D surface</option>
                <option value="heatmap">Space–time map</option>
              </select></label>
            <span class="toolbar-note">Drag to rotate · scroll to zoom · camera icon exports PNG</span>
          </div>
          <div id="plot-evo" class="plot tall"></div>
          <div class="pair">
            <figure><figcaption>|q₁(x,t)|² — space–time</figcaption><div id="plot-map1" class="plot short"></div></figure>
            <figure><figcaption>|q₂(x,t)|² — space–time</figcaption><div id="plot-map2" class="plot short"></div></figure>
          </div>
        </div>

        <div class="tabpanel" data-panel="final" role="tabpanel" hidden>
          <figure><figcaption id="cap-final">Moduli at the final time</figcaption><div id="plot-final" class="plot"></div></figure>
          <div class="pair">
            <figure><figcaption>Total intensity at the final time (log scale — reveals radiation)</figcaption><div id="plot-final-log" class="plot short"></div></figure>
            <figure><figcaption>Real and imaginary parts at the final time</figcaption><div id="plot-final-reim" class="plot short"></div></figure>
          </div>
          <div id="final-table"></div>
        </div>

        <div class="tabpanel" data-panel="cons" role="tabpanel" hidden>
          <div class="pair">
            <figure><figcaption>Component powers \(P_j(t)=\int|q_j|^2\,dx\)</figcaption><div id="plot-power" class="plot short"></div></figure>
            <figure><figcaption>Power balance \(P(t)\,/\,\big(P(0)\,e^{2\int_0^t g}\big)-1\) — zero when the power law holds</figcaption><div id="plot-balance" class="plot short"></div></figure>
            <figure><figcaption>Peak intensity \(\max_x\,(|q_1|^2+|q_2|^2)\)</figcaption><div id="plot-peak" class="plot short"></div></figure>
            <figure><figcaption>Hamiltonian, relative change \((H-H_0)/|H_0|\) — invariant only for constant \(\beta,\gamma\), \(g=0\), time-independent \(V\) and \(\nu=1\)</figcaption><div id="plot-ham" class="plot short"></div></figure>
          </div>
        </div>

        <div class="tabpanel" data-panel="err" role="tabpanel" hidden>
          <div id="err-empty" class="empty small" hidden>
            <p>Errors need a known exact solution. Enable <strong>Exact-solution mode</strong> (and optionally the
              <strong>Convergence study</strong>) and run again.</p>
          </div>
          <div id="err-content">
            <div class="pair">
              <figure><figcaption>Error \(\max_x|q-q_{\text{exact}}|\) versus time</figcaption><div id="plot-err" class="plot short"></div></figure>
              <figure><figcaption id="cap-conv">Convergence study</figcaption><div id="plot-conv" class="plot short"></div></figure>
            </div>
            <div id="conv-table"></div>
          </div>
        </div>

        <div class="tabpanel" data-panel="data" role="tabpanel" hidden>
          <div class="pair">
            <div><h3 class="sub">Run summary</h3><div id="run-summary"></div></div>
            <div><h3 class="sub">Coefficients used</h3><div id="run-coeffs"></div>
              <h3 class="sub">Export</h3>
              <div class="exports">
                <button type="button" class="btn ghost small" data-export="final">Final profile (CSV)</button>
                <button type="button" class="btn ghost small" data-export="diag">Diagnostics vs t (CSV)</button>
                <button type="button" class="btn ghost small" data-export="surface">Space–time field (CSV)</button>
                <button type="button" class="btn ghost small" data-export="json">Full result (JSON)</button>
                <button type="button" class="btn ghost small" data-export="config">Problem setup (JSON)</button>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div id="busy" class="busy" hidden>
        <div class="spinner" aria-hidden="true"></div>
        <div><strong>Solving…</strong> <span id="busy-time">0.0 s</span><div class="muted" id="busy-detail"></div></div>
      </div>
    </section>

    <section class="card theory" id="theory" aria-labelledby="th-title">
      <div class="card-head"><h2 id="th-title">Theory &amp; numerical methods</h2></div>

      <details open>
        <summary>The three models</summary>
        <p>All three act on a two-component field \(\mathbf q=(q_1,q_2)^{\mathsf T}\) with the Manakov (SU(2)-symmetric) cubic nonlinearity, \(j=1,2\):</p>
        <p><strong>Variable-coefficient Manakov</strong></p>
        \[ i\,\partial_t q_j+\beta(t)\,\partial_x^2 q_j+2\gamma(t)\big(|q_1|^2+|q_2|^2\big)q_j+V(x,t)\,q_j=i\,g(t)\,q_j \]
        <p><strong>Space-fractional (Riesz) Manakov</strong>, \(0.5\le\alpha\le 2\) (tested at \(\alpha=1.2,\,1.5,\,2\) and recommended for \(1<\alpha\le2\); \(\alpha=1\) is \(L^2\)-critical and \(\alpha<1\) supercritical, so solitons may be unstable and every run with \(\alpha\le1\) carries a warning), where \((-\partial_x^2)^{\alpha/2}\) has Fourier symbol \(|k|^\alpha\):</p>
        \[ i\,\partial_t q_j-\beta(t)\,(-\partial_x^2)^{\alpha/2} q_j+2\gamma(t)\big(|q_1|^2+|q_2|^2\big)q_j+V(x,t)\,q_j=i\,g(t)\,q_j \]
        <p><strong>Time-fractional (Caputo) Manakov</strong>, \(0<\nu\le 1\), with classical dispersion (\(\alpha=2\)):</p>
        \[ i\,{}^{C}\!D_t^{\nu} q_j+\beta(t)\,\partial_x^2 q_j+2\gamma(t)\big(|q_1|^2+|q_2|^2\big)q_j+V(x,t)\,q_j=i\,g(t)\,q_j,
           \qquad {}^{C}\!D_t^{\nu}u(t)=\frac{1}{\Gamma(1-\nu)}\int_0^t\frac{u'(s)}{(t-s)^{\nu}}\,ds \]
        <p>For \(\alpha=2\) and \(\nu=1\) the three models coincide.</p>
      </details>

      <details>
        <summary>Split-step Fourier method (classical and space-fractional)</summary>
        <p>On a periodic grid the equation is split into a linear part, solved exactly in Fourier space,</p>
        \[ \hat q_j(k,t_1)=\hat q_j(k,t_0)\,\exp\!\Big(-i\,|k|^{\alpha}\!\int_{t_0}^{t_1}\!\beta(s)\,ds\Big), \]
        <p>and a nonlinear part, solved exactly pointwise because \(|\mathbf q|^2\) evolves only through the gain:</p>
        \[ q_j(x,t_1)=q_j(x,t_0)\,e^{\int_{t_0}^{t_1}g}\exp\!\Big(i\!\int_{t_0}^{t_1}\!\big[2\gamma(s)\,|\mathbf q(x,t_0)|^2e^{2\int_{t_0}^{s}g}+V(x,s)\big]ds\Big). \]
        <p>The time integrals use Simpson's rule on half steps. If \(\Delta t\) does not resolve a coefficient, the run follows the discretized
          coefficient (for instance, the power follows \(e^{2\sum g}\) rather than \(e^{2\int g}\)); each run compares the quadrature it applies
          (Simpson step sums for \(\beta\) and \(g\), one gain-weighted Simpson panel per step for the nonlinear phase, one Simpson panel per
          step for a time-dependent \(V\), checked at 16 grid points and the maxima of \(|V|\)) with a converged reference
          integral, itself cross-checked with Gauss&ndash;Legendre, and warns when they differ by more than \(10^{-6}\). These checks are estimates:
          a coefficient that oscillates in phase with every sampling grid can escape them. Strang composition \(\mathcal L_{\Delta t/2}\,\mathcal N_{\Delta t}\,\mathcal L_{\Delta t/2}\) is
          second-order. Both sub-steps conserve \(\int|\mathbf q|^2e^{-2\int g}\) exactly, so the discrete solution stays bounded in \(L^2\)
          for every \(\Delta t\); this does not guarantee accuracy for large \(\Delta t\), where split-step methods for the NLS can show
          instabilities (Weideman &amp; Herbst, 1986).</p>
      </details>

      <details>
        <summary>L1 scheme on a graded mesh (Caputo)</summary>
        <p>On the mesh \(t_n=T\,(n/N)^r\), the Caputo derivative is approximated by</p>
        \[ {}^{C}\!D_t^{\nu}u(t_n)\approx\sum_{k=1}^{n}a_{n,k}\,(u^k-u^{k-1}),\qquad
           a_{n,k}=\frac{(t_n-t_{k-1})^{1-\nu}-(t_n-t_k)^{1-\nu}}{\Gamma(2-\nu)\,(t_k-t_{k-1})}. \]
        <p>Caputo solutions behave like \(c_0+c_1t^{\nu}\) near \(t=0\), so a uniform mesh gives only order \(\nu\); the grading
          \(r=(2-\nu)/\nu\) restores the optimal global order \(2-\nu\) (Stynes, O'Riordan &amp; Gracia, 2017).
          Dispersion and gain are treated implicitly in Fourier space; the nonlinear, potential and forcing terms implicitly, by fixed-point
          (Picard) iteration accelerated with Anderson mixing, which converges also where plain Picard is not a contraction (small \(\nu\),
          larger amplitudes). The full memory sum is kept, so the work grows like \(N^2N_x\) and the memory like \(N\,N_x\). At \(\nu=1\) the L1 scheme would reduce to the
          dissipative backward Euler method, so runs with \(\nu=1\) use the split-step solver; close to \(\nu=1\) the L1 scheme still adds
          numerical dissipation and \(P(T)\) should be checked under refinement of \(N\).</p>
      </details>

      <details>
        <summary>Exact solutions used for validation</summary>
        <p><strong>Gauge reduction</strong> (classical and space-fractional, \(V=0\)). If \(\gamma(t)=\beta(t)\,e^{-2G(t)}\) with \(G=\int_0^tg\), then</p>
        \[ \mathbf q(x,t)=e^{G(t)}\,\mathbf p\big(x,\tau(t)\big),\qquad \tau(t)=\int_0^t\beta(s)\,ds, \]
        <p>where \(\mathbf p\) solves the constant-coefficient system. For \(\alpha=2\) this gives the Manakov bright soliton</p>
        \[ \mathbf p=\mathbf c\,2\eta\,\operatorname{sech}\!\big(2\eta(x-x_0+4\xi\tau)\big)\,e^{-2i\xi x-4i(\xi^2-\eta^2)\tau},\qquad|\mathbf c|=1; \]
        <p>for \(\alpha<2\) the stationary fractional soliton \(\mathbf p=\mathbf c\,\phi(x)\,e^{i\mu\tau}\), \(\big(\mu+(-\partial_x^2)^{\alpha/2}\big)\phi=2\phi^3\),
          computed by Petviashvili iteration (\(\mu=4\eta^2\)); a warning reports non-convergence, or a box too small for a localized state
          (homogeneous or periodic stationary branch), and the residual of the profile actually used enters the reference floor.</p>
        <p><strong>Periodic box.</strong> The classical reference is the soliton periodized over the box, \(\sum_m\mathbf p(x+mL,\tau)\), so a soliton
          may cross \(x=\pm L/2\). Its interaction with the images, the Fourier truncation on the grid and the profile residual form a
          <em>reference floor</em>; observed orders of errors below ten times that floor (or near round-off) are not reported. The image and
          spatial parts are estimates, so the convergence study also fits \(e=Ch^p+S\) to three consecutive levels and withholds orders
          whose errors stop decreasing as a power of \(\Delta t\) or exceed 10&nbsp;% of the peak; a study with only two levels publishes no order, and the error is measured at every step when the output times under-sample a coefficient.</p>
        <p><strong>Manufactured solution</strong> (Caputo). \(\mathbf q^*=\mathbf c\,(1+t^{\nu}+i t^2)\,\psi(x)\) with
          \({}^{C}\!D_t^{\nu}\big(1+t^\nu+it^2\big)=\Gamma(1+\nu)+\frac{2i\,t^{2-\nu}}{\Gamma(3-\nu)}\); the forcing that makes \(\mathbf q^*\) exact is computed
          automatically, for any \(\beta,\gamma,g,V\).</p>
      </details>

      <details>
        <summary>References</summary>
        <ul class="refs">
          <li>S. V. Manakov, <em>Sov. Phys. JETP</em> 38 (1974) 248.</li>
          <li>M. J. Ablowitz, B. Prinari, A. D. Trubatch, <em>Discrete and Continuous Nonlinear Schrödinger Systems</em>, Cambridge (2004).</li>
          <li>R. Radhakrishnan, M. Lakshmanan, J. Hietarinta, <em>Phys. Rev. E</em> 56 (1997) 2213.</li>
          <li>N. Laskin, Fractional Schrödinger equation, <em>Phys. Rev. E</em> 66 (2002) 056108.</li>
          <li>M. Stynes, E. O'Riordan, J. L. Gracia, <em>SIAM J. Numer. Anal.</em> 55 (2017) 1057.</li>
          <li>V. I. Petviashvili, <em>Sov. J. Plasma Phys.</em> 2 (1976) 257.</li>
          <li>D. E. Pelinovsky, Y. A. Stepanyants, <em>SIAM J. Numer. Anal.</em> 42 (2004) 1110.</li>
          <li>G. Strang, <em>SIAM J. Numer. Anal.</em> 5 (1968) 506.</li>
          <li>J. A. C. Weideman, B. M. Herbst, <em>SIAM J. Numer. Anal.</em> 23 (1986) 485.</li>
          <li>Y. Lin, C. Xu, <em>J. Comput. Phys.</em> 225 (2007) 1533.</li>
          <li>V. N. Serkin, A. Hasegawa, T. L. Belyaeva, <em>Phys. Rev. Lett.</em> 98 (2007) 074102.</li>
          <li>C. Klein, C. Sparber, P. Markowich, <em>Proc. R. Soc. A</em> 470 (2014) 20140364.</li>
          <li>D. G. Anderson, <em>J. ACM</em> 12 (1965) 547; H. F. Walker, P. Ni, <em>SIAM J. Numer. Anal.</em> 49 (2011) 1715.</li>
        </ul>
      </details>
    </section>

    <footer class="foot">Manakov Lab v<?= $h(MANAKOV_LAB_VERSION) ?> · MIT licence ·
      <a href="<?= $h($cfg['REPOSITORY_URL']) ?>">source code</a> ·
      <a href="<?= $h($cfg['CITATION_URL']) ?>">how to cite</a> ·
      numerical core in pure PHP · math typeset with KaTeX · charts by Plotly</footer>
  </div>
</main>

<template id="soliton-row">
  <div class="soliton">
    <div class="soliton-head"><span class="soliton-title"></span><button type="button" class="link remove">Remove</button></div>
    <div class="grid3">
      <label class="field"><span class="label">\(\eta\)</span><input data-k="eta" type="number" step="any" min="0.01" max="5"></label>
      <label class="field"><span class="label">\(\xi\)</span><input data-k="xi" type="number" step="any" min="-5" max="5"></label>
      <label class="field"><span class="label">\(x_0\)</span><input data-k="x0" type="number" step="any"></label>
    </div>
    <div class="grid4">
      <label class="field"><span class="label">\(\operatorname{Re}c_1\)</span><input data-k="p1re" type="number" step="any"></label>
      <label class="field"><span class="label">\(\operatorname{Im}c_1\)</span><input data-k="p1im" type="number" step="any"></label>
      <label class="field"><span class="label">\(\operatorname{Re}c_2\)</span><input data-k="p2re" type="number" step="any"></label>
      <label class="field"><span class="label">\(\operatorname{Im}c_2\)</span><input data-k="p2im" type="number" step="any"></label>
    </div>
  </div>
</template>

<script>window.MANAKOV_PRESETS = <?= json_encode($presets, $jsonFlags) ?>;</script>
</body>
</html>
