<?php
declare(strict_types=1);

/* Home page: what the laboratory solves and how it is validated. The laboratory itself is lab.php. */

require __DIR__ . '/src/bootstrap.php';

$cfg = manakov_config();
$h = static fn($s): string => htmlspecialchars((string)$s, ENT_QUOTES | ENT_SUBSTITUTE, 'UTF-8');
// a link that still carries the release placeholder (xxxx) is shown as pending instead of as a dead link
$ready = static fn(string $v): bool => $v !== '' && !str_contains($v, 'xxxx');
$repo = (string)$cfg['REPOSITORY_URL'];
$cite = (string)$cfg['CITATION_URL'];
$doi = (string)preg_replace('~^https?://(dx\.)?doi\.org/~i', '', (string)$cfg['DOI']);
$preview = @file_get_contents(__DIR__ . '/assets/img/collision-preview.svg');
header('X-Content-Type-Options: nosniff');
header('Referrer-Policy: same-origin');
?>
<!doctype html>
<html lang="en" data-theme="light">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Manakov Lab — vector solitons with variable coefficients and fractional derivatives</title>
  <meta name="description" content="Interactive laboratory for vector solitons of the variable-coefficient Manakov system and its Riesz space-fractional and Caputo time-fractional variants. Runs in the browser; the numerical core is pure PHP.">
  <script>
    // apply the saved theme before first paint (shared with the laboratory)
    try { const t = localStorage.getItem('manakov-theme'); if (t) document.documentElement.dataset.theme = t;
          else if (matchMedia('(prefers-color-scheme: dark)').matches) document.documentElement.dataset.theme = 'dark'; } catch (e) {}
  </script>
  <link rel="stylesheet" href="assets/vendor/katex/katex.min.css">
  <link rel="stylesheet" href="assets/css/app.css">
  <link rel="stylesheet" href="assets/css/home.css">
  <script defer src="assets/vendor/katex/katex.min.js"></script>
  <script defer src="assets/vendor/katex/contrib/auto-render.min.js"></script>
</head>
<body class="home">
<a class="skip" href="#model">Skip to the model</a>

<header class="topbar">
  <a class="brand" href="./">
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
    <a class="nav-section" href="#model">Model</a>
    <a class="nav-section" href="#validation">Validation</a>
    <a class="nav-section" href="#cite">Cite</a>
    <a class="btn primary small nav-cta" href="lab.php">Open the laboratory</a>
    <button id="theme-toggle" class="icon-btn" type="button" aria-label="Toggle dark mode" title="Toggle dark mode">
      <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>
    </button>
  </nav>
</header>

<main class="home-main">
  <!-- ================================================================== HERO -->
  <section class="hero">
    <div class="hero-text">
      <h1>Manakov Lab</h1>
      <p class="lead">An interactive laboratory for vector solitons of the Manakov system with variable coefficients
        and its fractional variants: set up an experiment in the browser, run it on the server and check the result
        against exact solutions.</p>
      <div class="hero-actions">
        <a class="btn primary" href="lab.php">Open the laboratory</a>
        <a class="btn ghost" href="#model">What it solves</a>
      </div>
      <p class="hero-note">Free and open source (MIT licence). Nothing to install: the numerical core is written in plain PHP.</p>
    </div>
    <figure class="hero-figure card">
      <?= $preview !== false ? $preview : '' ?>
      <figcaption>Two vector solitons collide and exchange energy between the components \(q_1\) and \(q_2\)
        (preset <em>Classical: shape-changing collision</em>, computed by the laboratory's own solver).</figcaption>
    </figure>
  </section>

  <!-- ================================================================= MODEL -->
  <section class="home-section" id="model" aria-labelledby="model-title">
    <h2 id="model-title">The model</h2>
    <p>The unknown is a two-component field \(\mathbf q=(q_1,q_2)\), for example two polarizations of light in an optical
      fibre or two species in a Bose&ndash;Einstein condensate. Both components share one cubic nonlinearity, which gives the
      system its vector (Manakov) solitons:</p>
    <div class="eq-block">\[ i\,\partial_t q_j+\beta(t)\,\partial_x^2 q_j+2\gamma(t)\big(|q_1|^2+|q_2|^2\big)q_j+V(x,t)\,q_j=i\,g(t)\,q_j,\qquad j=1,2. \]</div>
    <p>The coefficients are formulas chosen by the user: \(\beta(t)\) is the dispersion, \(\gamma(t)\) the strength of the
      nonlinearity, \(g(t)\) a gain (\(g&gt;0\)) or loss (\(g&lt;0\)), and \(V(x,t)\) an external potential, such as a trap.
      With constant coefficients and \(g=V=0\) the system is integrable; letting them vary in time models dispersion
      management, amplification and trapping.</p>

    <div class="variants">
      <article class="card variant">
        <h3>Classical</h3>
        <div class="variant-eq">\( i\,\partial_t q_j+\beta(t)\,\partial_x^2 q_j+\cdots \)</div>
        <p>The equation above. Solved by the Strang split-step Fourier method, second order in time.</p>
      </article>
      <article class="card variant">
        <h3>Space-fractional (Riesz)</h3>
        <div class="variant-eq">\( \partial_x^2 \;\longrightarrow\; -(-\partial_x^2)^{\alpha/2} \)</div>
        <p>Fractional dispersion of order \(\alpha\) (validated for \(1&lt;\alpha\le2\)): solitons with algebraic instead of
          exponential tails. Split-step with symbol \(|k|^\alpha\); solitons by Petviashvili iteration.</p>
      </article>
      <article class="card variant">
        <h3>Time-fractional (Caputo)</h3>
        <div class="variant-eq">\( i\,\partial_t \;\longrightarrow\; i\,{}^{C}\!D_t^{\nu} \)</div>
        <p>A memory of order \(0&lt;\nu\le1\) under which power is no longer conserved. L1 scheme on a graded mesh,
          order \(2-\nu\), with Anderson-accelerated nonlinear solves.</p>
      </article>
    </div>
  </section>

  <!-- ============================================================= FEATURES -->
  <section class="home-section" aria-labelledby="feat-title">
    <h2 id="feat-title">What you can do</h2>
    <ul class="features">
      <li><strong>Type the coefficients as formulas.</strong> \(\beta\), \(\gamma\) and \(g\) in \(t\), \(V\) in \(x\) and \(t\),
        checked and typeset as you write them; up to three initial solitons, each with its own amplitude, velocity, position
        and polarization.</li>
      <li><strong>Compare with exact solutions.</strong> Nonautonomous solitons, fractional solitons and manufactured Caputo
        solutions give the true error of a run, not only its plausibility.</li>
      <li><strong>Measure convergence automatically.</strong> On request, the run is repeated with smaller steps and the
        observed order is reported next to the theoretical one.</li>
      <li><strong>Export everything.</strong> Profiles, diagnostics and the space&ndash;time field as CSV, the full result and
        the problem setup as JSON, and a printable run report.</li>
    </ul>
    <p class="cta-row"><a class="btn primary" href="lab.php">Open the laboratory</a>
      <span class="muted">Eight ready-made experiments are included as starting points.</span></p>
  </section>

  <!-- =========================================================== VALIDATION -->
  <section class="home-section" id="validation" aria-labelledby="val-title">
    <h2 id="val-title">How it is validated</h2>
    <p>Every scheme is checked against a solution known in closed form: solitons carried through variable coefficients by a
      gauge reduction, fractional solitons computed to residuals near round-off, manufactured solutions of the Caputo model and
      a Mittag-Leffler mode of its linear part. The convergence studies test their own errors before publishing an order: an
      order is withheld when the error has reached the accuracy of the reference or has stopped decreasing, and the run says
      why. A second implementation in Python (NumPy), ported from the PHP core, reproduces the PHP powers and Hamiltonians to
      about \(10^{-12}\), and the scripts distributed with the code regenerate every reported number.</p>
  </section>

  <!-- ================================================================= CITE -->
  <section class="home-section cite" id="cite" aria-labelledby="cite-title">
    <h2 id="cite-title">Cite and source code</h2>
    <p>Manakov Lab is described in an article submitted to <em>SoftwareX</em>. Until it is published, please cite the
      software release:</p>
    <blockquote class="citation">A. H. Salas, C. D. Acosta M., E. Solorzano Tovar. <em>Manakov Lab</em>, version
      <?= $h(MANAKOV_LAB_VERSION) ?>. Zenodo<?php if ($ready($doi)): ?>,
      <a href="https://doi.org/<?= $h($doi) ?>">doi:<?= $h($doi) ?></a><?php endif; ?>.</blockquote>
    <ul class="links">
      <li><span class="links-k">Archived release (Zenodo)</span>
        <?php if ($ready($doi)): ?><a href="https://doi.org/<?= $h($doi) ?>">https://doi.org/<?= $h($doi) ?></a><?php else: ?><span class="muted">DOI assigned with the first public release</span><?php endif; ?></li>
      <li><span class="links-k">Source code (GitHub)</span>
        <?php if ($ready($repo)): ?><a href="<?= $h($repo) ?>"><?= $h(preg_replace('~^https?://~', '', $repo)) ?></a><?php else: ?><span class="muted">published with the first public release</span><?php endif; ?></li>
      <li><span class="links-k">How to cite</span>
        <?php if ($ready($cite)): ?><a href="<?= $h($cite) ?>"><?= str_contains($cite, 'doi.org') ? 'Zenodo record (citation export in BibTeX and other formats)' : 'citation file' ?></a><?php else: ?><span class="muted">citation file published with the first public release</span><?php endif; ?></li>
    </ul>
  </section>
</main>

<footer class="home-foot">
  <div><strong>Manakov Lab v<?= $h(MANAKOV_LAB_VERSION) ?></strong> · MIT licence ·
    Alvaro H. Salas, Carlos D. Acosta M., Elizabeth Solorzano Tovar</div>
  <div>Department of Mathematics and Statistics, Universidad Nacional de Colombia, Sede Manizales</div>
  <div class="muted">Math typeset with KaTeX · charts in the laboratory by Plotly</div>
</footer>

<script>
  document.addEventListener('DOMContentLoaded', function () {
    if (window.renderMathInElement) {
      renderMathInElement(document.body, { delimiters: [
        { left: '\\[', right: '\\]', display: true }, { left: '\\(', right: '\\)', display: false }], throwOnError: false });
    }
    document.getElementById('theme-toggle').addEventListener('click', function () {
      var root = document.documentElement;
      root.dataset.theme = root.dataset.theme === 'dark' ? 'light' : 'dark';
      try { localStorage.setItem('manakov-theme', root.dataset.theme); } catch (e) { /* storage unavailable */ }
    });
  });
</script>
</body>
</html>
