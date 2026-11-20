<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

/**
 * Sets up and runs one experiment: initial data, exact solutions (validation mode),
 * diagnostics (powers, peak, Hamiltonian, error), convergence study and the JSON-ready result.
 *
 * Input $P (already validated by the API layer):
 *   model 'classic'|'space'|'caputo', L, Nx, T, dt (split-step) or N, r (Caputo), alpha, nu,
 *   ast => [beta, gamma, g, V] (ExprParser ASTs), solitons => list of
 *   [eta, xi, x0, p1re, p1im, p2re, p2im], exact (bool), convergence (bool)
 */
final class Simulation
{
    private const N_SAVE = 400;

    private array $P;
    private array $x;
    private array $k;
    private float $dx;
    private \Closure $beta;
    private \Closure $gamma;
    private \Closure $g;
    private ?\Closure $V;
    private \Closure $G;
    private \Closure $tau;
    private array $warnings = [];
    private ?float $deadline;
    /**
     * Exact mode: estimated size of the reference's own defects (absolute, in |q|), by source:
     * images (interaction of the periodized classical soliton with its periodic images), spatial (Fourier
     * truncation of the classical soliton on the grid), profile (residual of the fractional profile).
     * Errors below ten times their sum are not attributed to the time integrator.
     */
    private array $floor = ['images' => 0.0, 'spatial' => 0.0, 'profile' => 0.0, 'roundoff' => 0.0];
    /** Largest relative residual of the fractional profiles actually used (null if none). */
    private ?float $profileResidual = null;
    /** Description of the fractional stationary solution used (branch, edge/peak ratio, iterations). */
    private ?array $profileInfo = null;
    private array $Ginfo;
    private array $tauInfo;
    /** Antiderivative of the effective gamma (exact mode overrides it), built on demand. */
    private ?\Closure $Ggamma = null;
    /** Time-dependent V(x,t) evaluated at the sample points of its quadrature diagnostic (null otherwise). */
    private ?\Closure $Vsample = null;
    /** Wall-clock time at construction, the origin of runtime_total. */
    private float $tCreated;
    /** Number of grid points at which the quadrature of V(x,t) is checked. */
    private int $Vpoints = 0;
    /** @var array<int, array<string, float>> step-quadrature discrepancies by number of steps */
    private array $quadCache = [];

    /** Relative overlap with the periodic images above which the classical exact reference is refused. */
    private const REF_REFUSE = 1e-3;
    /** Absolute reference defect below which nothing is reported. */
    private const FLOOR_REPORT = 1e-10;
    /** Relative residual above which a fractional profile is not treated as stationary. */
    private const PROFILE_TOL = 1e-9;
    /** Edge/peak ratio above which a fractional stationary solution is reported as not localized. */
    private const LOCALIZED = 1e-2;
    /** Relative discrepancy of the split-step quadrature of a coefficient above which dt does not resolve it. */
    private const QUAD_TOL = 1e-6;
    /** Dimensionless step-independent fraction of the error above which an observed order is withheld. */
    private const SATURATION = 0.1;
    /** Error, relative to the peak modulus, above which an observed order is not in the asymptotic regime. */
    private const PREASYMPTOTIC = 0.1;
    /** Relative interpolation defect of a coefficient between output times above which outputs under-sample it. */
    private const OUTPUT_ALIAS = 0.05;
    /** Rounding level of an exact-mode comparison, relative to the peak modulus of the initial field. */
    private const ROUNDOFF = 1e-12;

    public function __construct(array $P, ?float $deadline = null)
    {
        $this->tCreated = microtime(true);   // runtime_total includes the reference tables built here
        $this->P = $P;
        $this->deadline = $deadline;
        $this->warnings = array_values($P['warnings'] ?? []);
        $this->domainChecks();
        $Nx = $P['Nx'];
        $this->dx = $P['L'] / $Nx;
        $x = [];
        for ($i = 0; $i < $Nx; $i++) {
            $x[] = ($i - intdiv($Nx, 2)) * $this->dx;
        }
        $this->x = $x;
        $this->k = FFT::wavenumbers($Nx, $P['L']);
        $this->beta = ExprCompiler::timeFunction($P['ast']['beta']);
        $this->gamma = ExprCompiler::timeFunction($P['ast']['gamma']);
        $this->g = ExprCompiler::timeFunction($P['ast']['g']);
        $this->V = ExprCompiler::fieldFunction($P['ast']['V'], $x);
        $this->checkCoefficients();
        if ($this->V !== null && in_array('t', ExprCompiler::variables($P['ast']['V']), true)) {
            // 16 equally spaced grid points and the largest |V| at t = 0 and t = T
            $idx = [];
            for ($j = 0; $j < min(16, $Nx); $j++) {
                $idx[] = intdiv($j * $Nx, min(16, $Nx));
            }
            foreach ([0.0, (float)$P['T']] as $t) {
                $abs = array_map('abs', ($this->V)($x, $t));
                $idx[] = (int)array_search(max($abs), $abs, true);
            }
            $idx = array_values(array_unique($idx));
            $this->Vpoints = count($idx);
            $xs = array_map(static fn($i) => $x[$i], $idx);
            $Vs = ExprCompiler::fieldFunction($P['ast']['V'], $xs);
            $this->Vsample = static fn(float $t): array => $Vs($xs, $t);
        }
        [$this->G, $this->Ginfo] = Integral::antiderivativeAdaptive($this->g, $P['T']);
        [$this->tau, $this->tauInfo] = Integral::antiderivativeAdaptive($this->beta, $P['T']);
        foreach (['g' => $this->Ginfo, 'beta' => $this->tauInfo] as $name => $info) {
            if (!$info['resolved']) {
                $this->warnings[] = "$name(t) varies too fast on [0, T] for the reference quadrature ({$info['m']} Simpson panels): "
                    . "its integral, used by the exact solutions and the gain diagnostic, may be inaccurate.";
            }
        }
    }

    /**
     * Validated domain of the models (see the theory panel of lab.php): warnings and refusals that do not
     * depend on the coefficients.
     */
    private function domainChecks(): void
    {
        $P = $this->P;
        $this->mathematicalDomain();
        if ($P['model'] === 'space' && $P['alpha'] <= 1.0) {
            if ($P['exact']) {
                throw new SimulationException('Exact mode is available only for 1 < alpha <= 2: for alpha <= 1 the cubic '
                    . 'fractional NLS is L2-critical (alpha = 1) or supercritical, and fractional solitons may be unstable.', 'alpha');
            }
            $this->warnings[] = 'alpha <= 1: the cubic fractional NLS is L2-critical (alpha = 1) or supercritical (alpha < 1); '
                . 'solitons may be unstable and solutions may collapse. Results are validated only for 1 < alpha <= 2.';
        }
        if ($P['model'] === 'caputo') {
            if ($P['nu'] >= 1.0) {
                if ($P['exact']) {
                    throw new SimulationException('Caputo exact mode requires nu < 1. For nu = 1 use the Classical model.', 'nu');
                }
                $this->warnings[] = 'nu = 1: the Caputo derivative is the ordinary derivative, so the run uses the conservative '
                    . 'split-step solver with dt = T/N (the L1 scheme would reduce to the dissipative backward Euler method).';
                if (isset($P['dt']) && abs($P['dt'] - $P['T'] / $P['N']) > 1e-12 * $P['T']) {
                    $this->warnings[] = sprintf('nu = 1: the step is set by N; dt = %.6g was ignored and dt = T/N = %.6g used.',
                        $P['dt'], $P['T'] / $P['N']);
                }
            } elseif ($P['nu'] >= 0.95) {
                $this->warnings[] = 'nu close to 1: the L1 scheme is close to backward Euler and adds numerical dissipation; '
                    . 'check P(T)/P(0) under refinement of N before interpreting the power loss.';
            }
        }
    }

    /**
     * Mathematical domain of every model, enforced in the core so that library use (tests, scripts and the
     * Python port) obeys the same contract as the web API. The API adds the work-budget limits on top.
     */
    private function mathematicalDomain(): void
    {
        $P = $this->P;
        $bad = static function (string $msg, string $field): never {
            throw new SimulationException($msg, $field);
        };
        $pos = static fn($v): bool => is_float($v) || is_int($v) ? is_finite((float)$v) && $v > 0 : false;
        if (!in_array($P['model'] ?? null, ['classic', 'space', 'caputo'], true)) {
            $bad("model must be 'classic', 'space' or 'caputo'.", 'model');
        }
        if (!$pos($P['L'] ?? null)) {
            $bad('L must be a positive finite number.', 'L');
        }
        $Nx = $P['Nx'] ?? 0;
        if (!is_int($Nx) || $Nx < 4 || ($Nx & ($Nx - 1)) !== 0) {
            $bad('Nx must be a power of two, at least 4 (radix-2 FFT).', 'Nx');
        }
        if (!$pos($P['T'] ?? null)) {
            $bad('T must be a positive finite number.', 'T');
        }
        if ($P['model'] === 'caputo') {
            $nu = $P['nu'] ?? null;
            if (!(is_float($nu) || is_int($nu)) || !is_finite((float)$nu) || $nu <= 0.0 || $nu > 1.0) {
                $bad('The Caputo order must satisfy 0 < nu <= 1.', 'nu');
            }
            if (!is_int($P['N'] ?? null) || $P['N'] < 1) {
                $bad('N (number of time steps) must be a positive integer.', 'N');
            }
            if ($nu < 1.0 && (!$pos($P['r'] ?? null) || $P['r'] < 1.0)) {
                $bad('The grading exponent must satisfy r >= 1.', 'r');
            }
        } else {
            $dt = $P['dt'] ?? null;
            if (!$pos($dt) || $dt > $P['T']) {
                $bad('dt must satisfy 0 < dt <= T.', 'dt');
            }
        }
        if ($P['model'] === 'space') {
            $a = $P['alpha'] ?? null;
            if (!(is_float($a) || is_int($a)) || !is_finite((float)$a) || $a <= 0.0 || $a > 2.0) {
                $bad('The fractional order must satisfy 0 < alpha <= 2.', 'alpha');
            }
        }
        if (empty($P['solitons']) || !is_array($P['solitons'])) {
            $bad('At least one soliton is required.', 'solitons');
        }
        foreach ($P['solitons'] as $s) {
            if (!$pos($s['eta'] ?? null)) {
                $bad('Every soliton needs eta > 0.', 'eta');
            }
            foreach (['xi', 'x0', 'p1re', 'p1im', 'p2re', 'p2im'] as $f) {
                if (!is_finite((float)($s[$f] ?? NAN))) {
                    $bad("Soliton parameter $f must be a finite number.", $f);
                }
            }
        }
    }

    private function checkCoefficients(): void
    {
        $T = $this->P['T'];
        foreach (['beta' => $this->beta, 'gamma' => $this->gamma, 'g' => $this->g] as $name => $f) {
            for ($i = 0; $i <= 40; $i++) {
                if (!is_finite($f($T * $i / 40))) {
                    throw new SimulationException("$name(t) is not finite on [0, T].", $name);
                }
            }
        }
        if ($this->V !== null) {
            for ($i = 0; $i <= 10; $i++) {
                foreach (($this->V)($this->x, $T * $i / 10) as $v) {
                    if (!is_finite($v)) {
                        throw new SimulationException('V(x,t) is not finite on the domain.', 'V');
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------ setup
    private function polarization(array $s): array
    {
        return Solitons::normalise([(float)$s['p1re'], (float)$s['p1im'], (float)$s['p2re'], (float)$s['p2im']]);
    }

    /**
     * Fractional profile centred at x0. Petviashvili computes it centred at a node; it is described (localized
     * soliton, periodic or homogeneous stationary state of the box, convergence) and its residual after the
     * shift is measured, since a shift by a non-integer number of cells does not commute with the cubic term.
     * In exact mode ($snap) an off-node centre with a residual above PROFILE_TOL is moved to the nearest node;
     * a residual still above PROFILE_TOL afterwards becomes a reference defect (floor['profile']).
     */
    private function fractionalProfile(float $eta, float $x0, bool $snap = false): array
    {
        $alpha = $this->P['alpha'];
        $mu = 4.0 * $eta * $eta;
        $maxit = 3000;
        [$phi, , $iters] = Solitons::fractional($this->x, $this->k, $alpha, $mu, 1e-13, $maxit);
        $this->describeProfile($phi, $iters > $maxit ? null : $iters, $eta);
        $prof = Solitons::shift($phi, $this->k, $x0);
        $rel = Solitons::relativeResidual($prof, $this->k, $alpha, $mu);
        $m = round($x0 / $this->dx);
        if ($rel > self::PROFILE_TOL && abs($x0 - $m * $this->dx) > 1e-12 * max(1.0, abs($x0))) {
            if ($snap) {
                $x0n = $m * $this->dx;
                $this->warnings[] = sprintf('x0 = %.6g lies between grid nodes: the Fourier-shifted fractional profile is not '
                    . 'a stationary solution on this grid (relative residual %.1e), and exact mode would report that defect '
                    . 'as time-integration error. The soliton was moved to the nearest node, x0 = %.6g.', $x0, $rel, $x0n);
                $prof = Solitons::shift($phi, $this->k, $x0n);
                $rel = Solitons::relativeResidual($prof, $this->k, $alpha, $mu);
            } else {
                $this->warnings[] = sprintf('The fractional soliton at x0 = %.6g lies between grid nodes and its shifted profile '
                    . 'is not stationary on this grid (relative residual %.1e). Increase Nx or place x0 on a node '
                    . '(a multiple of dx = %.6g).', $x0, $rel, $this->dx);
            }
        }
        if ($snap && $rel > self::PROFILE_TOL) {
            // the defect r of the stationary equation acts as a forcing: |error| <~ |r| |tau| e^G
            $c3 = 0.0;
            foreach ($prof as $v) {
                $c3 = max($c3, abs(2.0 * $v * $v * $v));
            }
            $this->floor['profile'] = $rel * $c3 * max(1.0, $this->tauInfo['maxAbs']) * exp($this->Ginfo['max']);
            $this->warnings[] = sprintf('Exact mode: the fractional profile satisfies its stationary equation only to a relative '
                . 'residual %.1e (Petviashvili %s). This defect alone can produce errors of about %.1e, which are not '
                . 'attributed to the time integrator: the corresponding observed orders are omitted.',
                $rel, $this->profileInfo['iterations'] === null ? 'did not converge' : 'converged to its tolerance',
                $this->floor['profile']);
        }
        $this->profileResidual = max($this->profileResidual ?? 0.0, $rel);
        return $prof;
    }

    /** Identifies the stationary solution computed on the box and warns if it is not a localized soliton. */
    private function describeProfile(array $phi, ?int $iters, float $eta): void
    {
        $mx = max(array_map('abs', $phi));
        $mn = min(array_map('abs', $phi));
        $edge = $mx > 0.0 ? abs($phi[0]) / $mx : 1.0;          // phi[0] is at x = -L/2
        $branch = ($mx > 0.0 && $mn / $mx > 1.0 - 1e-6) ? 'homogeneous' : ($edge > self::LOCALIZED ? 'periodic' : 'localized');
        $this->profileInfo = ['branch' => $branch, 'edge' => $edge, 'iterations' => $iters, 'amplitude' => $mx];
        if ($iters === null) {
            $this->warnings[] = 'Petviashvili iteration did not converge in 3000 iterations: the fractional profile is not a '
                . 'converged stationary solution on this grid.';
        }
        if ($branch === 'homogeneous') {
            $this->warnings[] = sprintf('The stationary solution found on this box is spatially constant, phi = sqrt(mu/2) = %.6g '
                . '(the homogeneous branch), not a soliton: the box is too small for a localized state with eta = %.4g. '
                . 'Exact mode then tests the time integration of a plane wave. Increase L.', sqrt(2.0 * $eta * $eta), $eta);
        } elseif ($branch === 'periodic') {
            $this->warnings[] = sprintf('The stationary solution found on this box is not localized: its value at x = +-L/2 is '
                . '%.2g of the peak, so it is a periodic stationary state of the box rather than an isolated soliton. '
                . 'Increase L for a localized soliton.', $edge);
        }
    }

    /**
     * Classical exact mode: the reference is the soliton periodized over the box, S = sum_m p(x + mL, tau). S solves
     * the periodic linear problem exactly; in the cubic term the cross products with the images leave a defect whose
     * leading part near the soliton is 6 |p|^2 |p_{-1} + p_{+1}| ~ 96 eta^3 sech(2 eta L). floor['images'] scales it by
     * max(1, |tau|): an estimate of the interaction, not a strict bound. Soliton and images must remain separated:
     * the reference is refused when the tail midway between them, 2 sech(eta L), exceeds REF_REFUSE. The Fourier
     * truncation of the initial soliton on the grid is estimated from its spectrum, pi sech(pi (k + 2 xi) / (4 eta))
     * (floor['spatial']); it does not grow with time, so convergence() also tests the errors a posteriori.
     */
    private function checkClassicalReference(float $eta, float $xi): void
    {
        $L = $this->P['L'];
        $eG = exp($this->Ginfo['max']);
        $overlap = 2.0 / cosh($eta * $L);
        if ($overlap > self::REF_REFUSE) {
            throw new SimulationException(sprintf('Exact mode: the soliton (eta = %.4g) is too wide for the periodic box: midway '
                . 'to its periodic image its tail is %.1e of the peak, so it interacts with its images and is no longer the '
                . 'reference solution. Use L >= %.4g.', $eta, $overlap, ceil(log(4.0 / self::REF_REFUSE) / $eta)), 'L');
        }
        $growth = 96.0 * $eta ** 3 * max(1.0, $this->tauInfo['maxAbs']) * $eG;
        $this->floor['images'] = $growth / cosh(2.0 * $eta * $L);
        $kN = M_PI * $this->P['Nx'] / $L;
        $d = $kN - 2.0 * abs($xi);
        if ($d <= 0.0) {
            throw new SimulationException(sprintf('Exact mode: the carrier wavenumber 2|xi| = %.4g exceeds the largest grid '
                . 'wavenumber pi Nx/L = %.4g. Increase Nx.', 2.0 * abs($xi), $kN), 'Nx');
        }
        $this->floor['spatial'] = 4.0 * $eta / M_PI * $eG * exp(-M_PI * $d / (4.0 * $eta));
        if ($this->floor['images'] > self::FLOOR_REPORT) {
            $this->warnings[] = sprintf('Exact mode: the soliton interacts weakly with its periodic images (tail %.1e of the peak '
                . 'midway); the periodized reference is estimated to be accurate to about %.1e (L >= %.4g makes this negligible). Observed '
                . 'orders of errors below ten times this value are omitted.', $overlap, $this->floor['images'],
                ceil(log(2.0 * $growth / self::FLOOR_REPORT) / (2.0 * $eta)));
        }
        if ($this->floor['spatial'] > self::FLOOR_REPORT) {
            $need = 2.0 * abs($xi) + 4.0 * $eta / M_PI * log(4.0 * $eta / M_PI * $eG / self::FLOOR_REPORT);
            $nxMin = 2 ** (int)ceil(log($need * $L / M_PI, 2));
            $this->warnings[] = sprintf('Exact mode: the grid does not resolve the soliton spectrum (dx = %.4g for a soliton of '
                . 'width 1/(2 eta) = %.4g). The truncation of the initial soliton alone is estimated at %.1e; the spatial error '
                . 'can grow with time beyond this estimate. Observed orders of errors below ten times it are omitted, and the '
                . 'convergence study also omits orders whose errors stop decreasing with dt. Nx >= %d makes it negligible.',
                $this->dx, 0.5 / $eta, $this->floor['spatial'], $nxMin);
        }
    }

    private function initialCondition(): array
    {
        $Q = Solitons::zeros($this->P['Nx']);
        foreach ($this->P['solitons'] as $s) {
            $c = $this->polarization($s);
            if ($this->P['model'] === 'space' && $this->P['alpha'] < 2.0) {
                $prof = $this->fractionalProfile((float)$s['eta'], (float)$s['x0']);
                $ph = array_map(static fn($xv) => -2.0 * $s['xi'] * $xv, $this->x);
                Solitons::addScaled($Q, $c, $prof, $ph);
            } else {
                Solitons::addManakov($Q, $this->x, 0.0, (float)$s['eta'], (float)$s['xi'], $c, (float)$s['x0']);
            }
        }
        return $Q;
    }

    /** @return array{0: \Closure, 1: ?\Closure, 2: ?\Closure} exact(t), gamma override, forcing */
    private function setupExact(): array
    {
        $P = $this->P;
        if (count($P['solitons']) > 1) {
            $this->warnings[] = 'Exact mode uses only the first soliton.';
        }
        $s = $P['solitons'][0];
        $eta = (float)$s['eta']; $xi = (float)$s['xi']; $x0 = (float)$s['x0'];
        $c = $this->polarization($s);
        $x = $this->x; $G = $this->G; $tau = $this->tau; $beta = $this->beta;

        if ($P['model'] !== 'caputo') {
            if ($this->V !== null) {
                $this->warnings[] = 'Exact mode requires V = 0; the potential was switched off.';
                $this->V = null;
            }
            $gammaEx = static fn(float $t): float => $beta($t) * exp(-2.0 * $G($t));
            if ($P['model'] === 'space' && $P['alpha'] < 2.0) {
                if ($xi != 0.0) {
                    $this->warnings[] = 'For alpha < 2 the exact soliton is stationary: xi was set to 0.';
                }
                $mu = 4.0 * $eta * $eta;
                $prof = $this->fractionalProfile($eta, $x0, true);
                $n = count($x);
                $exact = static function (float $t) use ($prof, $c, $mu, $G, $tau, $n): array {
                    $Q = Solitons::zeros($n);
                    $s = exp($G($t));
                    $env = array_map(static fn($v) => $s * $v, $prof);
                    Solitons::addScaled($Q, $c, $env, array_fill(0, $n, $mu * $tau($t)));
                    return $Q;
                };
            } else {
                $this->checkClassicalReference($eta, $xi);
                $L = $P['L'];
                $exact = static fn(float $t): array => Solitons::manakovPeriodic($x, $L, $tau($t), $eta, $xi, $c, $x0, exp($G($t)));
            }
            return [$exact, $gammaEx, null];
        }

        // Caputo: manufactured solution q* = c (1 + t^nu + i t^2) psi(x)
        $nu = $P['nu'];
        $psR = $psI = [];
        foreach ($x as $i => $xv) {
            $a = 2.0 * $eta / cosh(2.0 * $eta * ($xv - $x0));
            $psR[$i] = $a * cos(-2.0 * $xi * $xv);
            $psI[$i] = $a * sin(-2.0 * $xi * $xv);
        }
        $xr = $psR; $xim = $psI;
        FFT::transform($xr, $xim);
        foreach ($this->k as $i => $kv) {
            $xr[$i] *= -$kv * $kv;
            $xim[$i] *= -$kv * $kv;
        }
        FFT::transform($xr, $xim, true);
        $g1 = Gamma::gamma(1.0 + $nu);
        $g3 = Gamma::gamma(3.0 - $nu);
        // complex helpers on scalars: [re, im]
        $A = static fn(float $t): array => [1.0 + $t ** $nu, $t * $t];
        $dA = static fn(float $t): array => [$g1, 2.0 * $t ** (2.0 - $nu) / $g3];
        $exact = static function (float $t) use ($A, $c, $psR, $psI): array {
            [$ar, $ai] = $A($t);
            $Q = [[], [], [], []];
            for ($j = 0; $j < 2; $j++) {
                $cr = $c[2 * $j] * $ar - $c[2 * $j + 1] * $ai;   // c_j * A
                $ci = $c[2 * $j] * $ai + $c[2 * $j + 1] * $ar;
                foreach ($psR as $i => $pr) {
                    $Q[2 * $j][$i] = $cr * $pr - $ci * $psI[$i];
                    $Q[2 * $j + 1][$i] = $cr * $psI[$i] + $ci * $pr;
                }
            }
            return $Q;
        };
        $gam = $this->gamma; $gf = $this->g; $V = $this->V;
        $forcing = static function (float $t) use ($exact, $A, $dA, $c, $psR, $psI, $xr, $xim, $beta, $gam, $gf, $V, $x): array {
            $q = $exact($t);
            [$ar, $ai] = $A($t);
            [$dr, $di] = $dA($t);
            $bt = $beta($t); $gm = 2.0 * $gam($t); $gt = $gf($t);
            $Vx = $V !== null ? $V($x, $t) : null;
            $out = [[], [], [], []];
            for ($j = 0; $j < 2; $j++) {
                $cr = $c[2 * $j]; $cim = $c[2 * $j + 1];
                // i c dA  and  beta c A
                $u = $cr * $dr - $cim * $di; $v = $cr * $di + $cim * $dr;      // c*dA
                $iu = -$v; $iv = $u;                                            // i*c*dA
                $w = $bt * ($cr * $ar - $cim * $ai); $z = $bt * ($cr * $ai + $cim * $ar);
                foreach ($psR as $i => $pr) {
                    $I = $q[0][$i] ** 2 + $q[1][$i] ** 2 + $q[2][$i] ** 2 + $q[3][$i] ** 2;
                    $s = $gm * $I + ($Vx !== null ? $Vx[$i] : 0.0);
                    $qr = $q[2 * $j][$i]; $qi = $q[2 * $j + 1][$i];
                    $out[2 * $j][$i] = $iu * $pr - $iv * $psI[$i] + $w * $xr[$i] - $z * $xim[$i] + $s * $qr + $gt * $qi;
                    $out[2 * $j + 1][$i] = $iu * $psI[$i] + $iv * $pr + $w * $xim[$i] + $z * $xr[$i] + $s * $qi - $gt * $qr;
                }
            }
            return $out;
        };
        $this->warnings[] = 'Caputo exact mode uses a manufactured solution q* = c(1 + t^nu + i t^2) psi(x) with a computed forcing term.';
        return [$exact, null, $forcing];
    }

    /** Expected order of the time discretization: 2 for Strang (also nu = 1), min(r nu, 2 - nu) for L1. */
    public static function expectedOrder(array $P): float
    {
        return $P['model'] === 'caputo' && $P['nu'] < 1.0 ? min($P['r'] * $P['nu'], 2.0 - $P['nu']) : 2.0;
    }

    private function makeSolver(?\Closure $gammaOverride, ?\Closure $forcing): SplitStepSolver|CaputoSolver
    {
        $P = $this->P;
        $gam = $gammaOverride ?? $this->gamma;
        if ($P['model'] === 'caputo' && $P['nu'] >= 1.0) {
            return new SplitStepSolver($P['L'], $P['Nx'], 2.0, $this->beta, $gam, $this->g, $this->V, $this->deadline);
        }
        if ($P['model'] === 'caputo') {
            $s = new CaputoSolver($P['L'], $P['Nx'], $P['nu'], $this->beta, $gam, $this->g, $this->V, $forcing, $this->deadline);
            // Anderson depth; 0 gives plain fixed-point (Picard) iteration. Not exposed in the web API.
            $s->anderson = (int)($P['anderson'] ?? 10);
            return $s;
        }
        $alpha = $P['model'] === 'space' ? $P['alpha'] : 2.0;
        return new SplitStepSolver($P['L'], $P['Nx'], $alpha, $this->beta, $gam, $this->g, $this->V, $this->deadline);
    }

    private function runSolver(SplitStepSolver|CaputoSolver $s, array $Q0, int $nSave, callable $onSave, ?float $dt = null, ?int $N = null): array
    {
        $P = $this->P;
        if ($s instanceof CaputoSolver) {
            return $s->solve($Q0, $P['T'], $N ?? $P['N'], $P['r'], $nSave, $onSave);
        }
        if ($P['model'] === 'caputo') {                 // nu = 1: the step is always T/N
            return $s->solve($Q0, $P['T'], $P['T'] / ($N ?? $P['N']), $nSave, $onSave);
        }
        return $s->solve($Q0, $P['T'], $dt ?? $P['dt'], $nSave, $onSave);
    }

    /** Number of steps actually taken by the split-step solver for a requested dt (it uses T/N). */
    private function stepsFor(float $dt): int
    {
        return max(1, (int)round($this->P['T'] / $dt));
    }

    /** Total reference floor (absolute); see $floor. */
    private function refFloor(): float
    {
        return array_sum($this->floor);
    }

    /**
     * Relative discrepancy between the quadrature that the split-step solver applies (Simpson on half steps for
     * beta and g; the gain-weighted Simpson step of the nonlinear phase for gamma) and the reference integrals,
     * for a run with N steps. Large values mean that dt does not resolve the
     * coefficient, so the run follows the discretized coefficient (e.g. the power follows exp(2 sum g)).
     */
    private function stepQuadrature(int $N, \Closure $gamEff): array
    {
        if (isset($this->quadCache[$N])) {
            return $this->quadCache[$N];
        }
        if ($this->Ggamma === null) {
            $G = $this->G;
            [$this->Ggamma] = Integral::antiderivativeAdaptive(static fn(float $t): float => $gamEff($t) * exp(2.0 * $G($t)), $this->P['T']);
        }
        $T = $this->P['T'];
        $q = [
            'beta' => Integral::stepQuadratureError($this->beta, $this->tau, $T, $N),
            'gamma' => Integral::nonlinearQuadratureError($gamEff, $this->g, $this->Ggamma, $this->G, $T, $N),
            'g' => Integral::stepQuadratureError($this->g, $this->G, $T, $N),
        ];
        if ($this->V !== null && $this->Vsample !== null) {     // a time-independent V is integrated exactly
            $q['V'] = Integral::panelQuadratureErrorSampled($this->Vsample, $T, $N);
        }
        return $this->quadCache[$N] = $q;
    }

    /**
     * Coefficients that the output times (every k steps) under-sample: between consecutive outputs, the value at the
     * fraction theta = 0.382 of the interval differs from the linear interpolation of the outputs by more than
     * OUTPUT_ALIAS of max(1, max|f|). An irrational fraction avoids sampling an oscillation at its own nodes.
     * @return list<string>
     */
    private function outputAliasing(int $N, \Closure $gamEff): array
    {
        $k = max(1, intdiv($N, self::N_SAVE));
        if ($k === 1) {
            return [];
        }
        $T = $this->P['T'];
        $D = $k * $T / $N;
        $th = 0.5 * (3.0 - sqrt(5.0));
        $out = [];
        foreach (['beta' => $this->beta, 'gamma' => $gamEff, 'g' => $this->g] as $name => $f) {
            $d = 0.0;
            $scale = 1.0;
            $f0 = $f(0.0);
            for ($j = 0; ($j + 1) * $D <= $T * (1.0 + 1e-12); $j++) {
                $f1 = $f(($j + 1) * $D);
                $d = max($d, abs($f($j * $D + $th * $D) - ((1.0 - $th) * $f0 + $th * $f1)));
                $scale = max($scale, abs($f0), abs($f1));
                $f0 = $f1;
            }
            if ($d > self::OUTPUT_ALIAS * $scale) {
                $out[] = $name;
            }
        }
        return $out;
    }

    /** Maximum exact error over every time step of a split-step run (outputs at every step). */
    private function errorAtEveryStep(\Closure $exact, ?\Closure $gammaEx, ?\Closure $forcing, array $Q0, ?float $dt, ?int $N): float
    {
        $s = $this->makeSolver($gammaEx, $forcing);
        $e = 0.0;
        $this->runSolver($s, $Q0, PHP_INT_MAX, static function (float $t, array $Q) use (&$e, $exact): void {
            $e = max($e, self::maxDiff($Q, $exact($t)));
        }, $dt, $N);
        return $e;
    }

    /** max over x and components of the complex modulus |q_j|. */
    private static function peakModulus(array $A): float
    {
        $m = 0.0;
        for ($c = 0; $c < 4; $c += 2) {
            foreach ($A[$c] as $i => $v) {
                $m = max($m, sqrt($v * $v + $A[$c + 1][$i] * $A[$c + 1][$i]));
            }
        }
        return $m;
    }

    /** max over x and components of the complex modulus |q_j - p_j|. */
    private static function maxDiff(array $A, array $B): float
    {
        $m = 0.0;
        for ($c = 0; $c < 4; $c += 2) {
            $ar = $A[$c]; $ai = $A[$c + 1]; $br = $B[$c]; $bi = $B[$c + 1];
            foreach ($ar as $i => $v) {
                $dr = $v - $br[$i];
                $di = $ai[$i] - $bi[$i];
                $d = sqrt($dr * $dr + $di * $di);
                if (is_nan($d)) {
                    return INF;
                }
                if ($d > $m) {
                    $m = $d;
                }
            }
        }
        return $m;
    }

    private static function moduli(array $Q): array
    {
        $a1 = $a2 = [];
        foreach ($Q[0] as $i => $r) {
            $a1[$i] = sqrt($r * $r + $Q[1][$i] ** 2);
            $a2[$i] = sqrt($Q[2][$i] ** 2 + $Q[3][$i] ** 2);
        }
        return [$a1, $a2];
    }

    private static function r(array $a, int $digits = 7): array
    {
        $out = [];
        foreach ($a as $v) {
            $out[] = is_finite($v) ? round($v, $digits) : null;
        }
        return $out;
    }

    private static function sci(float $v): float
    {
        return is_finite($v) ? (float)sprintf('%.12e', $v) : 0.0;
    }

    private static function indices(int $m, int $n): array
    {
        if ($m <= $n) {
            return range(0, $m - 1);
        }
        $idx = [];
        for ($i = 0; $i < $n; $i++) {
            $idx[(int)round($i * ($m - 1) / ($n - 1))] = true;
        }
        return array_keys($idx);
    }

    // -------------------------------------------------------------------- run
    public function run(): array
    {
        $start = microtime(true);
        $P = $this->P;
        $model = $P['model'];
        $exact = $gammaEx = $forcing = null;
        if ($P['exact']) {
            [$exact, $gammaEx, $forcing] = $this->setupExact();
            $Q0 = $exact(0.0);
            // rounding level of the comparison: errors near it carry no information on the order
            $this->floor['roundoff'] = self::ROUNDOFF * self::peakModulus($Q0);
        } else {
            $Q0 = $this->initialCondition();
        }
        $solver = $this->makeSolver($gammaEx, $forcing);
        $alpha = $model === 'space' ? $P['alpha'] : 2.0;
        $Dk = array_map(static fn($kv) => abs($kv) ** $alpha, $this->k);
        $Nx = $P['Nx'];
        $dx = $this->dx;
        $ix = self::indices($Nx, 256);
        $beta = $this->beta; $G = $this->G; $V = $this->V;
        $gamEff = $gammaEx ?? $this->gamma;

        $ts = $P1 = $P2 = $peak = $H = $gain = $err = [];
        $frames1 = $frames2 = [];
        $initial = null;
        $final = null;
        $onSave = function (float $t, array $Q) use (&$ts, &$P1, &$P2, &$peak, &$H, &$gain, &$err, &$frames1, &$frames2,
                                                    &$initial, &$final, $exact, $Dk, $dx, $Nx, $ix, $beta, $G, $V, $gamEff) {
            [$r1, $i1, $r2, $i2] = $Q;
            $s1 = $s2 = $pk = $I2sum = $VI = 0.0;
            $f1 = $f2 = [];
            $Vx = $V !== null ? $V($this->x, $t) : null;
            foreach ($r1 as $i => $a) {
                $a1 = $a * $a + $i1[$i] * $i1[$i];
                $a2 = $r2[$i] * $r2[$i] + $i2[$i] * $i2[$i];
                $s1 += $a1; $s2 += $a2;
                $I = $a1 + $a2;
                if ($I > $pk) { $pk = $I; }
                $I2sum += $I * $I;
                if ($Vx !== null) { $VI += $Vx[$i] * $I; }
            }
            foreach ($ix as $j) {
                $f1[] = $r1[$j] * $r1[$j] + $i1[$j] * $i1[$j];
                $f2[] = $r2[$j] * $r2[$j] + $i2[$j] * $i2[$j];
            }
            $a = $r1; $b = $i1; $c = $r2; $d = $i2;
            FFT::transform($a, $b);
            FFT::transform($c, $d);
            $kin = 0.0;
            foreach ($Dk as $i => $w) {
                $kin += $w * ($a[$i] * $a[$i] + $b[$i] * $b[$i] + $c[$i] * $c[$i] + $d[$i] * $d[$i]);
            }
            $ts[] = $t;
            $P1[] = $s1 * $dx; $P2[] = $s2 * $dx; $peak[] = $pk;
            $H[] = $beta($t) * $kin * $dx / $Nx - $gamEff($t) * $I2sum * $dx - $VI * $dx;
            $gain[] = exp(2.0 * $G($t));
            $frames1[] = $f1; $frames2[] = $f2;
            if ($exact !== null) {
                $err[] = self::maxDiff($Q, $exact($t));
            }
            if ($initial === null) {
                $initial = $Q;
            }
            $final = $Q;
        };

        $quad = null;
        $Nmain = 0;
        if ($solver instanceof SplitStepSolver) {
            $Nmain = $model === 'caputo' ? $P['N'] : $this->stepsFor($P['dt']);
            $quad = $this->stepQuadrature($Nmain, $gamEff);
            foreach ($quad as $name => $d) {
                if ($d > self::QUAD_TOL) {
                    $this->warnings[] = sprintf('dt = %.4g does not resolve %s: the quadrature that the solver applies (%s) '
                        . 'differs from the integral by %.1e (relative). The run then follows the discretized coefficient%s. '
                        . 'Reduce dt.', $P['T'] / $Nmain, $name === 'V' ? 'V(x,t)' : "$name(t)",
                        match ($name) {
                            'gamma' => 'one gain-weighted Simpson panel per step for the nonlinear phase',
                            'V' => sprintf('one Simpson panel per step, checked at %d grid points', $this->Vpoints),
                            default => 'Simpson on half steps',
                        },
                        $d, $name === 'g' ? '; in particular the power follows exp(2 x step sum of g), not exp(2 int g)' : '');
                }
            }
        }
        $this->runSolver($solver, $Q0, self::N_SAVE, $onSave);
        foreach ($final as $col) {
            foreach ($col as $v) {
                if (!is_finite($v)) {
                    throw new SimulationException('The solution blew up (non-finite values). Reduce dt / increase N, or check the coefficients.', $model === 'caputo' ? 'N' : 'dt');
                }
            }
        }
        $tMain = microtime(true) - $start;

        // ---- assemble result
        $it = self::indices(count($ts), 120);
        $sx = array_map(fn($j) => $this->x[$j], $ix);
        $sI1 = $sI2 = [];
        foreach ($it as $j) {
            $sI1[] = self::r($frames1[$j], 6);
            $sI2[] = self::r($frames2[$j], 6);
        }
        [$fa1, $fa2] = self::moduli($final);
        [$ia1, $ia2] = self::moduli($initial);
        $res = [
            'model' => $model,
            'x' => self::r($this->x),
            't' => self::r($ts),
            'surface' => ['x' => self::r($sx), 't' => self::r(array_map(fn($j) => $ts[$j], $it)), 'I1' => $sI1, 'I2' => $sI2],
            'final' => [
                'T' => end($ts),
                'q1' => ['re' => self::r($final[0]), 'im' => self::r($final[1]), 'abs' => self::r($fa1)],
                'q2' => ['re' => self::r($final[2]), 'im' => self::r($final[3]), 'abs' => self::r($fa2)],
            ],
            'initial' => ['abs1' => self::r($ia1), 'abs2' => self::r($ia2)],
            'diagnostics' => ['P1' => self::r($P1, 14), 'P2' => self::r($P2, 14), 'peak' => self::r($peak, 8),
                'H' => self::r($H, 14), 'gainFactor' => self::r($gain, 12)],
            'warnings' => $this->warnings,
        ];

        $sampling = null;
        if ($exact !== null) {
            $errMax = max($err);
            $sampling = 'outputs';
            $aliased = $solver instanceof SplitStepSolver ? $this->outputAliasing($Nmain, $gamEff) : [];
            if ($aliased) {
                // the output times can fall where the error happens to vanish (e.g. at the periods of beta)
                $sampling = 'every step';
                $errMax = max($errMax, $this->errorAtEveryStep($exact, $gammaEx, $forcing, $Q0, null, null));
                $this->warnings[] = sprintf('The output times (every %d steps) under-sample %s, so the error between them is not '
                    . 'represented by the stored outputs. The exact error (and the convergence study) is therefore measured at '
                    . 'every time step.', max(1, intdiv($Nmain, self::N_SAVE)), implode(', ', array_map(static fn($n) => "$n(t)", $aliased)));
            }
            [$e1, $e2] = self::moduli($exact(end($ts)));
            $res['exact'] = ['q1abs' => self::r($e1), 'q2abs' => self::r($e2),
                'err_t' => array_map([self::class, 'sci'], $err), 'err_final' => self::sci(end($err)), 'err_max' => self::sci($errMax),
                'errorSampling' => $sampling,
                'refFloor' => self::sci($this->refFloor()),
                'refFloorParts' => array_map([self::class, 'sci'], $this->floor)];
        }

        if ($P['convergence'] && $exact !== null) {
            $peak = max(self::peakModulus($Q0), self::peakModulus($exact((float)end($ts))));
            $res['convergence'] = $this->convergence($exact, $gammaEx, $forcing, $Q0, $errMax, $sampling === 'every step', $peak);
        }
        $res['warnings'] = $this->warnings;

        $steps = $model === 'caputo' ? $P['N'] : (int)round($P['T'] / $P['dt']);
        $res['meta'] = [
            'Nx' => $Nx, 'L' => $P['L'], 'dx' => $dx, 'T' => $P['T'], 'steps' => $steps,
            'dt' => $P['dt'] ?? null, 'N' => $P['N'] ?? null, 'r' => $P['r'] ?? null,
            'dtEffective' => $solver instanceof SplitStepSolver ? $P['T'] / ($model === 'caputo' ? $P['N'] : $this->stepsFor($P['dt'])) : null,
            'alpha' => $alpha, 'nu' => $P['nu'] ?? 1.0, 'expectedOrder' => self::expectedOrder($P),
            'profileResidual' => $this->profileResidual,
            'profile' => $this->profileInfo,
            'quadrature' => $quad,
            'exact' => (bool)$P['exact'], 'gammaDerived' => $gammaEx !== null,
            'picardMax' => $solver instanceof CaputoSolver ? $solver->maxPicard : null,
            'picardMean' => $solver instanceof CaputoSolver && $solver->stepsDone > 0 ? round($solver->sumPicard / $solver->stepsDone, 2) : null,
            'picardSteps' => $solver instanceof CaputoSolver ? $solver->picardSteps : null,
            'andersonDepth' => $solver instanceof CaputoSolver ? $solver->anderson : null,
            'runtime_main' => round($tMain, 3), 'runtime_total' => round(microtime(true) - $this->tCreated, 3),
            'version' => \MANAKOV_LAB_VERSION,
            'engine' => 'PHP ' . PHP_VERSION . (function_exists('opcache_get_status') && (opcache_get_status(false)['jit']['on'] ?? false) ? ' + JIT' : ''),
            'scheme' => [
                'classic' => 'Strang split-step Fourier (2nd order)',
                'space' => 'Strang split-step Fourier, Riesz symbol |k|^alpha (2nd order)',
                'caputo' => $solver instanceof SplitStepSolver
                    ? 'nu = 1: Strang split-step Fourier (2nd order)'
                    : 'L1 scheme on graded mesh + Fourier spectral, implicit (Anderson-accelerated Picard)',
            ][$model],
        ];
        return $res;
    }

    private function convergence(\Closure $exact, ?\Closure $gammaEx, ?\Closure $forcing, array $Q0, float $errMain,
                                 bool $dense, float $peak): array
    {
        $P = $this->P;
        $conv = ['h' => [], 'err' => [], 'order' => [], 'label' => '', 'expected' => null, 'refFloor' => $this->refFloor(),
            'unresolved' => [], 'saturated' => [], 'preasymptotic' => [], 'untested' => [], 'dtIndependent' => null];
        $floor = 10.0 * $this->refFloor();   // errors below this are not attributable to the integrator
        $gamEff = $gammaEx ?? $this->gamma;
        $errorOf = function (array $runArgs) use ($exact, $gammaEx, $forcing, $Q0, $dense): ?float {
            $s = $this->makeSolver($gammaEx, $forcing);
            $e = 0.0;
            $cb = static function (float $t, array $Q) use (&$e, $exact): void {
                $e = max($e, self::maxDiff($Q, $exact($t)));
            };
            try {
                // same output times as the main run, so that every level is sampled alike (every step if $dense)
                $this->runSolver($s, $Q0, $dense ? PHP_INT_MAX : self::N_SAVE, $cb, $runArgs['dt'] ?? null, $runArgs['N'] ?? null);
            } catch (SimulationException $ex) {
                if (str_contains($ex->getMessage(), 'time limit')) {
                    throw $ex;
                }
                return null;
            }
            return is_finite($e) ? $e : null;
        };

        if ($P['model'] === 'caputo') {
            $N = $P['N'];
            $Ns = array_values(array_unique([max(16, intdiv($N, 8)), max(16, intdiv($N, 4)), max(16, intdiv($N, 2)), $N]));
            sort($Ns);
            $conv['label'] = 'N (time steps)';
            $conv['expected'] = self::expectedOrder($P);
            foreach ($Ns as $n) {
                $conv['h'][] = $n;
                $conv['err'][] = $n === $N ? $errMain : $errorOf(['N' => $n]);
            }
            foreach ($conv['h'] as $i => $h) {
                $e0 = $i ? $conv['err'][$i - 1] : null;
                $e1 = $conv['err'][$i];
                $conv['order'][] = ($i && $e0 && $e1 && $e1 > $floor) ? log($e0 / $e1) / log($h / $conv['h'][$i - 1]) : null;
            }
        } else {
            $dt = $P['dt'];
            $dts = array_values(array_filter([8 * $dt, 4 * $dt, 2 * $dt, $dt], static fn($d) => $d <= $P['T'] / 10)) ?: [$dt];
            $conv['label'] = 'dt';
            $conv['expected'] = 2.0;
            $bad = [];
            foreach ($dts as $j => $d) {
                $n = $this->stepsFor($d);
                $conv['h'][] = $P['T'] / $n;                 // the step actually taken (T/N), not the requested one
                $conv['err'][] = $d == $dt ? $errMain : $errorOf(['dt' => $d]);
                $bad[$j] = max($this->stepQuadrature($n, $gamEff)) > self::QUAD_TOL;
                if ($bad[$j]) {
                    $conv['unresolved'][] = $P['T'] / $n;
                }
            }
            foreach ($conv['h'] as $i => $h) {
                $e0 = $i ? $conv['err'][$i - 1] : null;
                $e1 = $conv['err'][$i];
                $ok = $i && $e0 && $e1 && $e1 > $floor && !$bad[$i - 1] && !$bad[$i];
                $conv['order'][] = $ok ? log($e0 / $e1) / log($conv['h'][$i - 1] / $h) : null;
            }
        }
        $this->checkAsymptotic($conv, $P['model'] === 'caputo', $peak);
        return $conv;
    }

    /**
     * A posteriori test of the observed orders. With e = C h^p + S, S an error that does not decrease with the
     * time step (spatial resolution, reference, round-off), three consecutive levels give p from the differences
     * of the errors, which do not contain S, and then S (fitPowerLaw, valid for unequal step ratios). An order is
     * withheld when S exceeds SATURATION of the error of the finer level, when the errors do not decrease
     * monotonically or no p > 0 fits them (saturation), or when that error exceeds PREASYMPTOTIC times the peak
     * modulus of the solution (not yet in the asymptotic regime). With only two levels S cannot be separated, so the
     * order is withheld: as saturated if the error does not decrease, otherwise as 'untested'.
     */
    /**
     * Fit e = C h^p + S to three levels h_a > h_b > h_c (ratios need not be equal). The differences eliminate S:
     * (e_a - e_b)/(e_b - e_c) = (h_a^p - h_b^p)/(h_b^p - h_c^p) =: f(p), and f increases strictly with p > 0 from
     * log(h_a/h_b)/log(h_b/h_c), so p follows by bisection; then C = (e_b - e_c)/(h_b^p - h_c^p), S = e_c - C h_c^p.
     * Returns null when the errors do not decrease or when no p > 0 fits them.
     *
     * @return array{p: float, S: float}|null
     */
    public static function fitPowerLaw(float $ha, float $hb, float $hc, float $ea, float $eb, float $ec): ?array
    {
        $d1 = $ea - $eb;
        $d2 = $eb - $ec;
        if (!($d1 > 0.0 && $d2 > 0.0 && $ha > $hb && $hb > $hc && $hc > 0.0)) {
            return null;
        }
        $lA = log($ha / $hb);
        $lB = log($hb / $hc);
        $target = log($d1 / $d2);
        $lnf = static fn(float $p): float => $p * $lB + log(expm1($p * $lA)) - log(expm1($p * $lB));
        if ($target <= log($lA / $lB)) {                    // f(0+) = lA/lB: only p <= 0 fits
            return null;
        }
        $lo = 0.0;
        $hi = 30.0;
        if ($lnf($hi) <= $target) {
            $lo = $hi;
        } else {
            for ($k = 0; $k < 200 && $hi - $lo > 1e-13 * $hi; $k++) {
                $mid = 0.5 * ($lo + $hi);
                if ($lnf($mid) < $target) {
                    $lo = $mid;
                } else {
                    $hi = $mid;
                }
            }
        }
        $p = 0.5 * ($lo + $hi);
        $C = $d2 / ($hb ** $p - $hc ** $p);
        return ['p' => $p, 'S' => $ec - $C * $hc ** $p];
    }

    private function checkAsymptotic(array &$conv, bool $caputo, float $peak): void
    {
        $e = $conv['err'];
        $n = count($e);
        $hs = array_map(static fn($h) => $caputo ? 1.0 / $h : (float)$h, $conv['h']);
        $Smax = null;
        for ($i = 1; $i < $n; $i++) {
            if ($conv['order'][$i] === null) {
                continue;
            }
            if ($e[$i] > self::PREASYMPTOTIC * $peak) {
                $conv['order'][$i] = null;
                $conv['preasymptotic'][] = $conv['h'][$i];
                continue;
            }
            if ($n < 3) {
                if ($e[$i] >= $e[$i - 1]) {                 // the error does not decrease: nothing to attribute
                    $conv['order'][$i] = null;
                    $conv['saturated'][] = $conv['h'][$i];
                    $Smax = max($Smax ?? 0.0, $e[$i]);
                } else {                                   // decreasing, but S cannot be separated from C h^p
                    $conv['order'][$i] = null;
                    $conv['untested'][] = $conv['h'][$i];
                }
                continue;
            }
            [$a, $b, $c] = $i + 1 < $n ? [$i - 1, $i, $i + 1] : [$i - 2, $i - 1, $i];
            $sat = true;
            if ($e[$a] !== null && $e[$b] !== null && $e[$c] !== null) {
                $fit = self::fitPowerLaw($hs[$a], $hs[$b], $hs[$c], $e[$a], $e[$b], $e[$c]);
                if ($fit !== null) {
                    $sat = abs($fit['S']) > self::SATURATION * $e[$i];
                    if ($sat) {
                        $Smax = max($Smax ?? 0.0, abs($fit['S']));
                    }
                } else {
                    $Smax = max($Smax ?? 0.0, $e[$c]);
                }
            }
            if ($sat) {
                $conv['order'][$i] = null;
                $conv['saturated'][] = $conv['h'][$i];
            }
        }
        $conv['dtIndependent'] = $Smax;
        if ($conv['saturated']) {
            $this->warnings[] = sprintf('Convergence study: the errors do not decrease as a power of the time step at %d level(s); '
                . 'a component that does not depend on the step%s dominates there (spatial resolution, the reference or '
                . 'round-off), so those observed orders are omitted. Increase Nx (or L) and check again.',
                count($conv['saturated']), $Smax !== null ? sprintf(', about %.1e,', $Smax) : '');
        }
        if ($conv['untested']) {
            $this->warnings[] = 'Convergence study: the study has only two levels, so the a posteriori saturation test, '
                . 'which needs three, cannot be applied and the observed order is omitted. '
                . ($caputo ? 'Use N >= 34 for three levels or N >= 68 for four.'
                           : 'Use dt <= T/40 for three levels or dt <= T/80 for four (levels need dt <= T/10).');
        }
        if ($conv['preasymptotic']) {
            $this->warnings[] = sprintf('Convergence study: at %d level(s) the error exceeds %d%% of the peak amplitude, so the '
                . 'run is not in the asymptotic regime and those observed orders are omitted. Reduce the time step.',
                count($conv['preasymptotic']), (int)round(100 * self::PREASYMPTOTIC));
        }
    }
}
