<?php
declare(strict_types=1);

/**
 * Dependency-free test suite for the numerical core.   php tests/run_tests.php
 * Exit status 0 when every test passes.
 */

if (PHP_SAPI !== 'cli') {
    http_response_code(404);
    exit;
}

require dirname(__DIR__) . '/src/bootstrap.php';
require __DIR__ . '/lib_problem.php';

use ManakovLab\ExprException;
use ManakovLab\ExprParser;
use ManakovLab\Numerics\CaputoSolver;
use ManakovLab\Numerics\FFT;
use ManakovLab\Numerics\Integral;
use ManakovLab\Numerics\Simulation;
use ManakovLab\Numerics\SimulationException;
use ManakovLab\Numerics\SplitStepSolver;
use ManakovLab\Numerics\Solitons;

$failures = 0;
$count = 0;

function check(string $name, callable $fn): void
{
    global $failures, $count;
    $count++;
    $t = microtime(true);
    try {
        $fn();
        printf("ok    %-70s %6.2f s\n", $name, microtime(true) - $t);
    } catch (\Throwable $e) {
        $failures++;
        printf("FAIL  %-70s\n      %s\n", $name, $e->getMessage());
    }
}

function assertTrue(bool $c, string $msg): void
{
    if (!$c) {
        throw new RuntimeException($msg);
    }
}

function expectSimulationError(array $problem, string $field): void
{
    try {
        (new Simulation($problem))->run();
    } catch (SimulationException $e) {
        assertTrue($e->field === $field, "expected field '$field', got '" . ($e->field ?? 'null') . "': " . $e->getMessage());
        return;
    }
    throw new RuntimeException('expected a SimulationException');
}

function hasWarning(array $res, string $needle): bool
{
    foreach ($res['warnings'] as $w) {
        if (str_contains($w, $needle)) {
            return true;
        }
    }
    return false;
}

check('FFT: forward + inverse is the identity', function (): void {
    $n = 64; $re = $im = [];
    for ($i = 0; $i < $n; $i++) { $re[] = sin(0.3 * $i) + 0.1 * $i; $im[] = cos(0.7 * $i); }
    $a = $re; $b = $im;
    FFT::transform($a, $b);
    FFT::transform($a, $b, true);
    $e = 0.0;
    foreach ($a as $i => $v) { $e = max($e, abs($v - $re[$i]), abs($b[$i] - $im[$i])); }
    assertTrue($e < 1e-12, "round-off error $e");
});

check('Expression parser: whitelist, empty input, errors with position', function (): void {
    $p = new ExprParser(['t']);
    assertTrue(ExprParser::toCanonical($p->parse('')) === '0', 'empty expression should be 0');
    foreach (["system('ls')", 'x + 1', 'foo(t)', '1 +', str_repeat('1+', 200) . '1'] as $bad) {
        try {
            $p->parse($bad);
            throw new RuntimeException("accepted '$bad'");
        } catch (ExprException) {
        }
    }
});

check('Solitons: convention exp(-2 i xi x), velocity -4 xi', function (): void {
    $r = (new Simulation(problemFrom('classic-exact', ['coeffs' => ['beta' => '1', 'g' => '0'],
        'solitons' => [['eta' => 0.5, 'xi' => 0.25, 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0, 'p2im' => 0]],
        'params' => ['T' => 4, 'dt' => 0.01], 'exact' => false, 'convergence' => false])))->run();
    $a = $r['final']['q1']['abs'];
    $imax = array_keys($a, max($a))[0];
    $xPeak = $r['x'][$imax];
    assertTrue(abs($xPeak + 4.0) < 0.2, "peak at x = $xPeak, expected -4 (= -4 xi T)");
});

check('Petviashvili: residual of the fractional soliton below 1e-12', function (): void {
    $L = 80.0; $Nx = 512; $x = [];
    for ($i = 0; $i < $Nx; $i++) { $x[] = ($i - $Nx / 2) * $L / $Nx; }
    [, $res] = Solitons::fractional($x, FFT::wavenumbers($Nx, $L), 1.5, 1.0);
    assertTrue($res < 1e-12, "residual $res");
});

check('Split-step: weighted power conserved to round-off (general coefficients)', function (): void {
    $r = (new Simulation(problemFrom('classic-trap', ['coeffs' => ['g' => '0.05 sin(t)'], 'params' => ['T' => 5]])))->run();
    $d = $r['diagnostics']; $m = 0.0;
    foreach ($d['P1'] as $i => $p1) {
        $m = max($m, abs(($p1 + $d['P2'][$i]) / $d['gainFactor'][$i] / ($d['P1'][0] + $d['P2'][0]) - 1.0));
    }
    assertTrue($m < 1e-11, "drift $m");
});

check('Classical exact solution: observed orders 2 (preset classic-exact)', function (): void {
    $c = (new Simulation(problemFrom('classic-exact')))->run()['convergence'];
    foreach (array_slice($c['order'], 1) as $o) { assertTrue(abs($o - 2.0) < 0.03, 'order ' . json_encode($c['order'])); }
});

check('Caputo manufactured solution: orders 2 - nu (graded) and nu (uniform)', function (): void {
    $c = (new Simulation(problemFrom('caputo-mms')))->run()['convergence'];
    foreach (array_slice($c['order'], 1) as $o) { assertTrue(abs($o - 1.4) < 0.03, 'graded ' . json_encode($c['order'])); }
    $c = (new Simulation(problemFrom('caputo-mms', ['params' => ['r' => 1]])))->run()['convergence'];
    foreach (array_slice($c['order'], 1) as $o) { assertTrue(abs($o - 0.6) < 0.03, 'uniform ' . json_encode($c['order'])); }
});

check('Caputo: nonlinear solve converges for small nu (0.1 and 0.05)', function (): void {
    foreach ([0.1, 0.05] as $nu) {
        $r = (new Simulation(problemFrom('caputo-decay', ['params' => ['nu' => $nu, 'T' => 1, 'N' => 64, 'Nx' => 128]])))->run();
        assertTrue($r['meta']['picardMax'] <= 40, "nu = $nu: max iterations {$r['meta']['picardMax']}");
    }
});

check('Caputo: manufactured solution with amplitude 1 converges at nu = 0.6', function (): void {
    $sol = [['eta' => 0.5, 'xi' => -0.25, 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0.5, 'p2im' => -0.5]];
    $r = (new Simulation(problemFrom('caputo-mms', ['solitons' => $sol, 'params' => ['N' => 64], 'convergence' => false])))->run();
    assertTrue((float)$r['exact']['err_max'] < 5e-2, 'error ' . $r['exact']['err_max']);
});

check('Caputo with nu = 1 runs the conservative split-step solver, with a warning', function (): void {
    $r = (new Simulation(problemFrom('caputo-decay', ['params' => ['nu' => 1.0, 'T' => 3, 'N' => 300]])))->run();
    $d = $r['diagnostics'];
    $P = (end($d['P1']) + end($d['P2'])) / ($d['P1'][0] + $d['P2'][0]);
    assertTrue(abs($P - 1.0) < 1e-10, "P(T)/P(0) = $P");
    assertTrue(hasWarning($r, 'nu = 1'), 'missing warning');
    expectSimulationError(problemFrom('caputo-mms', ['params' => ['nu' => 1.0]]), 'nu');
});

check('Caputo close to nu = 1 carries a dissipation warning', function (): void {
    $r = (new Simulation(problemFrom('caputo-decay', ['params' => ['nu' => 0.97, 'T' => 1, 'N' => 50, 'Nx' => 128]])))->run();
    assertTrue(hasWarning($r, 'nu close to 1'), 'missing warning');
});

check('alpha <= 1: warning on every run, exact mode refused', function (): void {
    $r = (new Simulation(problemFrom('space-collision', ['params' => ['alpha' => 0.8, 'T' => 1, 'Nx' => 256, 'L' => 80]])))->run();
    assertTrue(hasWarning($r, 'alpha <= 1'), 'missing warning');
    expectSimulationError(problemFrom('space-exact', ['params' => ['alpha' => 1.0]]), 'alpha');
});

check('Non-finite coefficient reported with its field', function (): void {
    expectSimulationError(problemFrom('classic-trap', ['coeffs' => ['beta' => '1/(t-t)']]), 'beta');
});

// ---------------------------------------------------------------- validity of the references and domain
function boundarySoliton(float $L, int $Nx, float $T): array
{
    return problemFrom('classic-exact', ['params' => ['L' => $L, 'Nx' => $Nx, 'T' => $T, 'dt' => 0.01],
        'coeffs' => ['beta' => '1', 'gamma' => '1', 'g' => '0', 'V' => '0'],
        'solitons' => [['eta' => 0.5, 'xi' => -0.5, 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0, 'p2im' => 0]],
        'convergence' => true]);
}

check('Exact mode: a soliton crossing x = +-L/2 is compared with the periodized reference, order 2', function (): void {
    foreach ([7.0, 14.0] as $T) {                                       // centre at x = 14 and x = 28 > L/2 = 20
        $r = (new Simulation(boundarySoliton(40.0, 512, $T)))->run();
        $o = $r['convergence']['order'];
        assertTrue(abs(end($o) - 2.0) < 0.02, "T = $T: order " . end($o));
        assertTrue($r['exact']['refFloorParts']['images'] < 1e-14, 'images floor ' . $r['exact']['refFloorParts']['images']);
    }
});

check('Exact mode refused when the soliton overlaps its periodic images', function (): void {
    expectSimulationError(boundarySoliton(12.0, 128, 1.0), 'L');       // 2 sech(eta L) = 9.9e-3 > 1e-3
});

check('Exact mode: spatial floor reported and orders below it omitted (Nx = 128); order 2 at Nx = 512', function (): void {
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.0025]])))->run();
    assertTrue(hasWarning($r, 'does not resolve the soliton spectrum'), 'missing warning');
    $C = $r['convergence'];
    assertTrue($C['refFloor'] > 1e-4, 'floor ' . $C['refFloor']);
    foreach ($C['order'] as $i => $o) {
        assertTrue($o === null || $C['err'][$i] > 10 * $C['refFloor'], "order $i not omitted");
    }
    assertTrue(end($C['order']) === null, 'finest order not omitted');
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['T' => 1, 'Nx' => 512, 'dt' => 0.0025]])))->run();
    assertTrue(abs(end($r['convergence']['order']) - 2.0) < 0.01, 'order ' . end($r['convergence']['order']));
});

check('Reference quadrature: adaptive table resolves beta = 300 sin(40 pi t) on [0, 100]', function (): void {
    $f = static fn(float $t): float => 300.0 * sin(40.0 * M_PI * $t);
    $F = static fn(float $t): float => 300.0 * (1.0 - cos(40.0 * M_PI * $t)) / (40.0 * M_PI);
    $fixed = Integral::antiderivative($f, 100.0);
    [$adapt, $info] = Integral::antiderivativeAdaptive($f, 100.0);
    $e0 = $e1 = 0.0;
    for ($i = 0; $i <= 4000; $i++) {
        $t = 99.9 * $i / 4000 + 0.0123;      // inside [0, T]
        $e0 = max($e0, abs($fixed($t) - $F($t)));
        $e1 = max($e1, abs($adapt($t) - $F($t)));
    }
    assertTrue($info['resolved'] && $info['m'] > 4096, 'm = ' . $info['m']);
    assertTrue($e0 > 1e-2 && $e1 < 1e-7, "fixed $e0, adaptive $e1");
    [, $info] = Integral::antiderivativeAdaptive(static fn(float $t): float => 1.0 + 0.5 * sin($t), 5.0);
    assertTrue($info['resolved'] && $info['m'] === 4096, 'smooth coefficient: m = ' . $info['m']);
});

check('Fractional exact mode: homogeneous and non-converged profiles reported, orders omitted', function (): void {
    $sol = static fn(float $eta, float $x0): array => [['eta' => $eta, 'xi' => 0, 'x0' => $x0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0.6, 'p2im' => 0.4]];
    $r = (new Simulation(problemFrom('space-exact', ['params' => ['L' => 5, 'Nx' => 128, 'T' => 0.5, 'dt' => 0.0025],
        'solitons' => $sol(0.25, 0.0)])))->run();
    assertTrue($r['meta']['profile']['branch'] === 'homogeneous' && hasWarning($r, 'spatially constant'), 'branch');
    assertTrue(count(array_filter($r['convergence']['order'], static fn($o) => $o !== null)) === 0, 'orders at round-off');
    $r = (new Simulation(problemFrom('space-exact', ['params' => ['L' => 10, 'Nx' => 128, 'T' => 1, 'dt' => 0.005],
        'solitons' => $sol(0.2495, 0.0390625)])))->run();
    assertTrue($r['meta']['profile']['iterations'] === null && hasWarning($r, 'did not converge'), 'convergence');
    assertTrue($r['meta']['profile']['branch'] === 'periodic', 'branch ' . $r['meta']['profile']['branch']);
    assertTrue($r['exact']['refFloorParts']['profile'] > 1e-7, 'profile floor');
    assertTrue(count(array_filter($r['convergence']['order'], static fn($o) => $o !== null)) === 0, 'orders not omitted');
});

check('Split-step: under-resolved gain g(t) reported; convergence uses the steps actually taken', function (): void {
    $c = ['beta' => '0', 'gamma' => '0', 'g' => 'cos(8 pi t)', 'V' => '0'];
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.5], 'coeffs' => $c])))->run();
    assertTrue(hasWarning($r, 'does not resolve g(t)') && $r['meta']['quadrature']['g'] > 1e-2, 'missing warning');
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.25], 'coeffs' => $c])))->run();
    assertTrue(!hasWarning($r, 'does not resolve'), 'spurious warning');
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['T' => 1, 'dt' => 0.0119]])))->run();
    foreach ($r['convergence']['h'] as $h) {
        $n = 1.0 / $h;
        assertTrue(abs($n - round($n)) < 1e-9, "h = $h is not T/N");
    }
    assertTrue(abs($r['meta']['dtEffective'] - 1.0 / 84) < 1e-15, 'dtEffective ' . $r['meta']['dtEffective']);
});

check('Solvers validate their own contract (direct library use)', function (): void {
    $one = static fn(float $t): float => 1.0;
    $zero = static fn(float $t): float => 0.0;
    $cases = [
        ['nu', static fn() => new CaputoSolver(2 * M_PI, 64, 1.1, $one, $zero, $zero, null)],
        ['nu', static fn() => new CaputoSolver(2 * M_PI, 64, 0.0, $one, $zero, $zero, null)],
        ['Nx', static fn() => new CaputoSolver(2 * M_PI, 60, 0.5, $one, $zero, $zero, null)],
        ['alpha', static fn() => new SplitStepSolver(2 * M_PI, 64, 2.5, $one, $zero, $zero, null)],
        ['L', static fn() => new SplitStepSolver(-1.0, 64, 2.0, $one, $zero, $zero, null)],
    ];
    foreach ($cases as [$field, $make]) {
        try {
            $make();
            assertTrue(false, "no exception for $field");
        } catch (SimulationException $e) {
            assertTrue($e->field === $field, "field {$e->field} instead of $field");
        }
    }
});

check('Fractional exact mode: off-node centre moved to a node, order 2 kept', function (): void {
    $sol = [['eta' => 0.5, 'xi' => 0, 'x0' => 0.15625, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0.6, 'p2im' => 0.4]];
    $r = (new Simulation(problemFrom('space-exact', ['params' => ['Nx' => 256, 'T' => 1], 'solitons' => $sol])))->run();
    assertTrue(hasWarning($r, 'nearest node'), 'missing warning');
    assertTrue($r['meta']['profileResidual'] < 1e-11, 'residual ' . $r['meta']['profileResidual']);
    $o = $r['convergence']['order'];
    assertTrue(abs(end($o) - 2.0) < 0.02, 'order ' . end($o));
});

check('Core rejects inputs outside the mathematical domain', function (): void {
    expectSimulationError(problemFrom('caputo-decay', ['params' => ['nu' => 1.1, 'N' => 50, 'T' => 1, 'Nx' => 128]]), 'nu');
    expectSimulationError(problemFrom('caputo-decay', ['params' => ['nu' => 0.0, 'N' => 50, 'T' => 1, 'Nx' => 128]]), 'nu');
    expectSimulationError(problemFrom('space-collision', ['params' => ['alpha' => 2.5, 'T' => 1]]), 'alpha');
    expectSimulationError(problemFrom('classic-collision', ['params' => ['dt' => -0.01, 'T' => 1]]), 'dt');
    expectSimulationError(problemFrom('classic-collision', ['params' => ['Nx' => 500, 'T' => 1]]), 'Nx');
});

check('Caputo linear mode against the Mittag-Leffler solution, nu = 1/2, order 3/2', function (): void {
    // i D^nu q + q_xx = 0 with q = c e^{ix}: q(t) = E_{1/2}(-i sqrt t) q(0) = e^{-t} (1 - i erfi(sqrt t)) q(0)
    $erfi = static function (float $x): float {
        $s = 0.0; $term = $x;                                  // x^{2n+1}/n!
        for ($n = 0; $n < 80; $n++) {
            $s += $term / (2 * $n + 1);
            $term *= $x * $x / ($n + 1);
        }
        return 2.0 / sqrt(M_PI) * $s;
    };
    $T = 0.25; $Nx = 64; $L = 2.0 * M_PI;
    $er = exp(-$T); $ei = -exp(-$T) * $erfi(sqrt($T));
    $one = static fn(float $t): float => 1.0;
    $zero = static fn(float $t): float => 0.0;
    $errs = [];
    foreach ([32, 64, 128, 256] as $N) {
        $s = new CaputoSolver($L, $Nx, 0.5, $one, $zero, $zero, null);
        $Q = [[], [], [], []];
        foreach ($s->x as $i => $xv) {
            $Q[0][$i] = cos($xv); $Q[1][$i] = sin($xv);
            $Q[2][$i] = 0.5 * cos($xv); $Q[3][$i] = 0.5 * sin($xv);
        }
        $Q0 = $Q;
        $last = null;
        $s->solve($Q, $T, $N, 3.0, 400, function (float $t, array $Qt) use (&$last): void { $last = $Qt; });
        $e = 0.0;
        for ($j = 0; $j < 4; $j += 2) {
            foreach ($Q0[$j] as $i => $a) {
                $b = $Q0[$j + 1][$i];
                $e = max($e, hypot($last[$j][$i] - ($er * $a - $ei * $b), $last[$j + 1][$i] - ($er * $b + $ei * $a)));
            }
        }
        $errs[] = $e;
    }
    $order = log($errs[2] / $errs[3]) / log(2.0);
    assertTrue(abs($errs[3] - 8.058e-5) < 1e-7, 'error at N = 256: ' . $errs[3]);
    assertTrue(abs($order - 1.5) < 0.05, "order $order");
});

check('Convergence study withholds orders whose errors saturate (Nx = 128, D2-H2)', function (): void {
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['Nx' => 128]])))->run();
    $c = $r['convergence'];
    assertTrue(count($c['saturated']) >= 2, 'saturated ' . json_encode($c['saturated']));
    foreach ($c['order'] as $o) {
        assertTrue($o === null, 'published order ' . json_encode($c['order']));
    }
    assertTrue(hasWarning($r, 'do not decrease as a power'), 'missing saturation warning');
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['T' => 4, 'Nx' => 128, 'dt' => 0.0025]])))->run();
    assertTrue(count(array_filter($r['convergence']['order'], fn($o) => $o !== null)) === 0, 'T = 4 orders');
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['Nx' => 256]])))->run();
    $c = $r['convergence'];
    assertTrue($c['saturated'] === [] && $c['preasymptotic'] === [], 'Nx = 256 flagged');
    assertTrue(abs(end($c['order']) - 2.0) < 0.01, 'order ' . end($c['order']));
});

check('Error is measured at every step when the outputs under-sample the coefficients (N1)', function (): void {
    $osc = ['beta' => '300 sin(40 pi t)', 'gamma' => '300 sin(40 pi t)', 'g' => '0', 'V' => '0'];
    $sol = [['eta' => 0.5, 'xi' => 2.5, 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0, 'p2im' => 0]];
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['L' => 80, 'Nx' => 256, 'T' => 10, 'dt' => 0.0025],
        'coeffs' => $osc, 'solitons' => $sol, 'convergence' => false])))->run();
    assertTrue($r['exact']['errorSampling'] === 'every step', 'sampling ' . $r['exact']['errorSampling']);
    assertTrue((float)$r['exact']['err_max'] > 0.1, 'err_max ' . $r['exact']['err_max']);
    assertTrue(hasWarning($r, 'under-sample'), 'missing aliasing warning');
    $r = (new Simulation(problemFrom('classic-exact', [])))->run();
    assertTrue($r['exact']['errorSampling'] === 'outputs', 'preset sampling');
});

check('Nonlinear-phase diagnostic reproduces the gain-weighted Simpson panel of the solver', function (): void {
    $co = ['beta' => '0', 'gamma' => 'cos(4 pi t)', 'g' => '0', 'V' => '0'];
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.5], 'coeffs' => $co])))->run();
    assertTrue(abs($r['meta']['quadrature']['gamma'] - 1.0 / 3.0) < 1e-6, 'gamma ' . $r['meta']['quadrature']['gamma']);
    assertTrue(hasWarning($r, 'gain-weighted Simpson'), 'missing warning');
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.25], 'coeffs' => $co])))->run();
    assertTrue($r['meta']['quadrature']['gamma'] < 1e-12, 'gamma at h = 0.25: ' . $r['meta']['quadrature']['gamma']);
    assertTrue(!hasWarning($r, 'does not resolve'), 'spurious warning');
});

check('Reference tables are cross-checked with Gauss-Legendre (aliased coefficient)', function (): void {
    $f = static fn(float $t): float => cos(32768 * M_PI * $t);
    [$F, $info] = Integral::antiderivativeAdaptive($f, 1.0);
    assertTrue($info['m'] > 4096 && $info['resolved'], 'm = ' . $info['m']);
    assertTrue(abs($F(1.0)) < 1e-6 && abs($F(0.3) - sin(32768 * M_PI * 0.3) / (32768 * M_PI)) < 1e-6, 'F(1) = ' . $F(1.0));
    $co = ['beta' => '0', 'gamma' => '0', 'g' => 'cos(32768 pi t)', 'V' => '0'];
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.5], 'coeffs' => $co])))->run();
    assertTrue($r['meta']['quadrature']['g'] > 0.5, 'g ' . $r['meta']['quadrature']['g']);
    assertTrue(hasWarning($r, 'does not resolve g(t)'), 'missing warning');
});

check('Saturation test fits unequal step ratios and needs three levels', function (): void {
    $h = [1 / 10, 1 / 20, 1 / 41, 1 / 81];
    $fit = Simulation::fitPowerLaw($h[1], $h[2], $h[3], $h[1] ** 2, $h[2] ** 2, $h[3] ** 2);
    assertTrue(abs($fit['p'] - 2.0) < 1e-9 && abs($fit['S']) < 1e-15, 'fit ' . json_encode($fit));
    $r = (new Simulation(problemFrom('caputo-mms', ['params' => ['N' => 36]])))->run();
    $c = $r['convergence'];                                   // levels N = 16, 18, 36
    assertTrue($c['saturated'] === [] && abs($c['order'][2] - 1.3336) < 1e-3, 'caputo ' . json_encode($c['order']));
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['T' => 0.25, 'dt' => 0.25 / 81]])))->run();
    $c = $r['convergence'];                                   // N = 10, 20, 41, 81
    assertTrue($c['saturated'] === [] && abs(end($c['order']) - 2.0) < 0.01, 'classic ' . json_encode($c['order']));
    $c0 = ['beta' => '1', 'gamma' => '1', 'g' => '0', 'V' => '0'];
    $sol = static fn($eta, $xi) => [['eta' => $eta, 'xi' => $xi, 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0.6, 'p2im' => 0.4]];
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['L' => 80, 'Nx' => 64, 'T' => 4, 'dt' => 0.16],
        'coeffs' => $c0, 'solitons' => $sol(0.35, 0)])))->run();
    $c = $r['convergence'];                                   // two levels, the error grows
    assertTrue(count($c['h']) === 2 && $c['order'] === [null, null] && $c['saturated'] === [0.16], 'growing ' . json_encode($c));
    $r = (new Simulation(problemFrom('classic-exact', ['params' => ['L' => 80, 'Nx' => 128, 'T' => 4, 'dt' => 0.16],
        'coeffs' => $c0, 'solitons' => $sol(0.5, -1)])))->run();
    $c = $r['convergence'];                                   // two levels, the error decreases
    assertTrue($c['order'] === [null, null] && $c['untested'] === [0.16], 'decreasing ' . json_encode($c['order']));
    assertTrue(hasWarning($r, 'only two levels'), 'missing two-level warning');
});

check('Potential quadrature: a time-dependent V(x,t) that dt does not resolve is reported', function (): void {
    $co = static fn(string $V) => ['beta' => '0', 'gamma' => '0', 'g' => '0', 'V' => $V];
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.5], 'coeffs' => $co('cos(8 pi t)')])))->run();
    assertTrue(abs($r['meta']['quadrature']['V'] - 1.0) < 1e-9 && hasWarning($r, 'does not resolve V(x,t)'), 'coarse');
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 1 / 64], 'coeffs' => $co('cos(8 pi t)')])))->run();
    assertTrue($r['meta']['quadrature']['V'] < 1e-6 && !hasWarning($r, 'does not resolve'), 'fine ' . $r['meta']['quadrature']['V']);
    $r = (new Simulation(problemFrom('classic-collision', ['params' => ['T' => 1, 'Nx' => 128, 'dt' => 0.5], 'coeffs' => $co('0.01 x^2')])))->run();
    assertTrue(!isset($r['meta']['quadrature']['V']), 'static V is integrated exactly');
});

check('Metadata state the scheme, grading and expected order actually used', function (): void {
    $r = (new Simulation(problemFrom('caputo-decay', ['params' => ['nu' => 0.05, 'r' => 10.0, 'N' => 32, 'T' => 0.5, 'Nx' => 128]])))->run();
    assertTrue(abs($r['meta']['expectedOrder'] - 0.5) < 1e-12 && $r['meta']['r'] == 10.0, 'nu = 0.05');
    $r = (new Simulation(problemFrom('caputo-decay', ['params' => ['nu' => 1.0, 'N' => 20, 'T' => 1, 'Nx' => 128]])))->run();
    assertTrue($r['meta']['expectedOrder'] == 2.0 && str_contains($r['meta']['scheme'], 'Strang') && abs($r['meta']['dtEffective'] - 0.05) < 1e-15, 'nu = 1');
    assertTrue($r['meta']['runtime_total'] >= $r['meta']['runtime_main'], 'timer');
});

printf("\n%d tests, %d failures\n", $count, $failures);
exit($failures ? 1 : 0);
