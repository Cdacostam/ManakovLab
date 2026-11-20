<?php
declare(strict_types=1);

/**
 * Regenerates every number reported in the article describing Manakov Lab (see the README).
 *
 *   php tests/reproduce.php [section ...] [--json=FILE]
 *
 * Sections: table3 mittag general fractional alpha table4 spatial nu1 solver collision timing   (default: all but timing)
 * With the OPcache JIT (recommended):
 *   php -d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=128M tests/reproduce.php
 *
 * This script is a command-line tool: it refuses to run inside a web server.
 */

if (PHP_SAPI !== 'cli') {
    http_response_code(404);
    exit;
}

require dirname(__DIR__) . '/src/bootstrap.php';
require __DIR__ . '/lib_problem.php';

use ManakovLab\ExprParser;
use ManakovLab\Presets;
use ManakovLab\Numerics\CaputoSolver;
use ManakovLab\Numerics\ExprCompiler;
use ManakovLab\Numerics\FFT;
use ManakovLab\Numerics\Simulation;
use ManakovLab\Numerics\Solitons;
use ManakovLab\Numerics\SplitStepSolver;

ini_set('memory_limit', '2G');

$args = array_slice($argv, 1);
$jsonFile = null;
$sections = [];
foreach ($args as $a) {
    if (str_starts_with($a, '--json=')) {
        $jsonFile = substr($a, 7);
    } else {
        $sections[] = $a;
    }
}
$all = ['table3', 'mittag', 'general', 'fractional', 'alpha', 'table4', 'spatial', 'nu1', 'solver', 'collision'];
if (!$sections || $sections === ['all']) {
    $sections = $all;
}
$out = ['version' => MANAKOV_LAB_VERSION, 'engine' => engine(), 'date' => gmdate('c')];

function engine(): string
{
    $jit = function_exists('opcache_get_status') && (opcache_get_status(false)['jit']['on'] ?? false);
    return 'PHP ' . PHP_VERSION . ($jit ? ' + JIT' : '');
}

function simulate(string $id, array $over = []): array
{
    return (new Simulation(problemFrom($id, $over)))->run();
}

function line(string $s = ''): void
{
    echo $s, "\n";
}

function convRow(string $label, array $conv): array
{
    $e = array_map(static fn($v) => sprintf('%.4e', $v), $conv['err']);
    $o = array_map(static fn($v) => $v === null ? '--' : sprintf('%.3f', $v), $conv['order']);
    line(sprintf('  %-46s h = %s', $label, implode(', ', $conv['h'])));
    line(sprintf('  %-46s err = %s', '', implode(', ', $e)));
    line(sprintf('  %-46s order = %s', '', implode(', ', $o)));
    return ['h' => $conv['h'], 'err' => $conv['err'], 'order' => $conv['order']];
}

function powerPeak(array $res): array
{
    $d = $res['diagnostics'];
    $P0 = $d['P1'][0] + $d['P2'][0];
    $PT = end($d['P1']) + end($d['P2']);
    return [$PT / $P0, end($d['peak'])];
}

/** Two-component field from a list of solitons, with the convention of Solitons::manakov. */
function solitonField(array $x, array $list): array
{
    $Q = Solitons::zeros(count($x));
    foreach ($list as $s) {
        $c = Solitons::normalise([$s['p1re'], $s['p1im'], $s['p2re'], $s['p2im']]);
        Solitons::addManakov($Q, $x, 0.0, $s['eta'], $s['xi'], $c, $s['x0']);
    }
    return $Q;
}

// ---------------------------------------------------------------- Table 3
if (in_array('table3', $sections, true)) {
    line('== Table 3: maximum error over x and output times, observed orders');
    $r = simulate('classic-exact');
    $out['table3']['classic_exact'] = convRow('Classical, exact (preset classic-exact)', $r['convergence']);
    $r = simulate('space-exact');
    $out['table3']['space_exact'] = convRow('Space-fractional alpha = 1.5, exact (space-exact)', $r['convergence']);
    $r = simulate('caputo-mms');
    $out['table3']['caputo_graded'] = convRow('Caputo nu = 0.6, manufactured, graded r = 7/3', $r['convergence']);
    $out['table3']['caputo_graded']['picard'] = [$r['meta']['picardMean'], $r['meta']['picardMax']];
    $r = simulate('caputo-mms', ['params' => ['r' => 1]]);
    $out['table3']['caputo_uniform'] = convRow('Caputo nu = 0.6, manufactured, uniform r = 1', $r['convergence']);
    $out['table3']['caputo_uniform']['picard'] = [$r['meta']['picardMean'], $r['meta']['picardMax']];
    line(sprintf('  nonlinear solve at N = 512 (mean / max iterations per step): graded %.2f / %d, uniform %.2f / %d',
        $out['table3']['caputo_graded']['picard'][0], $out['table3']['caputo_graded']['picard'][1], $r['meta']['picardMean'], $r['meta']['picardMax']));
    $g = $out['table3']['caputo_graded']['err'];
    $u = $out['table3']['caputo_uniform']['err'];
    $out['table3']['grading_gain_N512'] = end($u) / end($g);
    line(sprintf('  grading reduces the error at N = 512 by a factor %.1f', end($u) / end($g)));
    line();
}

// --------------------------- Table 3: Caputo linear mode against an analytic (Mittag-Leffler) solution
// i D^nu q + q_xx = 0, q = c e^{ix}, nu = 1/2: q(t) = E_{1/2}(-i sqrt t) q(0) = e^{-t}(1 - i erfi(sqrt t)) q(0).
// The reference is derived independently of the scheme (unlike the manufactured solution).
if (in_array('mittag', $sections, true)) {
    line('== Table 3: Caputo nu = 1/2, linear mode e^{ix}, exact E_{1/2}(-i sqrt t); T = 0.25, L = 2 pi, Nx = 64, r = 3');
    $erfi = static function (float $x): float {
        $s = 0.0; $term = $x;
        for ($n = 0; $n < 80; $n++) {
            $s += $term / (2 * $n + 1);
            $term *= $x * $x / ($n + 1);
        }
        return 2.0 / sqrt(M_PI) * $s;
    };
    $T = 0.25;
    $one = static fn(float $t): float => 1.0;
    $zero = static fn(float $t): float => 0.0;
    $Ns = [32, 64, 128, 256];
    $errs = [];
    foreach ($Ns as $N) {
        $s = new CaputoSolver(2.0 * M_PI, 64, 0.5, $one, $zero, $zero, null);
        $Q = [[], [], [], []];
        foreach ($s->x as $i => $xv) {
            $Q[0][$i] = cos($xv); $Q[1][$i] = sin($xv);
            $Q[2][$i] = 0.5 * cos($xv); $Q[3][$i] = 0.5 * sin($xv);
        }
        $e = 0.0;
        // maximum over the output times (as in the other rows of Table 3); it is attained at t = T
        $s->solve($Q, $T, $N, 3.0, 400, function (float $t, array $Qt) use (&$e, $Q, $erfi): void {
            $er = exp(-$t); $ei = -exp(-$t) * $erfi(sqrt($t));
            for ($j = 0; $j < 4; $j += 2) {
                foreach ($Q[$j] as $i => $a) {
                    $b = $Q[$j + 1][$i];
                    $e = max($e, hypot($Qt[$j][$i] - ($er * $a - $ei * $b), $Qt[$j + 1][$i] - ($er * $b + $ei * $a)));
                }
            }
        });
        $errs[] = $e;
    }
    $ord = [null];
    for ($i = 1; $i < count($Ns); $i++) {
        $ord[] = log($errs[$i - 1] / $errs[$i]) / log(2.0);
    }
    $out['table3']['mittag_leffler'] = convRow('Caputo nu = 0.5, linear, Mittag-Leffler (max error over outputs)', ['h' => $Ns, 'err' => $errs, 'order' => $ord]);
    line();
}

// ---------------------------------------------- Table 3: general coefficients
if (in_array('general', $sections, true)) {
    line('== Table 3: general coefficients, V != 0, self-convergence against dt = 7.8125e-4');
    // beta = 1 + 0.8 sin 4t, gamma = 1 + 0.5 cos t, g = 0.02 sin t, V = -0.01 x^2 (1 + 0.5 sin t)
    $L = 60.0; $Nx = 512; $T = 5.0;
    $f = static fn(string $e, array $v) => (new ExprParser($v))->parse($e);
    $solver = static function () use ($L, $Nx, $f): SplitStepSolver {
        $dx = $L / $Nx;
        $x = [];
        for ($i = 0; $i < $Nx; $i++) {
            $x[] = ($i - intdiv($Nx, 2)) * $dx;
        }
        return new SplitStepSolver($L, $Nx, 2.0,
            ExprCompiler::timeFunction($f('1 + 0.8 sin(4t)', ['t'])),
            ExprCompiler::timeFunction($f('1 + 0.5 cos(t)', ['t'])),
            ExprCompiler::timeFunction($f('0.02 sin(t)', ['t'])),
            ExprCompiler::fieldFunction($f('-0.01 x^2 (1 + 0.5 sin(t))', ['x', 't']), $x));
    };
    $sol = [
        ['eta' => 0.5, 'xi' => -0.25, 'x0' => -8.0, 'p1re' => 1.0, 'p1im' => 0.0, 'p2re' => 1.0, 'p2im' => 0.0],
        ['eta' => 0.6, 'xi' => 0.25, 'x0' => 8.0, 'p1re' => 1.0, 'p1im' => 0.0, 'p2re' => 0.0, 'p2im' => 0.0],
    ];
    $s0 = $solver();
    $Q0 = solitonField($s0->x, $sol);
    $nSave = 50;                    // output every 0.1: common to all levels
    $collect = static function (SplitStepSolver $s, float $dt) use ($Q0, $T, $nSave): array {
        $frames = [];
        $s->solve($Q0, $T, $dt, $nSave, static function (float $t, array $Q) use (&$frames): void {
            $frames[sprintf('%.6f', $t)] = $Q;
        });
        return $frames;
    };
    $ref = $collect($solver(), 7.8125e-4);
    $dts = [0.05, 0.025, 0.0125, 0.00625];
    $errs = [];
    foreach ($dts as $dt) {
        $fr = $collect($solver(), $dt);
        $e = 0.0;
        foreach ($fr as $t => $Q) {
            $R = $ref[$t] ?? throw new RuntimeException("output time $t missing in the reference");
            for ($c = 0; $c < 4; $c += 2) {
                foreach ($Q[$c] as $i => $v) {
                    $e = max($e, hypot($v - $R[$c][$i], $Q[$c + 1][$i] - $R[$c + 1][$i]));
                }
            }
        }
        $errs[] = $e;
    }
    $ord = [null];
    for ($i = 1; $i < count($errs); $i++) {
        $ord[] = log($errs[$i - 1] / $errs[$i]) / log(2.0);
    }
    $out['table3']['general'] = convRow('Classical, general coefficients, V != 0, T = 5', ['h' => $dts, 'err' => $errs, 'order' => $ord]);
    $out['table3']['general']['solitons'] = $sol;
    line();
}

// ------------------------------------- Spatial resolution and domain (alpha = 1.5)
if (in_array('fractional', $sections, true)) {
    line('== Fractional soliton alpha = 1.5, eta = 0.5 (mu = 1): tails, domain size, resolution');
    $prof = static function (float $L, int $Nx, float $alpha = 1.5): array {
        $dx = $L / $Nx;
        $x = [];
        for ($i = 0; $i < $Nx; $i++) {
            $x[] = ($i - intdiv($Nx, 2)) * $dx;
        }
        [$phi, $res, $it] = Solitons::fractional($x, FFT::wavenumbers($Nx, $L), $alpha, 1.0);
        return [$x, $phi, $res, $it];
    };
    [$xr, $pr, $resr] = $prof(320.0, 8192);
    $ampRef = $pr[4096];
    // tail exponent: least-squares slope of log phi against log|x| on 20 <= x <= 60 (reference profile)
    $sx = $sy = $sxx = $sxy = 0.0; $n = 0;
    foreach ($xr as $i => $xv) {
        if ($xv >= 20.0 && $xv <= 60.0) {
            $lx = log($xv); $ly = log($pr[$i]);
            $sx += $lx; $sy += $ly; $sxx += $lx * $lx; $sxy += $lx * $ly; $n++;
        }
    }
    $slope = ($n * $sxy - $sx * $sy) / ($n * $sxx - $sx * $sx);
    line(sprintf('  reference L = 320, Nx = 8192: amplitude %.10f, Petviashvili residual %.1e, tail exponent %.3f (theory -(1+alpha) = -2.5)',
        $ampRef, $resr, $slope));
    $out['fractional'] = ['amp_ref' => $ampRef, 'tail_exponent' => $slope, 'residual_ref' => $resr];
    foreach ([[40.0, 256], [80.0, 512], [160.0, 1024]] as [$L, $Nx]) {
        [$x, $phi, $res] = $prof($L, $Nx);
        $d = abs($phi[intdiv($Nx, 2)] - $ampRef);
        $out['fractional']['domain'][] = ['L' => $L, 'Nx' => $Nx, 'amp_diff' => $d];
        line(sprintf('  L = %5.0f, Nx = %4d (dx = 0.156): |amplitude - reference| = %.2e', $L, $Nx, $d));
    }
    [, $p1] = $prof(80.0, 512);
    [, $p2] = $prof(80.0, 1024);
    $d = abs($p1[256] - $p2[512]);
    $out['fractional']['resolution_L80'] = $d;
    line(sprintf('  L = 80: Nx 512 -> 1024 changes the amplitude by %.1e', $d));
    $dl = $out['fractional']['domain'];
    line(sprintf('  truncation exponents: %.2f (40 -> 80), %.2f (80 -> 160)',
        log($dl[0]['amp_diff'] / $dl[1]['amp_diff']) / log(2), log($dl[1]['amp_diff'] / $dl[2]['amp_diff']) / log(2)));
    line();
}

// ------------------------------------------- Verification near the lower end of alpha
if (in_array('alpha', $sections, true)) {
    line('== Exact fractional soliton near the lower end of the validated range');
    foreach ([[1.2, 80.0, 512], [1.2, 160.0, 1024]] as [$alpha, $L, $Nx]) {
        $r = simulate('space-exact', ['params' => ['alpha' => $alpha, 'L' => $L, 'Nx' => $Nx]]);
        $out['alpha'][] = ['alpha' => $alpha, 'L' => $L, 'Nx' => $Nx] + convRow(sprintf('alpha = %.1f, L = %.0f, Nx = %d', $alpha, $L, $Nx), $r['convergence']);
    }
    line();
}

// ---------------------------------------------------------------- Table 4
if (in_array('table4', $sections, true)) {
    line('== Table 4: Caputo memory (preset caputo-decay: eta 0.6, xi 0.1, x0 5, L = 160, Nx = 512, T = 15)');
    $r = simulate('caputo-decay', ['model' => 'classic', 'params' => ['dt' => 0.01]]);
    [$P, $pk] = powerPeak($r);
    $out['table4']['nu1'] = ['P' => $P, 'peak' => $pk, 'peak0' => $r['diagnostics']['peak'][0]];
    line(sprintf('  nu = 1 (split-step, dt = 0.01): P(T)/P(0) = %.4f, peak = %.4f (initial %.4f)', $P, $pk, $r['diagnostics']['peak'][0]));
    $r = simulate('caputo-decay', ['model' => 'classic', 'params' => ['dt' => 0.01, 'L' => 80, 'Nx' => 256]]);
    [, $pk80] = powerPeak($r);
    $out['table4']['nu1_L80_peak'] = $pk80;
    line(sprintf('  (same with L = 80, Nx = 256: peak = %.4f; radiation re-enters the periodic domain)', $pk80));
    foreach ([0.9, 0.8, 0.7, 0.6] as $nu) {
        $row = ['nu' => $nu];
        foreach ([600, 1200] as $N) {
            $r = simulate('caputo-decay', ['params' => ['nu' => $nu, 'N' => $N]]);
            [$P, $pk] = powerPeak($r);
            $row["N$N"] = ['P' => $P, 'peak' => $pk, 'picard' => [$r['meta']['picardMean'], $r['meta']['picardMax']]];
        }
        $r = simulate('caputo-decay', ['params' => ['nu' => $nu, 'L' => 80, 'Nx' => 256]]);
        [$P, $pk] = powerPeak($r);
        $row['L80_Nx256'] = ['P' => $P, 'peak' => $pk];
        $out['table4']['rows'][] = $row;
        line(sprintf('  nu = %.1f: N=600 %.4f / %.4f [Picard %.1f/%d] | N=1200 %.4f / %.4f [Picard %.1f/%d] | L=80,Nx=256: %.4f / %.4f',
            $nu, $row['N600']['P'], $row['N600']['peak'], $row['N600']['picard'][0], $row['N600']['picard'][1],
            $row['N1200']['P'], $row['N1200']['peak'], $row['N1200']['picard'][0], $row['N1200']['picard'][1], $P, $pk));
    }
    line();
}

// ------------------------------- Table 4: spatial refinement at fixed L (Nx only; the L/Nx test keeps dx)
if (in_array('spatial', $sections, true)) {
    line('== Table 4, spatial refinement at fixed L = 160 (caputo-decay, N = 600): P(T)/P(0) and peak for Nx = 256, 512, 1024');
    foreach ([0.9, 0.8, 0.7, 0.6] as $nu) {
        $row = ['nu' => $nu];
        $txt = [];
        foreach ([256, 512, 1024] as $Nx) {
            $r = simulate('caputo-decay', ['params' => ['nu' => $nu, 'N' => 600, 'Nx' => $Nx]]);
            [$P, $pk] = powerPeak($r);
            $row["Nx$Nx"] = ['P' => $P, 'peak' => $pk];
            $txt[] = sprintf('Nx=%d %.6f / %.6f', $Nx, $P, $pk);
        }
        $out['spatial'][] = $row;
        line(sprintf('  nu = %.1f: %s', $nu, implode(' | ', $txt)));
    }
    line();
}

// --------------------------------------------- L1 scheme as nu -> 1 (numerical dissipation)
if (in_array('nu1', $sections, true)) {
    line('== L1 scheme close to nu = 1: P(T)/P(0) under refinement of N (caputo-decay)');
    foreach ([0.95, 0.99] as $nu) {
        $row = ['nu' => $nu];
        $txt = [];
        foreach ([600, 1200, 2400] as $N) {
            $r = simulate('caputo-decay', ['params' => ['nu' => $nu, 'N' => $N]]);
            [$P] = powerPeak($r);
            $row["N$N"] = $P;
            $txt[] = sprintf('N=%d %.4f', $N, $P);
        }
        $out['nu1'][] = $row;
        line(sprintf('  nu = %.2f: %s', $nu, implode(' | ', $txt)));
    }
    line();
}

// ------------------------------------- Nonlinear solve: Anderson versus plain Picard
if (in_array('solver', $sections, true)) {
    line('== Nonlinear solve of the L1 scheme: Anderson (depth 10) versus plain Picard');
    $cases = [
        ['caputo-decay, nu = 0.1, T = 1, Nx = 128', 'caputo-decay', ['params' => ['nu' => 0.1, 'T' => 1, 'Nx' => 128]]],
        ['caputo-decay, nu = 0.05, T = 1, Nx = 128', 'caputo-decay', ['params' => ['nu' => 0.05, 'T' => 1, 'Nx' => 128]]],
        ['caputo-mms, amplitude 1 (eta = 0.5), nu = 0.6', 'caputo-mms',
            ['solitons' => [['eta' => 0.5, 'xi' => -0.25, 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0.5, 'p2im' => -0.5]],
             'convergence' => false]],
    ];
    foreach ($cases as [$label, $id, $over]) {
        $row = ['case' => $label];
        $txt = [];
        foreach ([64, 256, 1024] as $N) {
            $o = $over;
            $o['params'] = array_merge($o['params'] ?? [], ['N' => $N]);
            if ($id === 'caputo-mms' && $N !== 256) {
                continue;
            }
            $r = simulate($id, $o);
            $v = $id === 'caputo-mms' ? (float)$r['exact']['err_max'] : powerPeak($r)[0];
            $row["N$N"] = ['value' => $v, 'picard' => [$r['meta']['picardMean'], $r['meta']['picardMax']]];
            $txt[] = sprintf('N=%d %s %.6g [%.1f/%d]', $N, $id === 'caputo-mms' ? 'err' : 'P(T)/P(0)', $v,
                $r['meta']['picardMean'], $r['meta']['picardMax']);
        }
        $o = $over;
        $o['params'] = array_merge($o['params'] ?? [], ['N' => $id === 'caputo-mms' ? 256 : 64]);
        $o['anderson'] = 0;
        try {
            simulate($id, $o);
            $row['picard'] = 'converged';
        } catch (Throwable $e) {
            preg_match('/at step (\d+)/', $e->getMessage(), $m);
            $row['picard'] = 'fails at step ' . ($m[1] ?? '?');
        }
        $out['solver'][] = $row;
        line(sprintf('  %s: %s | plain Picard (N=%d): %s', $label, implode(' | ', $txt), $o['params']['N'], $row['picard']));
    }
    line();
}

// ----------------------------------------------------------------- Collision
if (in_array('collision', $sections, true)) {
    line('== Shape-changing collision (preset classic-collision), powers per soliton divided by exp(2 int g)');
    $r = simulate('classic-collision');
    $x = $r['x'];
    $dx = $x[1] - $x[0];
    $a1 = $r['final']['q1']['abs'];
    $a2 = $r['final']['q2']['abs'];
    $gain = end($r['diagnostics']['gainFactor']);
    $I = [];
    foreach ($a1 as $i => $v) {
        $I[$i] = $v * $v + $a2[$i] * $a2[$i];
    }
    // the two highest local maxima, and the intensity minimum between them
    $pk = [];
    for ($i = 1; $i < count($I) - 1; $i++) {
        if ($I[$i] > $I[$i - 1] && $I[$i] >= $I[$i + 1]) {
            $pk[$i] = $I[$i];
        }
    }
    arsort($pk);
    $two = array_slice(array_keys($pk), 0, 2);
    sort($two);
    $cut = $two[0];
    for ($i = $two[0]; $i <= $two[1]; $i++) {
        if ($I[$i] < $I[$cut]) {
            $cut = $i;
        }
    }
    $pw = static function (array $a, int $from, int $to) use ($dx, $gain): float {
        $s = 0.0;
        for ($i = $from; $i < $to; $i++) {
            $s += $a[$i] * $a[$i];
        }
        return $s * $dx / $gain;
    };
    $left = [$pw($a1, 0, $cut), $pw($a2, 0, $cut)];
    $right = [$pw($a1, $cut, count($I)), $pw($a2, $cut, count($I))];
    // soliton 1 starts on the left and moves right (velocity -4 xi with xi = -0.25), so it ends on the right
    $out['collision'] = ['soliton1' => $right, 'soliton2' => $left, 'gain' => $gain];
    line(sprintf('  soliton 1: (1.000, 1.000) -> (%.3f, %.3f), total %.3f', $right[0], $right[1], $right[0] + $right[1]));
    line(sprintf('  soliton 2: (2.400, 0.000) -> (%.3f, %.3f), total %.3f', $left[0], $left[1], $left[0] + $left[1]));
    line();
}

// -------------------------------------------------------------------- Timing
if (in_array('timing', $sections, true)) {
    line('== Table 5: run time per preset (runtime_total, median of 3), ' . engine());
    foreach (Presets::all() as $p) {
        $t = [];
        for ($k = 0; $k < 3; $k++) {
            $t[] = simulate($p['id'])['meta']['runtime_total'];
        }
        sort($t);
        $out['timing'][$p['id']] = ['median' => $t[1], 'min' => $t[0], 'max' => $t[2]];
        line(sprintf('  %-62s %6.2f s  (range %.2f-%.2f)', $p['title'], $t[1], $t[0], $t[2]));
    }
    line();
}

if ($jsonFile !== null) {
    file_put_contents($jsonFile, json_encode($out, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE | JSON_PRESERVE_ZERO_FRACTION) . "\n");
    line("JSON written to $jsonFile");
}
