<?php
declare(strict_types=1);

/**
 * Sweep of the a posteriori checks of the exact mode (Section 3.1 of the article).
 *
 *   php tests/sweep.php [--json=FILE] > results/sweep.txt
 *
 * Runs classic-exact with beta = gamma = 1, g = V = 0 and one soliton c = (1, 0.6 + 0.4i), x0 = 0, over two grids:
 *   A (four-level studies): eta in {0.5, 1, 1.5, 2, 3}, xi in {0, -0.5, -2, -4}, Nx in {32, 64, 128, 256},
 *     L in {40, 80}, (T, dt) in {(1, 0.0025), (1, 0.01), (4, 0.01)};
 *   B (two-level studies): eta in {0.15, 0.25, 0.35, 0.5}, xi in {0, -0.5, -1}, Nx in {16, 32, 64, 128},
 *     L in {40, 80}, T in {2, 4, 8, 16}, dt = T/25 (dt = T/12 leaves one level and gives no order).
 * Configurations that exact mode refuses (overlapping images, carrier beyond the grid) are counted and skipped.
 * One line per admitted configuration lists the levels, the errors, the observed orders (null = withheld) and
 * the level lists saturated / preasymptotic / untested; a summary closes each grid. Every observed order
 * that is published is compared with the expected order 2.
 *
 * This script is a command-line tool: it refuses to run inside a web server.
 */

if (PHP_SAPI !== 'cli') {
    http_response_code(404);
    exit;
}

require dirname(__DIR__) . '/src/bootstrap.php';
require __DIR__ . '/lib_problem.php';

use ManakovLab\Numerics\Simulation;

$jsonFile = null;
foreach (array_slice($argv, 1) as $a) {
    if (str_starts_with($a, '--json=')) {
        $jsonFile = substr($a, 7);
    }
}

$coeffs = ['beta' => '1', 'gamma' => '1', 'g' => '0', 'V' => '0'];
$grids = [
    'A' => (static function (): \Generator {
        foreach ([0.5, 1.0, 1.5, 2.0, 3.0] as $eta) foreach ([0.0, -0.5, -2.0, -4.0] as $xi)
        foreach ([32, 64, 128, 256] as $Nx) foreach ([40.0, 80.0] as $L) foreach ([[1.0, 0.0025], [1.0, 0.01], [4.0, 0.01]] as [$T, $dt]) {
            yield compact('eta', 'xi', 'Nx', 'L', 'T', 'dt');
        }
    })(),
    'B' => (static function (): \Generator {
        foreach ([0.15, 0.25, 0.35, 0.5] as $eta) foreach ([0.0, -0.5, -1.0] as $xi)
        foreach ([16, 32, 64, 128] as $Nx) foreach ([40.0, 80.0] as $L) foreach ([2.0, 4.0, 8.0, 16.0] as $T) {
            $dt = $T / 25;
            yield compact('eta', 'xi', 'Nx', 'L', 'T', 'dt');
        }
    })(),
];

$fmt = static fn(array $a, string $f): string => '[' . implode(', ', array_map(
    static fn($v) => $v === null ? 'null' : sprintf($f, $v), $a)) . ']';
$out = ['version' => MANAKOV_LAB_VERSION, 'engine' => 'PHP ' . PHP_VERSION, 'grids' => []];
printf("Manakov Lab sweep of the exact-mode checks (classic-exact, beta = gamma = 1, g = V = 0)\n");
foreach ($grids as $name => $gen) {
    $tot = $adm = $pub = 0;
    $dev = 0.0;
    $within = ['0.05' => 0, '0.1' => 0];
    $cnt = ['saturated' => 0, 'preasymptotic' => 0, 'untested' => 0, 'belowFloor' => 0];
    $levels = [];
    $rows = [];
    printf("\n== grid %s ==\n", $name);
    foreach ($gen as $c) {
        $tot++;
        $over = ['params' => ['L' => $c['L'], 'Nx' => $c['Nx'], 'T' => $c['T'], 'dt' => $c['dt']], 'coeffs' => $coeffs,
            'solitons' => [['eta' => $c['eta'], 'xi' => $c['xi'], 'x0' => 0, 'p1re' => 1, 'p1im' => 0, 'p2re' => 0.6, 'p2im' => 0.4]],
            'exact' => true, 'convergence' => true];
        try {
            $r = (new Simulation(problemFrom('classic-exact', $over)))->run();
        } catch (\Throwable $e) {
            continue;                                   // refused by exact mode
        }
        $adm++;
        $cv = $r['convergence'];
        $n = count($cv['h']);
        $levels[$n] = ($levels[$n] ?? 0) + 1;
        foreach (['saturated', 'preasymptotic', 'untested'] as $k) {
            $cnt[$k] += count($cv[$k]);
        }
        foreach ($cv['order'] as $i => $o) {
            if ($o !== null) {
                $pub++;
                $dev = max($dev, abs($o - 2.0));
                $within['0.05'] += abs($o - 2.0) <= 0.05 ? 1 : 0;
                $within['0.1'] += abs($o - 2.0) <= 0.1 ? 1 : 0;
            } elseif ($i > 0 && $cv['err'][$i] !== null && $cv['err'][$i] <= 10 * $cv['refFloor']) {
                $cnt['belowFloor']++;
            }
        }
        printf("eta=%g xi=%g Nx=%d L=%g T=%g dt=%g | h=%s err=%s order=%s | floor=%.2e saturated=%d preasymptotic=%d untested=%d\n",
            $c['eta'], $c['xi'], $c['Nx'], $c['L'], $c['T'], $c['dt'], $fmt($cv['h'], '%.6g'), $fmt($cv['err'], '%.4e'),
            $fmt($cv['order'], '%.4f'), $cv['refFloor'], count($cv['saturated']), count($cv['preasymptotic']), count($cv['untested']));
        $rows[] = $c + ['h' => $cv['h'], 'err' => $cv['err'], 'order' => $cv['order'], 'refFloor' => $cv['refFloor'],
            'saturated' => $cv['saturated'], 'preasymptotic' => $cv['preasymptotic'], 'untested' => $cv['untested']];
    }
    ksort($levels);
    $sum = ['configurations' => $tot, 'admitted' => $adm, 'refused' => $tot - $adm, 'studiesByLevels' => $levels,
        'publishedOrders' => $pub, 'maxDeviationFrom2' => $pub ? $dev : null, 'publishedWithin' => $within,
        'withheld' => $cnt];
    printf("-- grid %s: %d configurations, %d admitted (%d refused); studies by number of levels %s; %d observed orders "
        . "published (%d within 0.05 of 2, %d within 0.1), max |order - 2| = %s; withheld: %d saturated, %d preasymptotic, %d untested (two levels), "
        . "%d below ten times the floor\n", $name, $tot, $adm, $tot - $adm, json_encode($levels), $pub,
        $within['0.05'], $within['0.1'], $pub ? sprintf('%.4f', $dev) : '-', $cnt['saturated'], $cnt['preasymptotic'], $cnt['untested'], $cnt['belowFloor']);
    $out['grids'][$name] = ['summary' => $sum, 'rows' => $rows];
}
if ($jsonFile !== null) {
    file_put_contents($jsonFile, json_encode($out, JSON_PRETTY_PRINT) . "\n");
}
