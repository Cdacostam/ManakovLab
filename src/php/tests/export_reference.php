<?php
declare(strict_types=1);

/**
 * Runs the cross-check cases (tests/crosscheck_cases.json) and writes the PHP diagnostics as JSON,
 * for comparison with the Python implementation (src/python/crosscheck.py).
 *   php tests/export_reference.php OUTPUT.json
 */

if (PHP_SAPI !== 'cli') {
    http_response_code(404);
    exit;
}

require dirname(__DIR__) . '/src/bootstrap.php';
require __DIR__ . '/lib_problem.php';

use ManakovLab\Numerics\Simulation;

$cases = json_decode((string)file_get_contents(__DIR__ . '/crosscheck_cases.json'), true, 32, JSON_THROW_ON_ERROR);
$out = ['version' => MANAKOV_LAB_VERSION, 'engine' => 'PHP ' . PHP_VERSION, 'cases' => []];
foreach ($cases as $c) {
    $r = (new Simulation(problemFrom($c['id'], $c['over'])))->run();
    $out['cases'][] = [
        'id' => $c['id'], 'over' => $c['over'],
        't' => $r['t'], 'diagnostics' => $r['diagnostics'],
        'final_abs' => [$r['final']['q1']['abs'], $r['final']['q2']['abs']],
        'final' => [[$r['final']['q1']['re'], $r['final']['q1']['im']], [$r['final']['q2']['re'], $r['final']['q2']['im']]],
        'err_t' => $r['exact']['err_t'] ?? null,
        'picard' => [$r['meta']['picardMean'], $r['meta']['picardMax']],
        'picardSteps' => $r['meta']['picardSteps'],
        'warnings' => $r['warnings'],
    ];
    fprintf(STDERR, "%s done\n", $c['id']);
}
file_put_contents($argv[1] ?? 'php_reference.json', json_encode($out, JSON_PRESERVE_ZERO_FRACTION) . "\n");
