<?php
declare(strict_types=1);

// Shared by tests/reproduce.php and tests/run_tests.php (command line only).
if (PHP_SAPI !== 'cli') {
    http_response_code(404);
    exit;
}

use ManakovLab\ExprParser;
use ManakovLab\Presets;

/** Build a validated problem from a preset id and overrides (same path as api/run.php). */
function problemFrom(string $id, array $over = []): array
{
    $p = Presets::byId($id) ?? throw new RuntimeException("unknown preset $id");
    foreach (['params', 'coeffs'] as $k) {
        if (isset($over[$k])) {
            $p[$k] = array_merge($p[$k], $over[$k]);
        }
    }
    foreach (['solitons', 'exact', 'convergence', 'model'] as $k) {
        if (array_key_exists($k, $over)) {
            $p[$k] = $over[$k];
        }
    }
    $ast = [];
    foreach (['beta' => ['t'], 'gamma' => ['t'], 'g' => ['t'], 'V' => ['x', 't']] as $n => $v) {
        $ast[$n] = (new ExprParser($v))->parse((string)$p['coeffs'][$n]);
    }
    $q = $p['params'];
    $pr = [
        'model' => $p['model'], 'L' => (float)$q['L'], 'Nx' => (int)$q['Nx'], 'T' => (float)$q['T'],
        'dt' => isset($q['dt']) ? (float)$q['dt'] : null, 'N' => isset($q['N']) ? (int)$q['N'] : null,
        'alpha' => (float)($q['alpha'] ?? 2.0), 'nu' => (float)($q['nu'] ?? 1.0),
        'ast' => $ast, 'solitons' => $p['solitons'], 'exact' => (bool)$p['exact'], 'convergence' => (bool)$p['convergence'],
    ];
    if (isset($over['anderson'])) {
        $pr['anderson'] = (int)$over['anderson'];
    }
    if ($p['model'] === 'caputo') {
        $r = $q['r'] ?? 'auto';
        $pr['r'] = $r !== 'auto' ? (float)$r
            : ($pr['nu'] > 0.0 && $pr['nu'] < 1.0 ? min(10.0, (2 - $pr['nu']) / $pr['nu']) : 1.0);
    }
    return $pr;
}

