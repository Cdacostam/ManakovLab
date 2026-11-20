<?php
declare(strict_types=1);

/**
 * POST /api/run.php   (Content-Type: application/json)
 *
 * {
 *   "model": "classic" | "space" | "caputo",
 *   "params": {"L":80,"Nx":512,"T":10,"dt":0.005,"alpha":1.5,"nu":0.8,"N":1500,"r":"auto"},
 *   "coeffs": {"beta":"1+0.5cos(2t)","gamma":"1","g":"0","V":"0"},
 *   "solitons": [{"eta":0.5,"xi":0,"x0":0,"p1re":1,"p1im":0,"p2re":0,"p2im":0}],
 *   "exact": false, "convergence": false
 * }
 */

use ManakovLab\ExprException;
use ManakovLab\ExprParser;
use ManakovLab\Http;
use ManakovLab\ConcurrencyGate;
use ManakovLab\Numerics\Simulation;
use ManakovLab\Numerics\SimulationException;

require dirname(__DIR__) . '/src/bootstrap.php';

if (($_SERVER['REQUEST_METHOD'] ?? '') !== 'POST') {
    header('Allow: POST');
    Http::error('Use POST with a JSON body.', 405);
}
$raw = file_get_contents('php://input', false, null, 0, 65536);
try {
    $in = json_decode((string)$raw, true, 16, JSON_THROW_ON_ERROR);
} catch (JsonException) {
    Http::error('The request body is not valid JSON.');
}
if (!is_array($in)) {
    Http::error('The request body must be a JSON object.');
}

$cfg = manakov_config();

// ------------------------------------------------------------------ helpers
$num = static function (array $src, string $key, float $default, float $lo, float $hi, string $label) {
    $v = $src[$key] ?? $default;
    if (!is_numeric($v)) {
        Http::error("$label must be a number.", 422, ['field' => $key]);
    }
    $v = (float)$v;
    if (!is_finite($v) || $v < $lo || $v > $hi) {
        Http::error("$label must be between $lo and $hi.", 422, ['field' => $key]);
    }
    return $v;
};

// ------------------------------------------------------------------- model
$model = $in['model'] ?? 'classic';
if (!in_array($model, ['classic', 'space', 'caputo'], true)) {
    Http::error('Unknown model.', 422, ['field' => 'model']);
}
$p = is_array($in['params'] ?? null) ? $in['params'] : [];
$params = [
    'L' => $num($p, 'L', 80, 5, 1000, 'Domain length L'),
    'Nx' => (int)$num($p, 'Nx', 512, 64, 2048, 'Grid points Nx'),
    'T' => $num($p, 'T', 10, 0.01, 200, 'Final time T'),
];
$notes = [];
$nxIn = $params['Nx'];
$params['Nx'] = 2 ** (int)round(log($params['Nx'], 2));
if ($params['Nx'] !== $nxIn) {
    $notes[] = "Nx = $nxIn is not a power of two (required by the FFT): Nx = {$params['Nx']} was used.";
}

if ($model === 'caputo') {
    $params['nu'] = $num($p, 'nu', 0.8, 0.05, 1.0, 'Caputo order ν');
    $params['N'] = (int)$num($p, 'N', 500, 16, 3000, 'Time steps N');
    if ($params['Nx'] > 1024) {
        Http::error('For the Caputo model Nx must be at most 1024.', 422, ['field' => 'Nx']);
    }
    $r = $p['r'] ?? 'auto';
    $params['r'] = ($r === 'auto' || $r === '' || $r === null) ? 'auto' : $num($p, 'r', 1, 1, 10, 'Grading exponent r');
    if ($params['r'] === 'auto' && $params['nu'] < 1.0 && (2 - $params['nu']) / $params['nu'] > 10) {
        $notes[] = sprintf('The optimal grading r = (2 - nu)/nu = %.3g exceeds the admissible maximum 10: r = 10 was used, '
            . 'so the expected order is min(r nu, 2 - nu) = %.3g.', (2 - $params['nu']) / $params['nu'], 10 * $params['nu']);
    }
} else {
    $params['dt'] = $num($p, 'dt', 0.005, 1e-5, 1.0, 'Time step dt');
    if ($params['dt'] > $params['T']) {
        Http::error('The time step dt must not exceed the final time T.', 422, ['field' => 'dt']);
    }
    if ($model === 'space') {
        if (!array_key_exists('alpha', $p)) {
            Http::error('The space-fractional model needs the order α.', 422, ['field' => 'alpha']);
        }
        $params['alpha'] = $num($p, 'alpha', 1.5, 0.5, 2.0, 'Fractional order α');
    }
}

$exact = (bool)($in['exact'] ?? false);
$convergence = $exact && (bool)($in['convergence'] ?? false);

// ------------------------------------------------------------ coefficients
$c = is_array($in['coeffs'] ?? null) ? $in['coeffs'] : [];
$coeffs = [];
$latex = [];
$asts = [];
foreach (['beta' => ['t'], 'gamma' => ['t'], 'g' => ['t'], 'V' => ['x', 't']] as $name => $vars) {
    $src = trim((string)($c[$name] ?? ($name === 'beta' || $name === 'gamma' ? '1' : '0')));
    if ($src === '' && ($name === 'beta' || $name === 'gamma')) {
        // an empty dispersion or nonlinearity would silently change the physics (beta = 0: no dispersion)
        Http::error("Coefficient $name must not be empty.", 422, ['field' => $name]);
    }
    try {
        $ast = (new ExprParser($vars))->parse($src);
    } catch (ExprException $e) {
        Http::error("Coefficient $name: " . $e->getMessage(), 422, ['field' => $name, 'position' => $e->position]);
    }
    $asts[$name] = $ast;
    $coeffs[$name] = ExprParser::toCanonical($ast);
    $latex[$name] = ExprParser::toLatex($ast);
}

// ---------------------------------------------------------------- solitons
$sIn = is_array($in['solitons'] ?? null) ? array_values($in['solitons']) : [];
if (count($sIn) < 1 || count($sIn) > 3) {
    Http::error('Provide between 1 and 3 initial solitons.', 422, ['field' => 'solitons']);
}
$half = $params['L'] / 2;
$solitons = [];
foreach ($sIn as $i => $s) {
    if (!is_array($s)) {
        Http::error('Invalid soliton definition.', 422);
    }
    $k = $i + 1;
    $row = [
        'eta' => $num($s, 'eta', 0.5, 0.01, 5, "Soliton $k: η"),
        'xi' => $num($s, 'xi', 0, -5, 5, "Soliton $k: ξ"),
        'x0' => $num($s, 'x0', 0, -$half, $half, "Soliton $k: x₀"),
    ];
    foreach (['p1re', 'p1im', 'p2re', 'p2im'] as $f) {
        $row[$f] = $num($s, $f, 0, -10, 10, "Soliton $k: polarization");
    }
    if (abs($row['p1re']) + abs($row['p1im']) + abs($row['p2re']) + abs($row['p2im']) == 0.0) {
        Http::error("Soliton $k: the polarization vector cannot be zero.", 422, ['field' => 'solitons']);
    }
    $solitons[] = $row;
}

// ------------------------------------------------------------- work budget
$Nx = $params['Nx'];
if ($model === 'caputo') {
    $work = $params['N'] ** 2 * $Nx * ($convergence ? 1.35 : 1.0);
    if ($work > $cfg['MANAKOV_MAX_WORK_CAPUTO']) {
        Http::error('This Caputo run is too large for the server budget (cost ∝ N²·Nx). Reduce N or Nx.', 422);
    }
} else {
    $steps = $params['T'] / $params['dt'];
    $work = $Nx * log($Nx, 2) * $steps * ($convergence ? 1.9 : 1.0);
    if ($work > $cfg['MANAKOV_MAX_WORK_SPLIT']) {
        Http::error('This run is too large for the server budget (cost ∝ Nx·log Nx·T/dt). Increase dt or reduce Nx / T.', 422);
    }
}

// --------------------------------------------------------------------- run
@set_time_limit($cfg['MANAKOV_TIMEOUT'] + 30);
@ini_set('memory_limit', (string)$cfg['MANAKOV_MEMORY_LIMIT']);
ignore_user_abort(false);

$gate = new ConcurrencyGate($cfg['MANAKOV_MAX_CONCURRENT']);
if (!$gate->acquire(30)) {
    Http::error('The server is busy with other simulations. Please try again in a moment.', 429);
}

$problem = [
    'model' => $model,
    'L' => $params['L'], 'Nx' => $params['Nx'], 'T' => $params['T'],
    'dt' => $params['dt'] ?? null,
    'N' => $params['N'] ?? null,
    'alpha' => $params['alpha'] ?? 2.0,
    'nu' => $params['nu'] ?? 1.0,
    'ast' => $asts, 'solitons' => $solitons,
    'exact' => $exact, 'convergence' => $convergence,
];
if ($model === 'caputo') {
    $problem['r'] = $params['r'] === 'auto' ? min(10.0, (2 - $params['nu']) / $params['nu']) : (float)$params['r'];
}
$problem['warnings'] = $notes;

try {
    $result = (new Simulation($problem, microtime(true) + $cfg['MANAKOV_TIMEOUT']))->run();
} catch (SimulationException $e) {
    Http::error($e->getMessage(), 422, $e->field !== null ? ['field' => $e->field] : []);
} catch (\Throwable $e) {
    error_log('[manakov-lab] ' . $e::class . ': ' . $e->getMessage() . ' @ ' . $e->getFile() . ':' . $e->getLine());
    Http::error('Internal solver error. See the server error log for details.', 500);
} finally {
    $gate->release();
}

Http::json([
    'ok' => true,
    'result' => $result,
    'input' => ['model' => $model, 'params' => $params, 'latex' => $latex, 'canonical' => $coeffs,
        'solitons' => $solitons, 'exact' => $exact, 'convergence' => $convergence],
]);
