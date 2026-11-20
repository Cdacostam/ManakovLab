<?php
declare(strict_types=1);

const MANAKOV_LAB_VERSION = '1.0.0';

require_once __DIR__ . '/ExprParser.php';
require_once __DIR__ . '/Presets.php';
require_once __DIR__ . '/Http.php';
require_once __DIR__ . '/ConcurrencyGate.php';
require_once __DIR__ . '/Numerics/Gamma.php';
require_once __DIR__ . '/Numerics/FFT.php';
require_once __DIR__ . '/Numerics/Integral.php';
require_once __DIR__ . '/Numerics/ExprCompiler.php';
require_once __DIR__ . '/Numerics/Solitons.php';
require_once __DIR__ . '/Numerics/SplitStepSolver.php';
require_once __DIR__ . '/Numerics/CaputoSolver.php';
require_once __DIR__ . '/Numerics/Simulation.php';

/** @return array<string,mixed> */
function manakov_config(): array
{
    static $cfg = null;
    return $cfg ??= require dirname(__DIR__) . '/config.php';
}
