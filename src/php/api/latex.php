<?php
declare(strict_types=1);

/**
 * POST /api/latex.php   expr=1%2B0.5cos(2t)&vars=t   (form-encoded; GET with the same parameters is also accepted)
 * -> {"ok":true,"latex":"1 + 0.5\\,\\cos\\left(2\\,t\\right)","canonical":"(1+(0.5*cos((2*t))))","uses":["t"]}
 */

use ManakovLab\ExprException;
use ManakovLab\ExprParser;
use ManakovLab\Http;

require dirname(__DIR__) . '/src/bootstrap.php';

// POST keeps formulas out of URLs (logs, length limits, firewall rules on query strings).
$in = $_SERVER['REQUEST_METHOD'] === 'POST' ? $_POST : $_GET;
$expr = (string)($in['expr'] ?? '');
if (strlen($expr) > 2000) {
    Http::json(['ok' => false, 'error' => 'Expression too long (at most 2000 characters).', 'position' => 2000], 200);
}
$vars = ($in['vars'] ?? 't') === 'xt' ? ['x', 't'] : ['t'];

try {
    $ast = (new ExprParser($vars))->parse($expr);
    Http::json([
        'ok' => true,
        'latex' => ExprParser::toLatex($ast),
        'canonical' => ExprParser::toCanonical($ast),
        'uses' => ExprParser::usedVariables($ast),
    ]);
} catch (ExprException $e) {
    Http::json(['ok' => false, 'error' => $e->getMessage(), 'position' => $e->position], 200);
}
