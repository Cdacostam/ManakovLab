<?php
declare(strict_types=1);

/**
 * Builds the static preview shown on the home page (assets/img/collision-preview.svg):
 * a waterfall plot of |q1|^2 and |q2|^2 for the preset classic-collision, computed by the
 * same solver the laboratory uses. Colours come from the page's CSS variables, so the
 * SVG follows the light and dark themes when it is inlined by index.php.
 *   php tests/make_home_preview.php
 */

if (PHP_SAPI !== 'cli') {
    http_response_code(404);
    exit;
}

require dirname(__DIR__) . '/src/bootstrap.php';
require __DIR__ . '/lib_problem.php';

use ManakovLab\Numerics\Simulation;

$r = (new Simulation(problemFrom('classic-collision')))->run();
$S = $r['surface'];
$xs = $S['x'];
$ts = $S['t'];

// geometry (SVG user units)
$panelW = 300; $panelH = 230; $gap = 40; $left = 34; $top = 30; $bottom = 38;
$W = $left + 2 * $panelW + $gap + 10;
$H = $top + $panelH + $bottom;
$xMin = -40.0; $xMax = 40.0;          // window shown (domain is [-50, 50))
$slices = 30;                          // time slices drawn
$rise = 5.4;                           // vertical offset between consecutive slices
$amp = $panelH - ($slices - 1) * $rise - 6;
$peak = 0.0;
foreach (['I1', 'I2'] as $f) {
    foreach ($S[$f] as $row) {
        $peak = max($peak, max($row));
    }
}

$fmt = static fn(float $v): string => rtrim(rtrim(sprintf('%.1f', $v), '0'), '.');
$out = [];
$out[] = sprintf('<svg class="preview-svg" viewBox="0 0 %d %d" role="img" aria-labelledby="pv-title pv-desc" xmlns="http://www.w3.org/2000/svg">', $W, $H);
$out[] = '<title id="pv-title">Shape-changing collision of two vector solitons</title>';
$out[] = '<desc id="pv-desc">Waterfall plot of the intensities of the two components over time for the preset classic-collision: '
    . 'the solitons approach, collide and exchange energy between the components.</desc>';

$nt = count($ts);
foreach (['I1' => ['var(--series-1)', '|q₁|²'], 'I2' => ['var(--series-2)', '|q₂|²']] as $f => [$col, $label]) {
    $x0 = $left + ($f === 'I1' ? 0 : $panelW + $gap);
    $out[] = sprintf('<g><text x="%s" y="20" class="pv-label" fill="%s">%s</text>', $fmt($x0), $col, $label);
    // back to front: later times are drawn higher and first, so earlier slices hide them
    for ($k = $slices - 1; $k >= 0; $k--) {
        $j = (int)round($k * ($nt - 1) / ($slices - 1));
        $base = $top + $panelH - $k * $rise;
        $pts = [];
        foreach ($xs as $i => $x) {
            if ($x < $xMin || $x > $xMax || $i % 2 === 1) {   // every second grid point is enough at this size
                continue;
            }
            $px = $x0 + ($x - $xMin) / ($xMax - $xMin) * $panelW;
            $py = $base - $S[$f][$j][$i] / $peak * $amp;
            $pts[] = $fmt($px) . ',' . $fmt($py);
        }
        $first = explode(',', $pts[0]);
        $last = explode(',', end($pts));
        $line = 'M' . implode('L', $pts);
        $out[] = sprintf('<path d="%sL%s,%sL%s,%sZ" fill="var(--surface)" stroke="none"/>', $line, $last[0], $fmt($base), $first[0], $fmt($base));
        $out[] = sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linejoin="round"/>', $line, $col, $k === 0 ? '1.6' : '1.1');
    }
    // axes hints
    $yb = $top + $panelH + 18;
    $out[] = sprintf('<text x="%s" y="%s" class="pv-axis" text-anchor="middle">x</text>', $fmt($x0 + $panelW / 2), $fmt($yb + 4));
    $out[] = sprintf('<text x="%s" y="%s" class="pv-axis">−40</text><text x="%s" y="%s" class="pv-axis" text-anchor="end">40</text>',
        $fmt($x0), $fmt($yb + 4), $fmt($x0 + $panelW), $fmt($yb + 4));
    $out[] = '</g>';
}
$out[] = sprintf('<text x="%s" y="%s" class="pv-axis" transform="rotate(-90 %s %s)" text-anchor="middle">time t → %s</text>',
    $fmt(12), $fmt($top + $panelH / 2), $fmt(12), $fmt($top + $panelH / 2), $fmt(end($ts)));
$out[] = '</svg>';

$file = dirname(__DIR__) . '/assets/img/collision-preview.svg';
@mkdir(dirname($file), 0755, true);
file_put_contents($file, implode("\n", $out) . "\n");
fprintf(STDERR, "wrote %s (%d bytes), peak |q|^2 = %.4f, T = %s\n", $file, filesize($file), $peak, end($ts));
