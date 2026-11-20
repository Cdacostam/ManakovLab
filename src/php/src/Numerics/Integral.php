<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

final class Integral
{
    /** Simpson's rule on [a, b]. */
    public static function simpson(\Closure $f, float $a, float $b): float
    {
        return ($b - $a) / 6.0 * ($f($a) + 4.0 * $f(0.5 * ($a + $b)) + $f($b));
    }

    /**
     * Antiderivative F(t) = int_0^t f(s) ds, tabulated with composite Simpson on
     * m panels of [0, T] and completed with one local Simpson panel, so each
     * evaluation costs three calls of f and is accurate to O(h^5).
     */
    public static function antiderivative(\Closure $f, float $T, int $m = 4096): \Closure
    {
        return self::closure($f, self::table($f, $T, $m), $T / $m, $m);
    }

    /**
     * Antiderivative with a controlled table: starting from m0 panels, the number of panels is doubled
     * until two successive tables agree at their common nodes to tol * max(1, max|F|), or until mMax.
     * The coarser of the two agreeing tables is used, and it must also agree with a composite three-point
     * Gauss-Legendre table on the same panels, whose nodes are not commensurate with Simpson's (this catches a
     * coefficient aliased at both Simpson grids, e.g. cos(32768 pi t) on [0, 1]). A coefficient that varies on a
     * scale shorter than the panel can otherwise give a wrong reference (e.g. beta = 300 sin(40 pi t) on [0, 100]).
     * Both tests are estimates of the quadrature error, not a certification for every formula.
     * @return array{0: \Closure, 1: array{m: int, resolved: bool, max: float, min: float, maxAbs: float}}
     */
    public static function antiderivativeAdaptive(\Closure $f, float $T, float $tol = 1e-8, int $m0 = 4096,
                                                  int $mMax = 1048576): array
    {
        $m = $m0;
        $F = self::table($f, $T, $m);
        $resolved = false;
        while (true) {
            if (2 * $m > $mMax) {
                break;
            }
            $F2 = self::table($f, $T, 2 * $m);
            $d = 0.0;
            $scale = 1.0;
            foreach ($F as $i => $v) {
                $d = max($d, abs($v - $F2[2 * $i]));
                $scale = max($scale, abs($F2[2 * $i]));
            }
            if ($d <= $tol * $scale && self::agreesWithGauss($f, $T, $m, $F, $tol)) {
                $resolved = true;
                break;
            }
            $m *= 2;
            $F = $F2;
        }
        $info = ['m' => $m, 'resolved' => $resolved, 'max' => max($F), 'min' => min($F),
            'maxAbs' => max(abs(max($F)), abs(min($F)))];
        return [self::closure($f, $F, $T / $m, $m), $info];
    }

    /** Whether the Simpson table F on m panels agrees with a composite 3-point Gauss-Legendre table to tol. */
    private static function agreesWithGauss(\Closure $f, float $T, int $m, array $F, float $tol): bool
    {
        $h = $T / $m;
        $r = 0.5 * sqrt(0.6);
        $S = 0.0;
        $d = 0.0;
        $scale = 1.0;
        for ($i = 0; $i < $m; $i++) {
            $c = ($i + 0.5) * $h;
            $S += $h / 18.0 * (5.0 * $f($c - $r * $h) + 8.0 * $f($c) + 5.0 * $f($c + $r * $h));
            $d = max($d, abs($S - $F[$i + 1]));
            $scale = max($scale, abs($S));
        }
        return $d <= $tol * $scale;
    }

    /** Composite Simpson table F[i] = int_0^{i h} f, h = T/m. */
    private static function table(\Closure $f, float $T, int $m): array
    {
        $h = $T / $m;
        $F = [0.0];
        for ($i = 0; $i < $m; $i++) {
            $F[] = $F[$i] + self::simpson($f, $i * $h, ($i + 1) * $h);
        }
        return $F;
    }

    private static function closure(\Closure $f, array $F, float $h, int $m): \Closure
    {
        return static function (float $t) use ($f, $F, $h, $m): float {
            if ($t <= 0.0) {
                return 0.0;
            }
            $i = min($m - 1, (int)floor($t / $h));
            return $F[$i] + self::simpson($f, $i * $h, $t);
        };
    }

    /**
     * Discrepancy of the nonlinear phase quadrature of a split-step run. Each step applies to |q(t0)|^2 the phase
     * 2 h/6 [gamma(t0) + 4 gamma(tm) e^{2 Gm} + gamma(t1) e^{2 G1}], with Gm, G1 the Simpson half-step integrals of g
     * from t0 (SplitStepSolver::nonlinear); the exact phase is 2 int_{t0}^{t1} gamma(s) e^{2 int_{t0}^s g} ds.
     * Weighting step n by e^{2 G(t0)} makes the exact sum telescope to Gam(t_n) = int_0^{t_n} gamma e^{2G}. Returns
     * max_n |sum - Gam(t_n)| / max(1, max|Gam|).
     */
    public static function nonlinearQuadratureError(\Closure $gamma, \Closure $g, \Closure $Gam, \Closure $G, float $T, int $N): float
    {
        $dt = $T / $N;
        $S = 0.0;
        $d = 0.0;
        $scale = 1.0;
        for ($n = 1; $n <= $N; $n++) {
            $t0 = ($n - 1) * $dt;
            $tm = $t0 + 0.5 * $dt;
            $t1 = $n * $dt;
            $Gm = self::simpson($g, $t0, $tm);
            $G1 = $Gm + self::simpson($g, $tm, $t1);
            $step = $dt / 6.0 * ($gamma($t0) + 4.0 * $gamma($tm) * exp(2.0 * $Gm) + $gamma($t1) * exp(2.0 * $G1));
            $S += exp(2.0 * $G($t0)) * $step;
            $ref = $Gam($t1);
            $d = max($d, abs($S - $ref));
            $scale = max($scale, abs($ref));
        }
        return $d / $scale;
    }

    /**
     * Largest discrepancy, relative to max(1, max|F|), between the step quadrature of a split-step run
     * (composite Simpson on half steps, as in SplitStepSolver) and the reference antiderivative F, at the
     * step nodes t_n = n T/N.
     */
    public static function stepQuadratureError(\Closure $f, \Closure $F, float $T, int $N): float
    {
        $dt = $T / $N;
        $S = 0.0;
        $d = 0.0;
        $scale = 1.0;
        $f0 = $f(0.0);
        for ($n = 1; $n <= $N; $n++) {
            $t0 = ($n - 1) * $dt;
            $tm = $t0 + 0.5 * $dt;
            $t1 = $n * $dt;
            // Simpson on [t0, tm] and [tm, t1], as in SplitStepSolver, with the shared endpoint values reused
            $fm = $f($tm);
            $f1 = $f($t1);
            $S += $dt / 12.0 * ($f0 + 4.0 * $f(0.5 * ($t0 + $tm)) + 2.0 * $fm + 4.0 * $f(0.5 * ($tm + $t1)) + $f1);
            $f0 = $f1;
            $ref = $F($t1);
            $d = max($d, abs($S - $ref));
            $scale = max($scale, abs($ref));
        }
        return $d / $scale;
    }

    /**
     * Discrepancy of the potential quadrature of a split-step run at K sample points of the grid. Each step adds to
     * the phase the Simpson panel h/6 [V(x,t0) + 4 V(x,tm) + V(x,t1)] (SplitStepSolver::nonlinear); the exact phase is
     * int_{t0}^{t1} V(x,s) ds. $Vs(t) returns the K values of V at the sample points. The reference integrates
     * every step with M composite Simpson panels, M doubled from ceil(m0/N) until it agrees with 2M panels and with
     * three-point Gauss-Legendre on M panels to tol * scale (or until M N reaches mMax). Returns
     * max_{n,k} |sum_n - ref_n| / scale, scale = max(1, max |ref|). Like the other tests, an estimate: V is checked
     * only at the sample points.
     */
    public static function panelQuadratureErrorSampled(\Closure $Vs, float $T, int $N, float $tol = 1e-8,
                                                       int $m0 = 4096, int $mMax = 1048576): float
    {
        $dt = $T / $N;
        $M = max(2, (int)ceil($m0 / $N));
        $r = 0.5 * sqrt(0.6);
        while (true) {
            $K = null;
            $AM = $A2 = $AG = $S = [];
            $dRef = $dGL = $dSol = $dSol2 = 0.0;
            $scale = 1.0;
            $hq = $dt / (4 * $M);
            $hp = $dt / $M;
            $prev = $Vs(0.0);
            for ($n = 1; $n <= $N; $n++) {
                $t0 = ($n - 1) * $dt;
                $v = [$prev];
                for ($j = 1; $j <= 4 * $M; $j++) {
                    $v[$j] = $Vs($t0 + $j * $hq);
                }
                $prev = $v[4 * $M];
                if ($K === null) {
                    $K = count($prev);
                    $AM = $A2 = $AG = $S = array_fill(0, $K, 0.0);
                }
                for ($k = 0; $k < $K; $k++) {
                    $iM = $i2 = $iG = 0.0;
                    for ($p = 0; $p < $M; $p++) {
                        $iM += $v[4 * $p][$k] + 4.0 * $v[4 * $p + 2][$k] + $v[4 * $p + 4][$k];
                    }
                    for ($p = 0; $p < 2 * $M; $p++) {
                        $i2 += $v[2 * $p][$k] + 4.0 * $v[2 * $p + 1][$k] + $v[2 * $p + 2][$k];
                    }
                    $AM[$k] += $hp / 6.0 * $iM;
                    $A2[$k] += $hp / 12.0 * $i2;
                    $S[$k] += $dt / 6.0 * ($v[0][$k] + 4.0 * $v[2 * $M][$k] + $v[4 * $M][$k]);
                }
                for ($p = 0; $p < $M; $p++) {
                    $c = $t0 + ($p + 0.5) * $hp;
                    $lo = $Vs($c - $r * $hp);
                    $hi = $Vs($c + $r * $hp);
                    for ($k = 0; $k < $K; $k++) {
                        $AG[$k] += $hp / 18.0 * (5.0 * $lo[$k] + 8.0 * $v[4 * $p + 2][$k] + 5.0 * $hi[$k]);
                    }
                }
                for ($k = 0; $k < $K; $k++) {
                    $dRef = max($dRef, abs($AM[$k] - $A2[$k]));
                    $dGL = max($dGL, abs($AM[$k] - $AG[$k]));
                    $dSol = max($dSol, abs($S[$k] - $AM[$k]));
                    $dSol2 = max($dSol2, abs($S[$k] - $A2[$k]));
                    $scale = max($scale, abs($A2[$k]));
                }
            }
            if ($dRef <= $tol * $scale && $dGL <= $tol * $scale) {
                return $dSol / $scale;
            }
            if (2 * $M * $N > $mMax) {
                return $dSol2 / $scale;
            }
            $M *= 2;
        }
    }
}
