<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

/**
 * Soliton profiles. A two-component field is stored as [re1, im1, re2, im2].
 * A polarization is [c1re, c1im, c2re, c2im] and is normalised to |c| = 1.
 */
final class Solitons
{
    public static function normalise(array $c): array
    {
        $n = sqrt($c[0] ** 2 + $c[1] ** 2 + $c[2] ** 2 + $c[3] ** 2);
        if ($n == 0.0) {
            throw new \InvalidArgumentException('polarization vector cannot be zero');
        }
        return [$c[0] / $n, $c[1] / $n, $c[2] / $n, $c[3] / $n];
    }

    public static function zeros(int $n): array
    {
        $z = array_fill(0, $n, 0.0);
        return [$z, $z, $z, $z];
    }

    /** Q += c * envelope(x) * exp(i phase(x)), with complex polarization c. */
    public static function addScaled(array &$Q, array $c, array $env, array $phase): void
    {
        foreach ($env as $i => $a) {
            $cr = $a * cos($phase[$i]);
            $ci = $a * sin($phase[$i]);
            $Q[0][$i] += $c[0] * $cr - $c[1] * $ci;
            $Q[1][$i] += $c[0] * $ci + $c[1] * $cr;
            $Q[2][$i] += $c[2] * $cr - $c[3] * $ci;
            $Q[3][$i] += $c[2] * $ci + $c[3] * $cr;
        }
    }

    /**
     * Bright Manakov soliton of i p_tau + p_xx + 2|p|^2 p = 0:
     *   p = c 2 eta sech(2 eta (x - x0 + 4 xi tau)) exp(-2 i xi x - 4 i (xi^2 - eta^2) tau)
     */
    public static function manakov(array $x, float $tau, float $eta, float $xi, array $c, float $x0, float $scale = 1.0): array
    {
        $Q = self::zeros(count($x));
        self::addManakov($Q, $x, $tau, $eta, $xi, $c, $x0, $scale);
        return $Q;
    }

    public static function addManakov(array &$Q, array $x, float $tau, float $eta, float $xi, array $c, float $x0, float $scale = 1.0): void
    {
        $env = [];
        $ph = [];
        $w = 2.0 * $eta;
        $p0 = -4.0 * ($xi * $xi - $eta * $eta) * $tau;
        foreach ($x as $i => $xv) {
            $env[$i] = $scale * $w / cosh($w * ($xv - $x0 + 4.0 * $xi * $tau));
            $ph[$i] = -2.0 * $xi * $xv + $p0;
        }
        self::addScaled($Q, $c, $env, $ph);
    }

    /**
     * The Manakov soliton periodized over a box of length L: sum_m p(x + mL), over the images whose tail
     * reaches the box above 1e-17 of the peak. It solves the periodic linear problem exactly and the
     * nonlinear one up to the interaction with the images, and stays valid when the soliton crosses x = +-L/2.
     */
    public static function manakovPeriodic(array $x, float $L, float $tau, float $eta, float $xi, array $c, float $x0,
                                           float $scale = 1.0): array
    {
        $R = log(2.0e17) / (2.0 * $eta);                 // beyond R from its centre a soliton is below 1e-17 of its peak
        $cc = $x0 - 4.0 * $xi * $tau;                    // centre of p(x); p(x + mL) is centred at cc - mL
        $lo = (int)ceil(($cc - 0.5 * $L - $R) / $L);
        $hi = (int)floor(($cc + 0.5 * $L + $R) / $L);
        $Q = self::zeros(count($x));
        self::addManakov($Q, $x, $tau, $eta, $xi, $c, $x0, $scale);
        for ($m = $lo; $m <= $hi; $m++) {
            if ($m === 0) {
                continue;
            }
            // an image contributes only within R of its centre, i.e. near one end of the box
            $idx = [];
            $xs = [];
            foreach ($x as $i => $v) {
                if (abs($v + $m * $L - $cc) < $R) {
                    $idx[] = $i;
                    $xs[] = $v + $m * $L;
                }
            }
            if (!$idx) {
                continue;
            }
            $P = self::zeros(count($xs));
            self::addManakov($P, $xs, $tau, $eta, $xi, $c, $x0, $scale);
            foreach ($idx as $j => $i) {
                $Q[0][$i] += $P[0][$j]; $Q[1][$i] += $P[1][$j];
                $Q[2][$i] += $P[2][$j]; $Q[3][$i] += $P[3][$j];
            }
        }
        return $Q;
    }

    /**
     * Real, even solution of  (mu + (-d_xx)^{alpha/2}) phi = 2 phi^3  by Petviashvili iteration.
     * @return array{0: list<float>, 1: float, 2: int} profile, residual, iterations
     */
    public static function fractional(array $x, array $k, float $alpha, float $mu, float $tol = 1e-13, int $maxit = 3000): array
    {
        $n = count($x);
        $Lk = [];
        foreach ($k as $i => $kv) {
            $Lk[$i] = $mu + abs($kv) ** $alpha;
        }
        $sq = sqrt($mu);
        $phi = [];
        foreach ($x as $i => $xv) {
            $phi[$i] = $sq / cosh($sq * $xv);
        }
        $zero = array_fill(0, $n, 0.0);
        $it = 0;
        for (; $it < $maxit; $it++) {
            $pr = $phi; $pi = $zero;
            FFT::transform($pr, $pi);
            $nr = []; foreach ($phi as $i => $v) { $nr[$i] = 2.0 * $v * $v * $v; }
            $ni = $zero;
            FFT::transform($nr, $ni);
            $num = 0.0; $den = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $num += $Lk[$i] * ($pr[$i] * $pr[$i] + $pi[$i] * $pi[$i]);
                $den += $pr[$i] * $nr[$i] + $pi[$i] * $ni[$i];
            }
            $M = $num / $den;
            $f = $M ** 1.5;
            for ($i = 0; $i < $n; $i++) {
                $nr[$i] = $f * $nr[$i] / $Lk[$i];
                $ni[$i] = $f * $ni[$i] / $Lk[$i];
            }
            FFT::transform($nr, $ni, true);
            $diff = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $d = abs($nr[$i] - $phi[$i]);
                if ($d > $diff) {
                    $diff = $d;
                }
            }
            $phi = $nr;
            if ($diff < $tol) {
                break;
            }
        }
        // residual
        $pr = $phi; $pi = $zero;
        FFT::transform($pr, $pi);
        for ($i = 0; $i < $n; $i++) { $pr[$i] *= $Lk[$i]; $pi[$i] *= $Lk[$i]; }
        FFT::transform($pr, $pi, true);
        $res = 0.0;
        foreach ($phi as $i => $v) {
            $res = max($res, abs($pr[$i] - 2.0 * $v * $v * $v));
        }
        return [$phi, $res, $it + 1];
    }

    /**
     * Relative residual max|(mu + |D|^alpha) phi - 2 phi^3| / max|2 phi^3| of a real profile on the grid.
     */
    public static function relativeResidual(array $phi, array $k, float $alpha, float $mu): float
    {
        $pr = $phi;
        $pi = array_fill(0, count($phi), 0.0);
        FFT::transform($pr, $pi);
        foreach ($k as $i => $kv) {
            $l = $mu + abs($kv) ** $alpha;
            $pr[$i] *= $l;
            $pi[$i] *= $l;
        }
        FFT::transform($pr, $pi, true);
        $res = 0.0;
        $scale = 0.0;
        foreach ($phi as $i => $v) {
            $c = 2.0 * $v * $v * $v;
            $res = max($res, abs($pr[$i] - $c));
            $scale = max($scale, abs($c));
        }
        return $scale > 0.0 ? $res / $scale : INF;
    }

    /** Periodic shift f(x - x0) by a Fourier phase (exact for band-limited data). */
    public static function shift(array $f, array $k, float $x0): array
    {
        if ($x0 == 0.0) {
            return $f;
        }
        $re = $f;
        $im = array_fill(0, count($f), 0.0);
        FFT::transform($re, $im);
        foreach ($k as $i => $kv) {
            $c = cos($kv * $x0);
            $s = -sin($kv * $x0);
            $r = $re[$i] * $c - $im[$i] * $s;
            $im[$i] = $re[$i] * $s + $im[$i] * $c;
            $re[$i] = $r;
        }
        FFT::transform($re, $im, true);
        return $re;
    }
}
