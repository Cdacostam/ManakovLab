<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

/**
 * Iterative radix-2 Cooley–Tukey FFT on split real/imaginary arrays.
 *
 * Conventions match NumPy: forward  X_k = sum_j x_j e^{-2 pi i jk/n},
 *                          inverse  x_j = (1/n) sum_k X_k e^{+2 pi i jk/n}.
 * Plans (bit-reversal permutation and twiddle tables) are cached per size.
 */
final class FFT
{
    /** @var array<int, array{0: list<int>, 1: list<float>, 2: list<float>}> */
    private static array $plans = [];

    /** @return array{0: list<int>, 1: list<float>, 2: list<float>} */
    private static function plan(int $n): array
    {
        if (isset(self::$plans[$n])) {
            return self::$plans[$n];
        }
        if ($n < 2 || ($n & ($n - 1)) !== 0) {
            throw new \InvalidArgumentException("FFT size must be a power of two, got $n");
        }
        $bits = (int)round(log($n, 2));
        $rev = [];
        for ($i = 0; $i < $n; $i++) {
            $r = 0;
            for ($b = 0, $v = $i; $b < $bits; $b++, $v >>= 1) {
                $r = ($r << 1) | ($v & 1);
            }
            $rev[] = $r;
        }
        $cos = [];
        $sin = [];
        $half = $n >> 1;
        for ($k = 0; $k < $half; $k++) {
            $a = 2.0 * M_PI * $k / $n;
            $cos[] = cos($a);
            $sin[] = sin($a);
        }
        return self::$plans[$n] = [$rev, $cos, $sin];
    }

    /**
     * In-place transform.
     * @param list<float> $re
     * @param list<float> $im
     */
    public static function transform(array &$re, array &$im, bool $inverse = false): void
    {
        $n = count($re);
        [$rev, $cosT, $sinT] = self::plan($n);

        for ($i = 0; $i < $n; $i++) {
            $j = $rev[$i];
            if ($j > $i) {
                $t = $re[$i]; $re[$i] = $re[$j]; $re[$j] = $t;
                $t = $im[$i]; $im[$i] = $im[$j]; $im[$j] = $t;
            }
        }

        $sign = $inverse ? 1.0 : -1.0;
        for ($size = 2; $size <= $n; $size <<= 1) {
            $half = $size >> 1;
            $step = intdiv($n, $size);
            for ($start = 0; $start < $n; $start += $size) {
                $k = 0;
                $end = $start + $half;
                for ($a = $start; $a < $end; $a++) {
                    $wr = $cosT[$k];
                    $wi = $sign * $sinT[$k];
                    $k += $step;
                    $b = $a + $half;
                    $xr = $re[$b];
                    $xi = $im[$b];
                    $tr = $wr * $xr - $wi * $xi;
                    $ti = $wr * $xi + $wi * $xr;
                    $ur = $re[$a];
                    $ui = $im[$a];
                    $re[$b] = $ur - $tr;
                    $im[$b] = $ui - $ti;
                    $re[$a] = $ur + $tr;
                    $im[$a] = $ui + $ti;
                }
            }
        }

        if ($inverse) {
            $s = 1.0 / $n;
            for ($i = 0; $i < $n; $i++) {
                $re[$i] *= $s;
                $im[$i] *= $s;
            }
        }
    }

    /** Angular wavenumbers in NumPy fftfreq order: 2 pi * fftfreq(n, L/n). */
    public static function wavenumbers(int $n, float $L): array
    {
        $k = [];
        $f = 2.0 * M_PI / $L;
        for ($m = 0; $m < $n; $m++) {
            $k[] = $f * ($m < intdiv($n, 2) ? $m : $m - $n);
        }
        return $k;
    }
}
