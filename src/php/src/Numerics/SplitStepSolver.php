<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

/**
 * Strang split-step Fourier solver for
 *
 *   i q_t - beta(t) (-d_xx)^{alpha/2} q + 2 gamma(t) (|q1|^2 + |q2|^2) q + V(x,t) q = i g(t) q
 *
 * (alpha = 2 is the classical variable-coefficient Manakov system, since -(-d_xx) = d_xx).
 *
 *  Linear part   : exact in Fourier space, q^ <- q^ exp(-i |k|^alpha int beta dt)
 *  Nonlinear part: exact pointwise; |q|^2 grows as exp(2 int g), so
 *                  q <- q exp(int g) exp(i int [2 gamma |q|^2 + V] dt)
 *  Time integrals: Simpson's rule (keeps 2nd order with time-dependent coefficients).
 *  Consecutive linear half-steps are fused between output times (halves the FFTs: two instead of four per component and step).
 *
 * Field layout: [re1, im1, re2, im2].
 */
final class SplitStepSolver
{
    public readonly array $x;
    public readonly array $k;
    public readonly array $Dk;
    public readonly float $dx;
    public readonly int $n;

    public function __construct(
        public readonly float $L,
        int $Nx,
        public readonly float $alpha,
        public readonly \Closure $beta,
        public readonly \Closure $gamma,
        public readonly \Closure $g,
        public readonly ?\Closure $V,
        private readonly ?float $deadline = null,
    ) {
        if (!is_finite($L) || $L <= 0.0) {
            throw new SimulationException('SplitStepSolver: L must be positive.', 'L');
        }
        if ($Nx < 4 || ($Nx & ($Nx - 1)) !== 0) {
            throw new SimulationException('SplitStepSolver: Nx must be a power of two, at least 4.', 'Nx');
        }
        if (!is_finite($alpha) || $alpha <= 0.0 || $alpha > 2.0) {
            throw new SimulationException('SplitStepSolver: the order must satisfy 0 < alpha <= 2.', 'alpha');
        }
        $this->n = $Nx;
        $this->dx = $L / $Nx;
        $x = [];
        for ($i = 0; $i < $Nx; $i++) {
            $x[] = ($i - intdiv($Nx, 2)) * $this->dx;
        }
        $this->x = $x;
        $this->k = FFT::wavenumbers($Nx, $L);
        $Dk = [];
        foreach ($this->k as $kv) {
            $Dk[] = abs($kv) ** $alpha;
        }
        $this->Dk = $Dk;
    }

    /** Apply exp(-i Dk B) to both components. */
    private function linear(array &$Q, float $B): void
    {
        $c = [];
        $s = [];
        foreach ($this->Dk as $i => $d) {
            $c[$i] = cos($d * $B);
            $s[$i] = -sin($d * $B);
        }
        for ($j = 0; $j < 4; $j += 2) {
            $re = $Q[$j];
            $im = $Q[$j + 1];
            FFT::transform($re, $im);
            foreach ($c as $i => $ci) {
                $r = $re[$i];
                $re[$i] = $r * $ci - $im[$i] * $s[$i];
                $im[$i] = $r * $s[$i] + $im[$i] * $ci;
            }
            FFT::transform($re, $im, true);
            $Q[$j] = $re;
            $Q[$j + 1] = $im;
        }
    }

    private function nonlinear(array &$Q, float $t0, float $t1): void
    {
        $g = $this->g;
        $gam = $this->gamma;
        $tm = 0.5 * ($t0 + $t1);
        $h = $t1 - $t0;
        $Gm = Integral::simpson($g, $t0, $tm);
        $G1 = $Gm + Integral::simpson($g, $tm, $t1);
        $cNl = $h / 6.0 * 2.0 * ($gam($t0) + 4.0 * $gam($tm) * exp(2.0 * $Gm) + $gam($t1) * exp(2.0 * $G1));
        $amp = exp($G1);
        [$r1, $i1, $r2, $i2] = $Q;
        $Vp = null;
        if ($this->V !== null) {
            $V = $this->V;
            $v0 = $V($this->x, $t0);
            $vm = $V($this->x, $tm);
            $v1 = $V($this->x, $t1);
            $Vp = [];
            foreach ($v0 as $i => $a) {
                $Vp[$i] = $h / 6.0 * ($a + 4.0 * $vm[$i] + $v1[$i]);
            }
        }
        foreach ($r1 as $i => $a) {
            $I = $a * $a + $i1[$i] * $i1[$i] + $r2[$i] * $r2[$i] + $i2[$i] * $i2[$i];
            $ph = $cNl * $I + ($Vp === null ? 0.0 : $Vp[$i]);
            $c = $amp * cos($ph);
            $s = $amp * sin($ph);
            $b = $i1[$i];
            $r1[$i] = $a * $c - $b * $s;
            $i1[$i] = $a * $s + $b * $c;
            $a2 = $r2[$i];
            $b2 = $i2[$i];
            $r2[$i] = $a2 * $c - $b2 * $s;
            $i2[$i] = $a2 * $s + $b2 * $c;
        }
        $Q = [$r1, $i1, $r2, $i2];
    }

    /**
     * Integrate from 0 to T with time step dt.
     * $onSave(float $t, array $Q) is called at t = 0, at ~nSave uniformly spaced times and at T.
     * @return array final field
     */
    public function solve(array $Q, float $T, float $dt, int $nSave, callable $onSave): array
    {
        if (!is_finite($T) || $T <= 0.0 || !is_finite($dt) || $dt <= 0.0) {
            throw new SimulationException('SplitStepSolver: T and dt must be positive.', 'dt');
        }
        $beta = $this->beta;
        $N = max(1, (int)round($T / $dt));
        $dt = $T / $N;
        $every = max(1, intdiv($N, max(1, $nSave)));
        $onSave(0.0, $Q);

        $this->linear($Q, Integral::simpson($beta, 0.0, 0.5 * $dt));
        for ($n = 1; $n <= $N; $n++) {
            $t0 = ($n - 1) * $dt;
            $t1 = $n * $dt;
            $this->nonlinear($Q, $t0, $t1);
            $Bsecond = Integral::simpson($beta, $t0 + 0.5 * $dt, $t1);
            if ($n % $every === 0 || $n === $N) {
                $this->linear($Q, $Bsecond);
                $onSave($t1, $Q);
                if ($n < $N) {
                    $this->linear($Q, Integral::simpson($beta, $t1, $t1 + 0.5 * $dt));
                }
            } else {
                $this->linear($Q, $Bsecond + Integral::simpson($beta, $t1, $t1 + 0.5 * $dt));
            }
            if (($n & 31) === 0 && $this->deadline !== null && microtime(true) > $this->deadline) {
                throw new SimulationException('The simulation exceeded the time limit. Reduce Nx, the number of time steps, or T.');
            }
        }
        return $Q;
    }
}
