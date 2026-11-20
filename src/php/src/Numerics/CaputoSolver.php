<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

/**
 * Time-fractional (Caputo) variable-coefficient Manakov system
 *
 *   i D_t^nu q + beta(t) q_xx + 2 gamma(t) (|q1|^2+|q2|^2) q + V(x,t) q = i g(t) q + f(x,t),   0 < nu <= 1
 *
 * L1 discretisation of the Caputo derivative on the graded mesh t_n = T (n/N)^r:
 *   D^nu u(t_n) ~ sum_{k=1}^{n} a_{n,k} (u^k - u^{k-1}),
 *   a_{n,k} = [(t_n - t_{k-1})^{1-nu} - (t_n - t_k)^{1-nu}] / (Gamma(2-nu) (t_k - t_{k-1})).
 * Fourier pseudo-spectral in space. Dispersion and gain implicit (diagonal in Fourier space);
 * nonlinear, potential and forcing terms implicit, solved by Anderson-accelerated fixed-point (Picard)
 * iteration. Full memory is kept.
 * With nu = 1 the L1 weights reduce to backward Euler on the graded mesh (first order, dissipative); this is
 * accepted for library use and limit checks, while Simulation runs nu = 1 with the split-step solver.
 */
final class CaputoSolver
{
    public readonly array $x;
    public readonly array $k;
    public readonly array $k2;
    public readonly float $dx;
    public readonly int $n;
    public int $maxPicard = 0;
    public int $sumPicard = 0;
    public int $stepsDone = 0;
    /** @var list<int> nonlinear iterations of each step of the last solve */
    public array $picardSteps = [];
    /** Depth of Anderson mixing (0 = plain Picard). */
    public int $anderson = 10;

    public function __construct(
        public readonly float $L,
        int $Nx,
        public readonly float $nu,
        public readonly \Closure $beta,
        public readonly \Closure $gamma,
        public readonly \Closure $g,
        public readonly ?\Closure $V,
        public readonly ?\Closure $forcing = null,
        private readonly ?float $deadline = null,
    ) {
        // Contract of the L1 scheme (also enforced by Simulation): 0 < nu <= 1; nu = 1 is backward Euler.
        if (!is_finite($nu) || $nu <= 0.0 || $nu > 1.0) {
            throw new SimulationException('CaputoSolver: the Caputo order must satisfy 0 < nu <= 1.', 'nu');
        }
        if (!is_finite($L) || $L <= 0.0) {
            throw new SimulationException('CaputoSolver: L must be positive.', 'L');
        }
        if ($Nx < 4 || ($Nx & ($Nx - 1)) !== 0) {
            throw new SimulationException('CaputoSolver: Nx must be a power of two, at least 4.', 'Nx');
        }
        $this->n = $Nx;
        $this->dx = $L / $Nx;
        $x = [];
        for ($i = 0; $i < $Nx; $i++) {
            $x[] = ($i - intdiv($Nx, 2)) * $this->dx;
        }
        $this->x = $x;
        $this->k = FFT::wavenumbers($Nx, $L);
        $this->k2 = array_map(static fn($v) => $v * $v, $this->k);
    }

    public static function mesh(float $T, int $N, float $r): array
    {
        $t = [];
        for ($i = 0; $i <= $N; $i++) {
            $t[] = $T * ($i / $N) ** $r;
        }
        return $t;
    }

    /** Nonlinear + potential + forcing terms, divided by i:  2i gamma |q|^2 q + i V q - i f. */
    private function rhs(array $Q, float $t): array
    {
        [$r1, $i1, $r2, $i2] = $Q;
        $gam2 = 2.0 * ($this->gamma)($t);
        $V = $this->V !== null ? ($this->V)($this->x, $t) : null;
        $o1 = $o2 = $o3 = $o4 = [];
        foreach ($r1 as $i => $a) {
            $w = $gam2 * ($a * $a + $i1[$i] * $i1[$i] + $r2[$i] * $r2[$i] + $i2[$i] * $i2[$i]);
            if ($V !== null) {
                $w += $V[$i];
            }
            $o1[$i] = -$w * $i1[$i];
            $o2[$i] = $w * $a;
            $o3[$i] = -$w * $i2[$i];
            $o4[$i] = $w * $r2[$i];
        }
        if ($this->forcing !== null) {
            [$f1, $g1, $f2, $g2] = ($this->forcing)($t);
            foreach ($o1 as $i => $_) {
                $o1[$i] += $g1[$i];
                $o2[$i] -= $f1[$i];
                $o3[$i] += $g2[$i];
                $o4[$i] -= $f2[$i];
            }
        }
        return [$o1, $o2, $o3, $o4];
    }

    private static function fft4(array $Q, bool $inverse = false): array
    {
        [$a, $b, $c, $d] = $Q;
        FFT::transform($a, $b, $inverse);
        FFT::transform($c, $d, $inverse);
        return [$a, $b, $c, $d];
    }

    /**
     * One application of the fixed-point map  q -> F^{-1}[(known + F[rhs(q)]) / (a_nn - lambda)].
     * @return array{0: array, 1: array} the new field in physical and in Fourier space
     */
    private function fixedPointMap(array $Qn, float $tn, array $known, array $invR, array $invI): array
    {
        $Nh = self::fft4($this->rhs($Qn, $tn));
        $QhNew = [];
        for ($c = 0; $c < 4; $c += 2) {
            $re = []; $im = [];
            $kr = $known[$c]; $ki = $known[$c + 1]; $nr = $Nh[$c]; $ni = $Nh[$c + 1];
            foreach ($kr as $i => $v) {
                $ar = $v + $nr[$i];
                $ai = $ki[$i] + $ni[$i];
                $re[$i] = $ar * $invR[$i] - $ai * $invI[$i];
                $im[$i] = $ar * $invI[$i] + $ai * $invR[$i];
            }
            $QhNew[$c] = $re;
            $QhNew[$c + 1] = $im;
        }
        return [self::fft4($QhNew, true), $QhNew];
    }

    /**
     * Anderson (type II) update from the iterate history X and residual history F (flat real vectors):
     * x_new = x_k + f_k - (dX + dF) gamma,  gamma = argmin |f_k - dF gamma|  (regularised normal equations).
     */
    private static function andersonStep(array $X, array $F, int $Nx): array
    {
        $k = count($X) - 1;
        $len = count($X[0]);
        $dF = $dX = [];
        for ($j = 0; $j < $k; $j++) {
            $a = $F[$j + 1]; $b = $F[$j]; $u = $X[$j + 1]; $w = $X[$j];
            $df = $dx = [];
            for ($i = 0; $i < $len; $i++) {
                $df[$i] = $a[$i] - $b[$i];
                $dx[$i] = $u[$i] - $w[$i];
            }
            $dF[$j] = $df;
            $dX[$j] = $dx;
        }
        $fk = $F[$k];
        $A = [];
        $rhs = [];
        $trace = 0.0;
        for ($p = 0; $p < $k; $p++) {
            for ($q = $p; $q < $k; $q++) {
                $sum = 0.0;
                $dp = $dF[$p]; $dq = $dF[$q];
                for ($i = 0; $i < $len; $i++) {
                    $sum += $dp[$i] * $dq[$i];
                }
                $A[$p][$q] = $A[$q][$p] = $sum;
            }
            $trace += $A[$p][$p];
            $sum = 0.0;
            $dp = $dF[$p];
            for ($i = 0; $i < $len; $i++) {
                $sum += $dp[$i] * $fk[$i];
            }
            $rhs[$p] = $sum;
        }
        $reg = 1e-12 * ($trace > 0.0 ? $trace / $k : 1.0);
        for ($p = 0; $p < $k; $p++) {
            $A[$p][$p] += $reg;
        }
        $g = self::solveSmall($A, $rhs);
        $xk = $X[$k];
        $out = [];
        for ($i = 0; $i < $len; $i++) {
            $v = $xk[$i] + $fk[$i];
            for ($j = 0; $j < $k; $j++) {
                $v -= ($dX[$j][$i] + $dF[$j][$i]) * $g[$j];
            }
            $out[$i] = $v;
        }
        return [array_slice($out, 0, $Nx), array_slice($out, $Nx, $Nx), array_slice($out, 2 * $Nx, $Nx), array_slice($out, 3 * $Nx, $Nx)];
    }

    /** Gaussian elimination with partial pivoting for the small (at most 10x10, the Anderson depth) least-squares system. */
    private static function solveSmall(array $A, array $b): array
    {
        $n = count($b);
        for ($c = 0; $c < $n; $c++) {
            $piv = $c;
            for ($r = $c + 1; $r < $n; $r++) {
                if (abs($A[$r][$c]) > abs($A[$piv][$c])) {
                    $piv = $r;
                }
            }
            if ($piv !== $c) {
                [$A[$c], $A[$piv]] = [$A[$piv], $A[$c]];
                [$b[$c], $b[$piv]] = [$b[$piv], $b[$c]];
            }
            $d = $A[$c][$c];
            if ($d == 0.0) {
                continue;
            }
            for ($r = $c + 1; $r < $n; $r++) {
                $f = $A[$r][$c] / $d;
                for ($j = $c; $j < $n; $j++) {
                    $A[$r][$j] -= $f * $A[$c][$j];
                }
                $b[$r] -= $f * $b[$c];
            }
        }
        $x = array_fill(0, $n, 0.0);
        for ($r = $n - 1; $r >= 0; $r--) {
            $s = $b[$r];
            for ($j = $r + 1; $j < $n; $j++) {
                $s -= $A[$r][$j] * $x[$j];
            }
            $x[$r] = $A[$r][$r] != 0.0 ? $s / $A[$r][$r] : 0.0;
        }
        return $x;
    }

    /**
     * @param callable(float, array): void $onSave
     * @return array final field [re1, im1, re2, im2]
     */
    public function solve(array $Q, float $T, int $N, float $r, int $nSave, callable $onSave,
                          float $tol = 1e-12, int $maxit = 100): array
    {
        if (!is_finite($T) || $T <= 0.0 || $N < 1 || !is_finite($r) || $r < 1.0) {
            throw new SimulationException('CaputoSolver: requires T > 0, N >= 1 and r >= 1.', $N < 1 ? 'N' : ($r < 1.0 ? 'r' : 'T'));
        }
        $nu = $this->nu;
        $t = self::mesh($T, $N, $r);
        $G2 = Gamma::gamma(2.0 - $nu);
        $Nx = $this->n;
        $e = 1.0 - $nu;

        // save at (approximately) uniformly spaced times
        $save = [];
        $m = min($nSave, $N);
        for ($s = 0, $j = 1; $s <= $m; $s++) {
            $target = $T * $s / $m;
            while ($j < $N && $t[$j] < $target) {
                $j++;
            }
            $save[max(1, $j)] = true;
        }
        $save[$N] = true;

        $onSave(0.0, $Q);
        $Qh = self::fft4($Q);
        $dU = [];                                  // dU[k] = flat (4*Nx) Fourier increment
        $this->maxPicard = 0;
        $this->sumPicard = 0;
        $this->stepsDone = 0;
        $this->picardSteps = [];

        for ($n = 1; $n <= $N; $n++) {
            $tn = $t[$n];
            // L1 weights
            $a = [];
            // (t_n - t_{k-1})^e - (t_n - t_k)^e = w^e expm1(e log1p(tau_k / w)),  w = t_n - t_k:
            // the direct difference loses all digits when tau_k << w (first steps of a strongly graded mesh)
            for ($kk = 1; $kk <= $n; $kk++) {
                $tau = $t[$kk] - $t[$kk - 1];
                $w = $tn - $t[$kk];
                $d = $w > 0 ? $w ** $e * expm1($e * log1p($tau / $w)) : $tau ** $e;
                $a[$kk] = $d / ($G2 * $tau);
            }
            $ann = $a[$n];

            // history sum H = sum_{k<n} a_k dU_k  (flat)
            $H = array_fill(0, 4 * $Nx, 0.0);
            for ($kk = 1; $kk < $n; $kk++) {
                $w = $a[$kk];
                foreach ($dU[$kk] as $i => $v) {
                    $H[$i] += $w * $v;
                }
            }
            $known = [];
            for ($j = 0; $j < 4; $j++) {
                $off = $j * $Nx;
                $col = $Qh[$j];
                $kn = [];
                foreach ($col as $i => $v) {
                    $kn[$i] = $ann * $v - $H[$off + $i];
                }
                $known[$j] = $kn;
            }
            // 1 / (ann - lambda),  lambda = -i beta k^2 + g
            $bt = ($this->beta)($tn);
            $gt = ($this->g)($tn);
            $invR = $invI = [];
            foreach ($this->k2 as $i => $k2) {
                $dr = $ann - $gt;
                $di = $bt * $k2;
                $D = $dr * $dr + $di * $di;
                $invR[$i] = $dr / $D;
                $invI[$i] = -$di / $D;
            }

            // Nonlinear solve: fixed-point (Picard) map accelerated by Anderson mixing.
            // The plain Picard map is a contraction only when Lip(N) < |a_nn - g|; Anderson
            // acceleration converges far beyond that regime (small nu, larger amplitudes)
            // and returns the same fixed point. Stopping test: max|G(q) - q| < tol * max(1, max|G(q)|).
            $Qn = $Q;
            $QhNew = null;
            $converged = false;
            $diff = INF;
            $hX = [];
            $hF = [];
            for ($it = 1; $it <= $maxit; $it++) {
                [$GQ, $GQh] = $this->fixedPointMap($Qn, $tn, $known, $invR, $invI);
                $diff = 0.0;
                $scale = 1.0;
                $xv = [];
                $fv = [];
                for ($c = 0; $c < 4; $c++) {
                    $p = $Qn[$c];
                    foreach ($GQ[$c] as $i => $v) {
                        $d = $v - $p[$i];
                        $xv[] = $p[$i];
                        $fv[] = $d;
                        $ad = abs($d);
                        if ($ad > $diff) {
                            $diff = $ad;
                        }
                        $av = abs($v);
                        if ($av > $scale) {
                            $scale = $av;
                        }
                    }
                }
                if (!is_finite($diff)) {
                    break;
                }
                if ($diff < $tol * $scale) {
                    $Qn = $GQ;
                    $QhNew = $GQh;
                    $converged = true;
                    break;
                }
                if ($this->anderson < 1) {
                    $Qn = $GQ;
                    continue;
                }
                $hX[] = $xv;
                $hF[] = $fv;
                if (count($hX) > $this->anderson + 1) {
                    array_shift($hX);
                    array_shift($hF);
                }
                $Qn = count($hX) === 1 ? $GQ : self::andersonStep($hX, $hF, $Nx);
            }
            if (!$converged) {
                throw new SimulationException(sprintf(
                    'The implicit nonlinear solve did not converge at step %d (t = %.4g, last change %.2e). '
                    . 'The nonlinearity is too strong for this step: increase N, reduce the amplitude, or increase nu.',
                    $n, $tn, $diff), 'N');
            }
            $this->maxPicard = max($this->maxPicard, $it);
            $this->sumPicard += $it;
            $this->picardSteps[] = $it;

            $flat = [];
            for ($c = 0; $c < 4; $c++) {
                $old = $Qh[$c];
                foreach ($QhNew[$c] as $i => $v) {
                    $flat[] = $v - $old[$i];
                }
            }
            $dU[$n] = $flat;
            $this->stepsDone = $n;
            $Qh = $QhNew;
            $Q = $Qn;
            if (isset($save[$n])) {
                $onSave($tn, $Q);
            }
            if (($n & 7) === 0 && $this->deadline !== null && microtime(true) > $this->deadline) {
                throw new SimulationException('The simulation exceeded the time limit. Reduce N or Nx.');
            }
        }
        return $Q;
    }
}
