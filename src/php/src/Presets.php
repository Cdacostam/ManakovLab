<?php
declare(strict_types=1);

namespace ManakovLab;

/** Ready-made experiments shown in the preset menu. */
final class Presets
{
    public static function all(): array
    {
        $sol = static fn(float $eta, float $xi, float $x0, float $p1re, float $p1im, float $p2re, float $p2im): array =>
            compact('eta', 'xi', 'x0', 'p1re', 'p1im', 'p2re', 'p2im');

        return [
            [
                'id' => 'classic-exact',
                'group' => 'Variable-coefficient Manakov',
                'title' => 'Classical: nonautonomous soliton, exact + convergence',
                'description' => 'Periodic dispersion and oscillating gain. γ(t) is derived from the integrability condition, so the exact solution is known and the error and 2nd-order convergence can be measured.',
                'model' => 'classic',
                'params' => ['L' => 80, 'Nx' => 512, 'T' => 10, 'dt' => 0.005],
                'coeffs' => ['beta' => '1 + 0.5 cos(2t)', 'gamma' => '1', 'g' => '0.05 sin(t)', 'V' => '0'],
                'solitons' => [$sol(0.6, 0.15, 10, 1, 0, 0.6, 0.4)],
                'exact' => true, 'convergence' => true,
            ],
            [
                'id' => 'classic-collision',
                'group' => 'Variable-coefficient Manakov',
                'title' => 'Classical: shape-changing collision',
                'description' => 'Two vector solitons exchange energy between components. The coefficients satisfy γ = β·exp(−2∫g), so the system is gauge-equivalent to the integrable Manakov system.',
                'model' => 'classic',
                'params' => ['L' => 100, 'Nx' => 512, 'T' => 30, 'dt' => 0.01],
                'coeffs' => ['beta' => '1 + 0.5 cos(2t)', 'gamma' => '(1 + 0.5 cos(2t)) exp(-0.1 (1 - cos(t)))',
                    'g' => '0.05 sin(t)', 'V' => '0'],
                'solitons' => [$sol(0.5, -0.25, -15, 1, 0, 1, 0), $sol(0.6, 0.25, 15, 1, 0, 0, 0)],
                'exact' => false, 'convergence' => false,
            ],
            [
                'id' => 'classic-trap',
                'group' => 'Variable-coefficient Manakov',
                'title' => 'Classical: dispersion management in a harmonic trap',
                'description' => 'Non-integrable case: periodically modulated dispersion β(t) ∈ [0.2, 1.8] (no sign change) and a confining potential V = −ω²x².',
                'model' => 'classic',
                'params' => ['L' => 120, 'Nx' => 512, 'T' => 40, 'dt' => 0.01],
                'coeffs' => ['beta' => '1 + 0.8 sin(4t)', 'gamma' => '1', 'g' => '0', 'V' => '-0.01 x^2'],
                'solitons' => [$sol(0.5, 0, -12, 1, 0, 0, 0.5)],
                'exact' => false, 'convergence' => false,
            ],
            [
                'id' => 'space-exact',
                'group' => 'Space-fractional (Riesz)',
                'title' => 'Space-fractional: soliton α = 1.5, exact + convergence',
                'description' => 'Stationary fractional vector soliton computed by Petviashvili iteration, carried through variable β(t), g(t) by the gauge reduction. Exact error and convergence.',
                'model' => 'space',
                'params' => ['L' => 80, 'Nx' => 512, 'T' => 10, 'dt' => 0.005, 'alpha' => 1.5],
                'coeffs' => ['beta' => '1 + 0.5 cos(2t)', 'gamma' => '1', 'g' => '0.05 sin(t)', 'V' => '0'],
                'solitons' => [$sol(0.5, 0, 0, 1, 0, 0.6, 0.4)],
                'exact' => true, 'convergence' => true,
            ],
            [
                'id' => 'space-collision',
                'group' => 'Space-fractional (Riesz)',
                'title' => 'Space-fractional: inelastic collision, α = 1.5',
                'description' => 'Two kicked fractional solitons. Without integrability the collision radiates (inspect the log-scale final profile).',
                'model' => 'space',
                'params' => ['L' => 160, 'Nx' => 1024, 'T' => 25, 'dt' => 0.01, 'alpha' => 1.5],
                'coeffs' => ['beta' => '1', 'gamma' => '1', 'g' => '0', 'V' => '0'],
                'solitons' => [$sol(0.5, -0.5, -15, 1, 0, 1, 0), $sol(0.5, 0.5, 15, 1, 0, 0, 0)],
                'exact' => false, 'convergence' => false,
            ],
            [
                'id' => 'caputo-mms',
                'group' => 'Time-fractional (Caputo)',
                'title' => 'Caputo: manufactured solution ν = 0.6, exact + convergence',
                'description' => 'All coefficients vary (including V). The forcing is computed so q* = c(1+t^ν+it²)ψ(x) is exact. The graded mesh restores the optimal order 2−ν.',
                'model' => 'caputo',
                'params' => ['L' => 60, 'Nx' => 256, 'T' => 1, 'N' => 512, 'nu' => 0.6, 'r' => 'auto'],
                'coeffs' => ['beta' => '1 + 0.5 cos(2t)', 'gamma' => '1 + 0.3 sin(t)', 'g' => '0.1 cos(t)',
                    'V' => '0.05 cos(0.5x)(1 + t)'],
                'solitons' => [$sol(0.25, -0.25, 0, 1, 0, 0.5, -0.5)],
                'exact' => true, 'convergence' => true,
            ],
            [
                'id' => 'caputo-decay',
                'group' => 'Time-fractional (Caputo)',
                'title' => 'Caputo: soliton decay under memory, ν = 0.8',
                'description' => 'A Manakov soliton with periodic dispersion. For ν < 1 the memory term breaks power conservation and the soliton decays and spreads. L = 160 keeps the radiation of the breathing soliton from re-entering the periodic domain.',
                'model' => 'caputo',
                'params' => ['L' => 160, 'Nx' => 512, 'T' => 15, 'N' => 600, 'nu' => 0.8, 'r' => 'auto'],
                'coeffs' => ['beta' => '1 + 0.5 cos(2t)', 'gamma' => '1', 'g' => '0', 'V' => '0'],
                'solitons' => [$sol(0.6, 0.1, 5, 1, 0, 0, 0.7)],
                'exact' => false, 'convergence' => false,
            ],
            [
                'id' => 'caputo-collision',
                'group' => 'Time-fractional (Caputo)',
                'title' => 'Caputo: two-soliton interaction, ν = 0.9',
                'description' => 'Two vector solitons approaching each other under fractional time memory.',
                'model' => 'caputo',
                'params' => ['L' => 100, 'Nx' => 256, 'T' => 20, 'N' => 700, 'nu' => 0.9, 'r' => 'auto'],
                'coeffs' => ['beta' => '1', 'gamma' => '1', 'g' => '0', 'V' => '0'],
                'solitons' => [$sol(0.6, -0.4, -12, 1, 0, 1, 0), $sol(0.6, 0.4, 12, 1, 0, 0, 0)],
                'exact' => false, 'convergence' => false,
            ],
        ];
    }

    public static function byId(string $id): ?array
    {
        foreach (self::all() as $p) {
            if ($p['id'] === $id) {
                return $p;
            }
        }
        return null;
    }
}
