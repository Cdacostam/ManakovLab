<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

final class SimulationException extends \RuntimeException
{
    /** @param ?string $field input field the error refers to (shown by the interface), if any */
    public function __construct(string $message, public readonly ?string $field = null)
    {
        parent::__construct($message);
    }
}

/** Gamma function via the Lanczos approximation (g = 7, n = 9), ~15 significant digits. */
final class Gamma
{
    private const G = 7;
    private const C = [
        0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
        -176.61502916214059, 12.507343278686905, -0.13857109526572012,
        9.9843695780195716e-6, 1.5056327351493116e-7,
    ];

    public static function gamma(float $z): float
    {
        if ($z < 0.5) {
            return M_PI / (sin(M_PI * $z) * self::gamma(1.0 - $z));
        }
        $z -= 1.0;
        $x = self::C[0];
        for ($i = 1; $i < self::G + 2; $i++) {
            $x += self::C[$i] / ($z + $i);
        }
        $t = $z + self::G + 0.5;
        return sqrt(2 * M_PI) * $t ** ($z + 0.5) * exp(-$t) * $x;
    }
}
