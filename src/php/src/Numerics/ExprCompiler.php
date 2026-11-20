<?php
declare(strict_types=1);

namespace ManakovLab\Numerics;

/**
 * Compiles a validated AST from ExprParser into fast PHP closures.
 *
 * The compiled closure takes an environment ['t' => float, 'x' => float|list<float>]
 * and returns a float, or a list<float> when an array variable is involved
 * (vectorised: loops run in native PHP, one closure call per AST node, not per point).
 * Sub-trees without variables are folded to constants at compile time.
 * IEEE semantics: division uses fdiv(), so 1/0 = INF and 0/0 = NAN, never an exception.
 */
final class ExprCompiler
{
    /** Scalar time function f(t). */
    public static function timeFunction(array $ast): \Closure
    {
        $f = self::compile($ast);
        return static fn(float $t): float => (float)$f(['t' => $t]);
    }

    /**
     * Field function V(x, t) -> list<float>, or null when V is identically zero.
     * Time-independent fields are evaluated once and cached.
     */
    public static function fieldFunction(array $ast, array $x): ?\Closure
    {
        $f = self::compile($ast);
        $n = count($x);
        $vars = self::variables($ast);
        $expand = static function ($v) use ($n): array {
            return is_array($v) ? $v : array_fill(0, $n, (float)$v);
        };
        if (!in_array('x', $vars, true) && !in_array('t', $vars, true)) {
            $c = (float)$f([]);
            if ($c == 0.0) {
                return null;
            }
            $arr = array_fill(0, $n, $c);
            return static fn(array $xx, float $t): array => $arr;
        }
        if (!in_array('t', $vars, true)) {
            $arr = $expand($f(['x' => $x, 't' => 0.0]));
            return static fn(array $xx, float $t): array => $arr;
        }
        return static fn(array $xx, float $t): array => $expand($f(['x' => $xx, 't' => $t]));
    }

    public static function variables(array $n): array
    {
        $out = [];
        $walk = static function (array $n) use (&$walk, &$out): void {
            switch ($n['type']) {
                case 'var': $out[$n['name']] = true; break;
                case 'call': $walk($n['arg']); break;
                case 'neg': $walk($n['e']); break;
                case 'bin': $walk($n['l']); $walk($n['r']); break;
            }
        };
        $walk($n);
        return array_keys($out);
    }

    public static function compile(array $n): \Closure
    {
        if (!self::variables($n)) {
            $raw = self::compileNode($n);
            $v = $raw([]);
            return static fn(array $e) => $v;
        }
        return self::compileNode($n);
    }

    private static function compileNode(array $n): \Closure
    {
        switch ($n['type']) {
            case 'num':
                $v = (float)$n['v'];
                return static fn(array $e): float => $v;
            case 'const':
                $v = $n['name'] === 'pi' ? M_PI : M_E;
                return static fn(array $e): float => $v;
            case 'var':
                $name = $n['name'];
                return static fn(array $e) => $e[$name];
            case 'neg':
                $f = self::compile($n['e']);
                return static function (array $e) use ($f) {
                    $a = $f($e);
                    if (is_array($a)) {
                        foreach ($a as $i => $v) {
                            $a[$i] = -$v;
                        }
                        return $a;
                    }
                    return -$a;
                };
            case 'call':
                $f = self::compile($n['arg']);
                $fn = self::scalarFunction($n['fn']);
                return static function (array $e) use ($f, $fn) {
                    $a = $f($e);
                    return is_array($a) ? array_map($fn, $a) : $fn((float)$a);
                };
            case 'bin':
                $l = self::compile($n['l']);
                $r = self::compile($n['r']);
                $op = $n['op'];
                return static fn(array $e) => self::binop($op, $l($e), $r($e));
        }
        throw new \LogicException('unknown AST node');
    }

    /** @return callable(float): float */
    private static function scalarFunction(string $name): callable
    {
        return match ($name) {
            'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'exp', 'log', 'sqrt', 'atan' => $name,
            'abs' => static fn(float $v): float => abs($v),
            'sech' => static fn(float $v): float => 1.0 / cosh($v),
            'heaviside' => static fn(float $v): float => $v > 0 ? 1.0 : ($v < 0 ? 0.0 : 0.5),
            'sign' => static fn(float $v): float => (float)($v <=> 0),
            default => throw new \LogicException("function $name not allowed"),
        };
    }

    private static function scalar(string $op, float $a, float $b): float
    {
        return match ($op) {
            '+' => $a + $b,
            '-' => $a - $b,
            '*' => $a * $b,
            '/' => fdiv($a, $b),
            '^' => self::pow($a, $b),
        };
    }

    private static function pow(float $a, float $b): float
    {
        if ($a == 0.0 && $b < 0) {
            return INF;
        }
        return $a ** $b;
    }

    /** @param float|list<float> $a @param float|list<float> $b @return float|list<float> */
    private static function binop(string $op, $a, $b)
    {
        $aa = is_array($a);
        $ba = is_array($b);
        if (!$aa && !$ba) {
            return self::scalar($op, (float)$a, (float)$b);
        }
        $out = [];
        if ($aa && $ba) {
            switch ($op) {
                case '+': foreach ($a as $i => $v) { $out[$i] = $v + $b[$i]; } break;
                case '-': foreach ($a as $i => $v) { $out[$i] = $v - $b[$i]; } break;
                case '*': foreach ($a as $i => $v) { $out[$i] = $v * $b[$i]; } break;
                case '/': foreach ($a as $i => $v) { $out[$i] = fdiv($v, $b[$i]); } break;
                default: foreach ($a as $i => $v) { $out[$i] = self::pow($v, $b[$i]); }
            }
        } elseif ($aa) {
            $b = (float)$b;
            switch ($op) {
                case '+': foreach ($a as $i => $v) { $out[$i] = $v + $b; } break;
                case '-': foreach ($a as $i => $v) { $out[$i] = $v - $b; } break;
                case '*': foreach ($a as $i => $v) { $out[$i] = $v * $b; } break;
                case '/': foreach ($a as $i => $v) { $out[$i] = fdiv($v, $b); } break;
                default: foreach ($a as $i => $v) { $out[$i] = self::pow($v, $b); }
            }
        } else {
            $a = (float)$a;
            switch ($op) {
                case '+': foreach ($b as $i => $v) { $out[$i] = $a + $v; } break;
                case '-': foreach ($b as $i => $v) { $out[$i] = $a - $v; } break;
                case '*': foreach ($b as $i => $v) { $out[$i] = $a * $v; } break;
                case '/': foreach ($b as $i => $v) { $out[$i] = fdiv($a, $v); } break;
                default: foreach ($b as $i => $v) { $out[$i] = self::pow($a, $v); }
            }
        }
        return $out;
    }
}
