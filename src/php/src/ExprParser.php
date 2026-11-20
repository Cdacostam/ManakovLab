<?php
declare(strict_types=1);

namespace ManakovLab;

/**
 * Recursive-descent parser for coefficient expressions such as
 *     1 + 0.5 cos(2t),   -0.01 x^2,   beta(t) = 0.05*sin(t)*exp(-t/10)
 *
 * Produces an AST that is
 *   - validated against a whitelist (variables, constants, functions),
 *   - rendered to LaTeX for KaTeX,
 *   - compiled to closures by Numerics\ExprCompiler and serialised to a canonical form.
 *
 * Grammar (implicit multiplication allowed: "2t", "0.5cos(t)", "2(x+1)"):
 *   expr    := term (('+'|'-') term)*
 *   term    := unary (('*'|'/') unary | <implicit> unary)*
 *   unary   := ('+'|'-') unary | power
 *   power   := primary (('^'|'**') unary)?          (right associative)
 *   primary := NUMBER | IDENT | IDENT '(' expr ')' | '(' expr ')'
 */
final class ExprParser
{
    public const FUNCTIONS = [
        'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'sech',
        'exp', 'log', 'sqrt', 'abs', 'atan', 'heaviside', 'sign',
    ];
    public const CONSTANTS = ['pi', 'e'];
    private const MAX_LENGTH = 300;
    private const MAX_DEPTH = 60;
    private const MAX_NODES = 250;

    /** @var list<array{type:string,value:string,pos:int}> */
    private array $tokens = [];
    private int $i = 0;
    private int $depth = 0;
    private int $nodes = 0;
    /** @var list<string> */
    private array $vars;

    /** @param list<string> $vars allowed independent variables, e.g. ['t'] or ['x','t'] */
    public function __construct(array $vars)
    {
        $this->vars = $vars;
    }

    /**
     * Parse and validate. Returns the AST.
     * @throws ExprException
     */
    public function parse(string $src): array
    {
        $src = trim($src);
        if ($src === '') {
            $src = '0';
        }
        if (strlen($src) > self::MAX_LENGTH) {
            throw new ExprException('Expression is too long (max ' . self::MAX_LENGTH . ' characters).', 0);
        }
        $this->tokens = $this->tokenize($src);
        $this->i = 0;
        $this->depth = 0;
        $this->nodes = 0;
        $ast = $this->expr();
        if ($this->peek()['type'] !== 'eof') {
            $t = $this->peek();
            throw new ExprException("Unexpected '{$t['value']}'.", $t['pos']);
        }
        return $ast;
    }

    /** Variables actually used by an AST. */
    public static function usedVariables(array $n): array
    {
        $out = [];
        $walk = function (array $n) use (&$walk, &$out): void {
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

    // ------------------------------------------------------------------ lexer
    private function tokenize(string $s): array
    {
        $tokens = [];
        $len = strlen($s);
        $p = 0;
        while ($p < $len) {
            $c = $s[$p];
            if (ctype_space($c)) { $p++; continue; }
            if (preg_match('/\G(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?/A', $s, $m, 0, $p)) {
                $tokens[] = ['type' => 'num', 'value' => $m[0], 'pos' => $p];
                $p += strlen($m[0]);
                continue;
            }
            if (preg_match('/\G[A-Za-z_][A-Za-z_0-9]*/A', $s, $m, 0, $p)) {
                $tokens[] = ['type' => 'id', 'value' => $m[0], 'pos' => $p];
                $p += strlen($m[0]);
                continue;
            }
            if ($c === '*' && $p + 1 < $len && $s[$p + 1] === '*') {
                $tokens[] = ['type' => 'op', 'value' => '^', 'pos' => $p];
                $p += 2;
                continue;
            }
            if (strpos('+-*/^()', $c) !== false) {
                $tokens[] = ['type' => 'op', 'value' => $c, 'pos' => $p];
                $p++;
                continue;
            }
            throw new ExprException("Invalid character '" . $c . "'.", $p);
        }
        $tokens[] = ['type' => 'eof', 'value' => 'end of input', 'pos' => $len];
        return $tokens;
    }

    // ----------------------------------------------------------------- parser
    private function peek(): array { return $this->tokens[$this->i]; }
    private function next(): array { return $this->tokens[$this->i++]; }
    private function isOp(string $v): bool
    {
        $t = $this->peek();
        return $t['type'] === 'op' && $t['value'] === $v;
    }

    private function node(array $n): array
    {
        if (++$this->nodes > self::MAX_NODES) {
            throw new ExprException('Expression is too complex.', $this->peek()['pos']);
        }
        return $n;
    }

    private function enter(): void
    {
        if (++$this->depth > self::MAX_DEPTH) {
            throw new ExprException('Expression is nested too deeply.', $this->peek()['pos']);
        }
    }

    private function expr(): array
    {
        $this->enter();
        $l = $this->term();
        while ($this->isOp('+') || $this->isOp('-')) {
            $op = $this->next()['value'];
            $l = $this->node(['type' => 'bin', 'op' => $op, 'l' => $l, 'r' => $this->term()]);
        }
        $this->depth--;
        return $l;
    }

    private function term(): array
    {
        $l = $this->unary();
        while (true) {
            if ($this->isOp('*') || $this->isOp('/')) {
                $op = $this->next()['value'];
                $l = $this->node(['type' => 'bin', 'op' => $op, 'l' => $l, 'r' => $this->unary()]);
                continue;
            }
            $t = $this->peek();
            if ($t['type'] === 'num' || $t['type'] === 'id' || ($t['type'] === 'op' && $t['value'] === '(')) {
                $l = $this->node(['type' => 'bin', 'op' => '*', 'l' => $l, 'r' => $this->unary(), 'implicit' => true]);
                continue;
            }
            return $l;
        }
    }

    private function unary(): array
    {
        $this->enter();
        if ($this->isOp('-')) {
            $this->next();
            $r = $this->node(['type' => 'neg', 'e' => $this->unary()]);
        } elseif ($this->isOp('+')) {
            $this->next();
            $r = $this->unary();
        } else {
            $r = $this->power();
        }
        $this->depth--;
        return $r;
    }

    private function power(): array
    {
        $base = $this->primary();
        if ($this->isOp('^')) {
            $this->next();
            return $this->node(['type' => 'bin', 'op' => '^', 'l' => $base, 'r' => $this->unary()]);
        }
        return $base;
    }

    private function primary(): array
    {
        $t = $this->next();
        if ($t['type'] === 'num') {
            if (!is_finite((float)$t['value'])) {
                throw new ExprException('Number out of range.', $t['pos']);
            }
            return $this->node(['type' => 'num', 'v' => $t['value']]);
        }
        if ($t['type'] === 'id') {
            $name = $t['value'];
            if (in_array($name, self::FUNCTIONS, true)) {
                if (!$this->isOp('(')) {
                    throw new ExprException("Function '$name' needs parentheses, e.g. $name(t).", $t['pos']);
                }
                $this->next();
                $arg = $this->expr();
                if (!$this->isOp(')')) {
                    throw new ExprException("Missing ')' after argument of $name.", $this->peek()['pos']);
                }
                $this->next();
                return $this->node(['type' => 'call', 'fn' => $name, 'arg' => $arg]);
            }
            if (in_array($name, $this->vars, true)) {
                return $this->node(['type' => 'var', 'name' => $name]);
            }
            if (in_array($name, self::CONSTANTS, true)) {
                return $this->node(['type' => 'const', 'name' => $name]);
            }
            $allowed = implode(', ', $this->vars);
            throw new ExprException("Unknown symbol '$name'. Allowed variables: $allowed; constants: pi, e.", $t['pos']);
        }
        if ($t['type'] === 'op' && $t['value'] === '(') {
            $e = $this->expr();
            if (!$this->isOp(')')) {
                throw new ExprException("Missing ')'.", $this->peek()['pos']);
            }
            $this->next();
            return $e;
        }
        throw new ExprException("Unexpected '{$t['value']}'.", $t['pos']);
    }

    // ------------------------------------------------------------ serialisers
    /** Canonical, fully parenthesised form (shown in the run report). */
    public static function toCanonical(array $n): string
    {
        switch ($n['type']) {
            case 'num': return $n['v'];
            case 'var': return $n['name'];
            case 'const': return $n['name'];
            case 'call': return $n['fn'] . '(' . self::toCanonical($n['arg']) . ')';
            case 'neg': return '(-' . self::toCanonical($n['e']) . ')';
            case 'bin':
                $op = $n['op'] === '^' ? '**' : $n['op'];
                return '(' . self::toCanonical($n['l']) . $op . self::toCanonical($n['r']) . ')';
        }
        throw new \LogicException('bad node');
    }

    private static function prec(array $n): int
    {
        if ($n['type'] === 'bin') {
            return ['+' => 1, '-' => 1, '*' => 2, '/' => 5, '^' => 4][$n['op']];
        }
        if ($n['type'] === 'neg') {
            return 3;
        }
        if ($n['type'] === 'num' && (str_contains($n['v'], 'e') || str_contains($n['v'], 'E'))) {
            return 2; // a \times 10^{b}
        }
        return 6;
    }

    private static function wrap(string $s): string
    {
        return '\\left(' . $s . '\\right)';
    }

    private static function startsWithDigit(string $latex): bool
    {
        return $latex !== '' && (ctype_digit($latex[0]) || $latex[0] === '.');
    }

    public static function toLatex(array $n): string
    {
        switch ($n['type']) {
            case 'num':
                if (preg_match('/^([\d.]+)[eE]([+-]?\d+)$/', $n['v'], $m)) {
                    $exp = (string)(int)$m[2];
                    return ($m[1] === '1' ? '' : $m[1] . '\\times ') . '10^{' . $exp . '}';
                }
                return $n['v'];
            case 'var':
                return $n['name'];
            case 'const':
                return $n['name'] === 'pi' ? '\\pi' : '\\mathrm{e}';
            case 'call':
                $a = self::toLatex($n['arg']);
                switch ($n['fn']) {
                    case 'sqrt': return '\\sqrt{' . $a . '}';
                    case 'abs': return '\\left|' . $a . '\\right|';
                    case 'exp': return '\\mathrm{e}^{' . $a . '}';
                    case 'log': return '\\ln' . self::wrap($a);
                    case 'atan': return '\\arctan' . self::wrap($a);
                    case 'sech': return '\\operatorname{sech}' . self::wrap($a);
                    case 'sign': return '\\operatorname{sgn}' . self::wrap($a);
                    case 'heaviside': return 'H' . self::wrap($a);
                    default: return '\\' . $n['fn'] . self::wrap($a);
                }
            case 'neg':
                $e = self::toLatex($n['e']);
                return '-' . (self::prec($n['e']) <= 1 || $n['e']['type'] === 'neg' ? self::wrap($e) : $e);
            case 'bin':
                $l = self::toLatex($n['l']);
                $r = self::toLatex($n['r']);
                switch ($n['op']) {
                    case '+':
                        if ($n['r']['type'] === 'neg') {
                            $r = self::wrap($r);
                        }
                        return $l . ' + ' . $r;
                    case '-':
                        if (self::prec($n['r']) <= 1 || $n['r']['type'] === 'neg') {
                            $r = self::wrap($r);
                        }
                        return $l . ' - ' . $r;
                    case '*':
                        if (self::prec($n['l']) < 2) {
                            $l = self::wrap($l);
                        }
                        if (self::prec($n['r']) <= 3) {
                            $r = self::wrap($r);
                        }
                        $sep = self::startsWithDigit($r) ? ' \\cdot ' : '\\,';
                        return $l . $sep . $r;
                    case '/':
                        return '\\frac{' . $l . '}{' . $r . '}';
                    case '^':
                        $baseNeedsParens = self::prec($n['l']) < 6 || $n['l']['type'] === 'call'
                            || ($n['l']['type'] === 'num' && str_contains($l, '.'));
                        if ($n['l']['type'] === 'call' && in_array($n['l']['fn'], ['sqrt', 'abs'], true)) {
                            $baseNeedsParens = false;
                        }
                        return ($baseNeedsParens ? self::wrap($l) : '{' . $l . '}') . '^{' . $r . '}';
                }
        }
        throw new \LogicException('bad node');
    }
}

final class ExprException extends \RuntimeException
{
    public int $position;

    public function __construct(string $message, int $position)
    {
        parent::__construct($message);
        $this->position = $position;
    }
}
