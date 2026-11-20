<?php
declare(strict_types=1);

namespace ManakovLab;

final class Http
{
    public static function json(array $data, int $status = 200): never
    {
        http_response_code($status);
        header('Content-Type: application/json; charset=utf-8');
        header('Cache-Control: no-store');
        header('X-Content-Type-Options: nosniff');
        $body = json_encode($data, JSON_UNESCAPED_SLASHES | JSON_UNESCAPED_UNICODE | JSON_PRESERVE_ZERO_FRACTION);
        $acceptsGzip = str_contains($_SERVER['HTTP_ACCEPT_ENCODING'] ?? '', 'gzip');
        if ($acceptsGzip && strlen((string)$body) > 2048 && function_exists('gzencode') && !ini_get('zlib.output_compression')) {
            header('Content-Encoding: gzip');
            header('Vary: Accept-Encoding');
            $body = gzencode((string)$body, 6);
        }
        echo $body;
        exit;
    }

    public static function error(string $message, int $status = 400, array $extra = []): never
    {
        self::json(['ok' => false, 'error' => $message] + $extra, $status);
    }
}
