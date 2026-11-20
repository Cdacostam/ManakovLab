<?php
declare(strict_types=1);

namespace ManakovLab;

/**
 * Limits how many simulations run at the same time, using exclusive file locks
 * ("slots") in the system temp directory. A slot is released automatically when
 * the PHP process ends, even after a fatal error.
 */
final class ConcurrencyGate
{
    /** @var resource|null */
    private $handle = null;

    public function __construct(private int $slots)
    {
    }

    public function acquire(int $waitSeconds = 30): bool
    {
        $deadline = microtime(true) + $waitSeconds;
        do {
            for ($slot = 0; $slot < max(1, $this->slots); $slot++) {
                $fh = @fopen(sys_get_temp_dir() . "/manakov-lab-slot-$slot.lock", 'c');
                if ($fh && flock($fh, LOCK_EX | LOCK_NB)) {
                    $this->handle = $fh;
                    return true;
                }
                if ($fh) {
                    fclose($fh);
                }
            }
            usleep(250_000);
        } while (microtime(true) < $deadline);
        return false;
    }

    public function release(): void
    {
        if ($this->handle) {
            flock($this->handle, LOCK_UN);
            fclose($this->handle);
            $this->handle = null;
        }
    }

    public function __destruct()
    {
        $this->release();
    }
}
