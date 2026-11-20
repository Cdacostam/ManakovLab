<?php
declare(strict_types=1);

if (realpath($_SERVER['SCRIPT_FILENAME'] ?? '') === __FILE__) { http_response_code(404); exit; }

/*
 * Manakov Lab configuration. Every value can be overridden with an
 * environment variable of the same name (e.g. in the web-server config).
 *
 * The numerical core is pure PHP. Enable the OPcache JIT for a speed-up of about 3x (2.9-3.4x on the presets):
 *   opcache.enable=1, opcache.jit=tracing, opcache.jit_buffer_size=128M
 */
return [
    // Hard wall-clock limit for one simulation (seconds)
    'MANAKOV_TIMEOUT' => (int)(getenv('MANAKOV_TIMEOUT') ?: 150),
    // Number of simulations allowed to run at the same time
    'MANAKOV_MAX_CONCURRENT' => (int)(getenv('MANAKOV_MAX_CONCURRENT') ?: 2),
    // Work budgets that reject runs which would take too long.
    // Split-step: Nx * log2(Nx) * (T/dt)   (~3.5e6 per second without JIT)
    'MANAKOV_MAX_WORK_SPLIT' => (float)(getenv('MANAKOV_MAX_WORK_SPLIT') ?: 4.0e8),
    // Caputo L1: N^2 * Nx                   (~1e7 per second without JIT)
    'MANAKOV_MAX_WORK_CAPUTO' => (float)(getenv('MANAKOV_MAX_WORK_CAPUTO') ?: 8.0e8),
    // PHP memory limit applied to simulation requests
    'MANAKOV_MEMORY_LIMIT' => getenv('MANAKOV_MEMORY_LIMIT') ?: '512M',
    // Plotly: a local copy is used if this file exists, otherwise the CDN below
    'PLOTLY_LOCAL' => 'assets/vendor/plotly/plotly.min.js',
    'PLOTLY_CDN' => 'https://cdn.plot.ly/plotly-2.35.2.min.js',
    // Links shown in the page footers (a value still containing xxxx is treated as not yet available on the home page)
    'REPOSITORY_URL' => getenv('MANAKOV_REPOSITORY_URL') ?: 'https://github.com/Cdacostam/ManakovLab',
    'CITATION_URL' => getenv('MANAKOV_CITATION_URL') ?: 'https://doi.org/10.5281/zenodo.24406892',
    // Zenodo DOI of the archived release, shown on the home page (bare DOI, without https://doi.org/)
    'DOI' => getenv('MANAKOV_DOI') ?: '10.5281/zenodo.24406892',
];
