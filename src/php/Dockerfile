FROM php:8.3-apache

# Pure-PHP numerical core: no Python or other runtime needed.
# OPcache JIT gives the solvers a speed-up of about 3x.
RUN docker-php-ext-install opcache \
 && printf '%s\n' \
    'opcache.enable=1' \
    'opcache.jit=tracing' \
    'opcache.jit_buffer_size=128M' \
    'max_execution_time=180' \
    'memory_limit=512M' \
    > /usr/local/etc/php/conf.d/manakov.ini

ENV APACHE_DOCUMENT_ROOT=/var/www/manakov-lab
RUN sed -ri 's!/var/www/html!${APACHE_DOCUMENT_ROOT}!g' /etc/apache2/sites-available/*.conf /etc/apache2/apache2.conf \
 && sed -ri 's!AllowOverride None!AllowOverride All!g' /etc/apache2/apache2.conf

COPY . /var/www/manakov-lab
RUN chown -R www-data:www-data /var/www/manakov-lab

EXPOSE 80
