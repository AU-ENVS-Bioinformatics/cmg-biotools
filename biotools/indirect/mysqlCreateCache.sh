#!/bin/sh
mysqladmin -u root -pPrinter1 create cache
mysql -u root -pPrinter1 -q cache < /usr/biotools/indirect/cacher.sql 
mysql -u root -pPrinter1 -q cache -e "GRANT ALL ON cache.* to ''@'localhost'"