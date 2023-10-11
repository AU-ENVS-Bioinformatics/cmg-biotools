DROP TABLE IF EXISTS `cache`;
CREATE TABLE `cache` (
 `id` varchar(128) NOT NULL DEFAULT '',
 `source` varchar(256) NOT NULL DEFAULT '',
 `content` longblob,
 `expire` datetime DEFAULT NULL,
 `stamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
 `usg` bigint(20) DEFAULT '0',
 PRIMARY KEY (`id`,`source`),
 KEY `expire` (`expire`)
);