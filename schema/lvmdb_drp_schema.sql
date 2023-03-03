

CREATE SCHEMA drp;

SET search_path TO drp;


CREATE TABLE drp.pipeline (pk serial PRIMARY KEY NOT NULL,
    version TEXT, label TEXT, release TEXT);


CREATE TABLE drp.fibers (pk serial PRIMARY KEY NOT NULL, fiberid INTEGER, specid INTEGER,
    blockid TEXT, finblock INTEGER, targettype TEXT, ifulabel TEXT, finifu INTEGER,
    xpmm real, ypmm real, ringnum INTEGER, status INTEGER ,
    ifu_pk INTEGER);


CREATE TABLE drp.ifu (pk serial PRIMARY KEY NOT NULL, fiberid INTEGER, ringid INTEGER,
    ringfibnum INTEGER, specfibid INTEGER, specid INTEGER, label TEXT, standard TEXT,
    skyA TEXT, skyB TEXT, xpmm real, ypmm real, raoff double precision,
    decoff double precision);


CREATE TABLE drp.wavelength(pk serial PRIMARY KEY NOT NULL, index INTEGER,
    value real);


CREATE TABLE drp.header (pk serial PRIMARY KEY NOT NULL,
	index INTEGER, key TEXT, value TEXT, comment TEXT, rss_pk INTEGER);


CREATE TABLE drp.rssfiber (pk serial PRIMARY KEY NOT NULL, rss_pk INTEGER, fibers_pk INTEGER,
    obsinfo_pk INTEGER, yidx INTEGER, xpos real, ypos real, wave_pk INTEGER,
    ra double precision, dec double precision);

CREATE TABLE drp.obsinfo (pk serial PRIMARY KEY NOT NULL, rss_pk INTEGER,
    mjd INTEGER, expnum INTEGER, exptime real, nexp INTEGER, dpos INTEGER,
    dra double precision, ddec double precision, airmass real, seeing real,
    bluesn2 real, redsn2 real, nirsn2 real, drp2qual INTEGER);


CREATE TABLE drp.rss (pk serial PRIMARY KEY NOT NULL, tileid INTEGER, ifura double precision,
    ifudec double precision, pipeline_pk INTEGER, footprint polygon);


ALTER TABLE ONLY drp.obsinfo
    ADD CONSTRAINT rss_fk
    FOREIGN KEY (rss_pk) REFERENCES drp.rss(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.rssfiber
    ADD CONSTRAINT rss_fk
    FOREIGN KEY (rss_pk) REFERENCES drp.rss(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.rssfiber
    ADD CONSTRAINT obsinfo_fk
    FOREIGN KEY (obsinfo_pk) REFERENCES drp.obsinfo(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.rssfiber
    ADD CONSTRAINT fibers_fk
    FOREIGN KEY (fibers_pk) REFERENCES drp.fibers(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.rssfiber
    ADD CONSTRAINT wave_fk
    FOREIGN KEY (wave_pk) REFERENCES drp.wavelength(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.fibers
    ADD CONSTRAINT ifu_fk
    FOREIGN KEY (ifu_pk) REFERENCES drp.ifu(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.header
    ADD CONSTRAINT rss_fk
    FOREIGN KEY (rss_pk) REFERENCES drp.rss(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY drp.rss
    ADD CONSTRAINT pipeline_fk
    FOREIGN KEY (pipeline_pk) REFERENCES drp.pipeline(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;


/* Indexes */
CREATE INDEX CONCURRENTLY pipeline_pk_idx ON drp.rss using BTREE(pipeline_pk);
CREATE INDEX CONCURRENTLY ifu_pk_idx ON drp.fibers using BTREE(ifu_pk);
CREATE INDEX CONCURRENTLY obsinfo_rss_pk_idx ON drp.obsinfo using BTREE(rss_pk);
CREATE INDEX CONCURRENTLY rssfibers_rss_pk_idx ON drp.rssfiber using BTREE(rss_pk);
CREATE INDEX CONCURRENTLY rssfibers_fibers_pk_idx ON drp.rssfiber using BTREE(fibers_pk);
CREATE INDEX CONCURRENTLY rssfibers_wave_pk_idx ON drp.rssfiber using BTREE(wave_pk);
CREATE INDEX CONCURRENTLY rssfibers_obsinfo_pk_idx ON drp.rssfiber using BTREE(obsinfo_pk);
CREATE INDEX CONCURRENTLY fibers_pk_idx ON drp.rssfiber using BTREE(fibers_pk);
CREATE INDEX CONCURRENTLY header_rss_pk_idx ON drp.header using BTREE(rss_pk);
