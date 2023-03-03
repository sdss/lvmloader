# !/usr/bin/env python
# -*- coding: utf-8 -*-
#

import pathlib
import itertools
import time

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

import numpy as np
from sqlalchemy import select, insert

from sdssdb.sqlalchemy.lvmdb import database, drp
from bundle import Bundle

class Timer:

    def __init__(self, message: str = '', verbose: bool = None):
        self.message = message
        self.verbose = verbose
        self.start_time = None
        self.end_time = None
        self.elapsed = None

    def __enter__(self):
        self.start_time = time.perf_counter()

    def __exit__(self, *exc_info):
        self.end_time = time.perf_counter()
        self.elapsed = self.end_time - self.start_time
        if self.verbose:
            print(f'{self.message} - Elapsed time [sec]: {self.elapsed}')


class LVMLoader:

    def __init__(self, filename: str, release: str = None, verbose: bool = None,
                 manga: bool = False):
        self.filepath = pathlib.Path(filename).resolve()
        self.filename = self.filepath.name
        self.drpver = None
        self.dapver = None
        self.release = release
        self.session = database.Session()
        self.verbose = verbose


    def __repr__(self):
        return f'<LVMLoader(file="{self.filename}")'

    def _close(self):
        self.session.close()

    def __del__(self):
        self._close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self._close()

    def load_into_db(self):

        print(f'Starting load: {self.filename}')
        with Timer(f'Loading file: {self.filename}', verbose=self.verbose):
            with fits.open(self.filepath) as hdu:
                self.set_params(hdu)

                with database.Session() as session, session.begin():
                    self.session = session
                    with Timer('Loading Metadata', verbose=self.verbose):
                        self.load_rss()
                        self.load_wavelength(hdu)
                        self.load_header(hdu)
                        self.load_obsinfo(hdu)
                with Timer('Loading pixels', verbose=self.verbose):
                    self.load_fibers(hdu)
        self._close()

    def set_params(self, hdu):
        # temporarily using MaNGA
        self.drpver = hdu[0].header.get('VERSDRP2')
        self.dapver = hdu[0].header.get('VERSDAP')
        self.tileid = int(hdu[0].header['PLATEIFU'].replace('-', ''))
        self.ra = hdu[0].header.get('IFURA')
        self.dec = hdu[0].header.get('IFUDEC')
        self.ring = 7 if str(hdu[0].header['IFUDSGN']).startswith('127') else 6

    def get_footprint(self):
        b = Bundle(ra=self.ra, dec=self.dec, ring=self.ring)
        return b.hexagon.to_table().as_array().tolist()

    def load_rss(self):
        q = select(drp.RSS).where(drp.RSS.tileid == self.tileid)
        rss = self.session.scalars(q).one_or_none()

        if not rss:
            foot = self.get_footprint()
            rss = drp.RSS(tileid=self.tileid, ifura=self.ra, ifudec=self.dec, footprint=foot)
            rss.pipeline = self.load_pipeline()
            self.session.add(rss)

        return rss

    def load_wavelength(self, hdu):

        if ww := self.session.scalars(select(drp.Wavelength)).all():
            return ww

        wave = hdu['WAVE'].data.tolist()
        nwave = len(wave)
        ww = [dict(zip(['index', 'value'], row)) for row in zip(np.arange(nwave).tolist(), wave)]
        self.session.execute(insert(drp.Wavelength), ww)

    def get_header(self, hdu):
        hdr = hdu['PRIMARY'].header
        rss = self.load_rss()
        return [dict(zip(['index', 'key', 'value', 'comment', 'rss_pk'], row))
                for row in zip(np.arange(len(hdr)).tolist(), hdr.keys(), hdr.values(),
                               hdr.comments, itertools.repeat(rss.pk))]

    def load_header(self, hdu):

        rss = self.load_rss()
        if self.session.scalars(select(drp.Header).where(drp.Header.rss==rss)).first():
            return
        data = self.get_header(hdu)
        self.session.execute(insert(drp.Header), data)


    def load_pipeline(self, label: str = 'DRP'):

        version = self.drpver if label == 'DRP' else self.dapver if label == 'DAP' else None

        q = select(drp.Pipeline).where(drp.Pipeline.version == version)
        pipe = self.session.scalars(q).one_or_none()

        if not pipe:
            pipe = drp.Pipeline(version=self.drpver, label=label, release=self.release)
            self.session.add(pipe)

        return pipe


    def get_obsinfo(self, hdu):
        rss = self.load_rss()

        cols = ['EXPNUM', 'MJD', 'EXPTIME', 'MGDPOS', 'MGDRA', 'MGDDEC',
                'SEEING', 'BLUESN2', 'REDSN2', 'DRP2QUAL']

        tt = Table(hdu['OBSINFO'].data)
        df = tt[cols].to_pandas()
        df.columns = df.columns.str.lower()

        dmap = {'C':0, 'N':1, 'E': 2, 'S':3}
        df['mgdpos'] = df['mgdpos'].apply(lambda x: dmap[x.strip()])
        df['expnum'] = df['expnum'].astype(int)
        df.insert(0, 'rss_pk', rss.pk)

        # cols = ['rss_pk', 'expnum', 'mjd', 'exptime', 'nexp', 'dpos',
        #         'dra', 'ddec', 'airmass', 'seeing', 'bluesn2', 'redsn2'
        #         'nirsn2', 'drp2qual']

        df = df.rename(columns={'mgdpos': 'dpos', 'mgdra':'dra', 'mgddec':'ddec'})

        return df.to_dict(orient='records')

    def load_obsinfo(self, hdu):
        rss = self.load_rss()
        if self.session.scalars(select(drp.ObsInfo).where(drp.ObsInfo.rss==rss)).first():
            return
        data = self.get_obsinfo(hdu)
        self.session.execute(insert(drp.ObsInfo), data)

    def convert_xy_radec(self, x, y):
        s = SkyCoord(self.ra, self.dec, unit=u.degree)
        radecs = s.spherical_offsets_by(x * u.arcsec, y * u.arcsec)
        return zip(*radecs.to_table().as_array().tolist())

    def get_fiber_data(self, hdu):
        xp = hdu['XPOS'].data
        yp = hdu['YPOS'].data
        nrows, nwave = xp.shape
        nexp = len(hdu['OBSINFO'].data)
        nfiber = int(nrows / nexp)
        wrange = np.arange(nwave).tolist()
        widx = int(np.median(wrange))

        rss = self.load_rss()
        rsspk = itertools.repeat(rss.pk)

        # TODO look up real fibers pk
        fibpk = np.tile(np.arange(nfiber) + 1, nexp).tolist()

        # look up obsinfo
        q = select(drp.ObsInfo.pk).where(drp.ObsInfo.rss == rss)
        if obsinfo := self.session.scalars(q).all():
            obspk = np.repeat(obsinfo, nfiber).tolist()
        else:
            obspk = np.repeat(None, nrows).tolist()

        yidx = np.arange(nrows).tolist()
        xpos = xp[:, wrange].ravel().tolist()
        ypos = yp[:, wrange].ravel().tolist()

        # get wave_pk
        wp = self.session.scalars(select(drp.Wavelength).where(drp.Wavelength.index==widx)).one()
        wave_pk = np.tile(wp.pk, nrows).tolist()

        # get ra/dec
        ra, dec = self.convert_xy_radec(xpos, ypos)

        cols = ['rss_pk', 'fibers_pk', 'obsinfo_pk',
                'yidx', 'xpos', 'ypos', 'wave_pk', 'ra',
                'dec']

        return [dict(zip(cols, row))
                for row in zip(rsspk, fibpk, obspk, yidx, xpos,
                               ypos, wave_pk, ra, dec)]

    def load_fibers(self, hdu):

        rss = self.load_rss()
        if self.session.scalars(select(drp.RSSFiber).where(drp.RSSFiber.rss==rss)).first():
            return

        data = self.get_fiber_data(hdu)

        with database.engine.begin() as conn:
            ins = drp.RSSFiber.__table__.insert()
            conn.execute(ins, data)


# files = list(pathlib.Path('/Users/Brian/Work/Manga/redux/v3_1_1/').rglob("*[9,11]*/*/*RSS*"))
def load_lvm_files(files):
    for file in files:
        ll = LVMLoader(file, verbose=True)
        ll.load_into_db()
