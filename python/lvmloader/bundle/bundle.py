# !usr/bin/env python
# -*- coding: utf-8 -*-
#
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from regions import PolygonSkyRegion


class Bundle(object):
    """The location, size, and shape of a LVM IFU bundle.

    A bundle of a given size at the RA/Dec coordinates.

    Setting envelope creates the hexagon at the outer boundry of the fibers,
    otherwise the hexagon runs through the centers of the outer fibers.

    Parameters:
        ra (float):
            The Right Ascension of the target
        dec (float):
            The Declination of the target
        ring (int):
            The ring size of the IFU, e.g. 25 = 1801 fibers
        use_envelope (bool):
            Expands the hexagon area to include an envelope surrounding the hex border
        local (bool):
            If True, grabs the fiducial metrology file from a local LVMCORE

    Attributes:
        fibers (SkyCoord):
            A SkyCoord array of IFU fiber coordinates
        skies (SkyCoord):
            An Nx2 array of sky fiber coordinates.  Default is unloaded.  Use :meth:`get_sky_coordinates` to load.
        hexagon (SkyCoord):
            A SkyCoord array of coordinates defining the hexagon vertices
        sky_region (PolygonSkyRegion):
            An Astropy SkyRegion
    """

    def __init__(self, ra: float, dec: float, ring: int = 25, use_envelope: bool = True,
                 local: bool = None, manga: bool = None, **kwargs):
        self.ra = ra
        self.dec = dec
        self.coord = SkyCoord(self.ra, self.dec, unit=u.degree)
        self.fiber_diameter = 35.3 * u.arcsec
        self.angle = 0.  # original rotated angle
        self.rotated = False  # rotated originally by -90 degrees
        self.manga = manga

        self.simbmap = None
        self.plateholes_file = None

        self.size = self.total_fiber_num(ring)
        self.rings = np.arange(1, ring + 1)
        self.totals = self.total_fiber_num(self.rings)

        self.hexagon = None
        self.sky_region = None

        if self.size not in self.totals:
            return

        self.simbmap = self._get_simbmap_file(local=local)

        # get the fiber coordinates
        # self.skies = None
        self.fibers = self.get_fiber_coordinates(size=self.size)

        # get the hexagon coordinates and create a sky region
        self.hexagon = self._calculate_hexagon(use_envelope=use_envelope)
        self.sky_region = PolygonSkyRegion(vertices=self.hexagon)

    def __repr__(self):
        return '<Bundle (ra={0}, dec={1}, ifu={2})>'.format(self.ra, self.dec, self.size)

    @staticmethod
    def total_fiber_num(ring: int) -> int:
        """ compute the total fiber number within a ring id """
        return 3 * ((ring - 1) ** 2) + 3 * (ring - 1) + 1

    def _get_simbmap_file(self, local: bool = None) -> pd.DataFrame:
        ''' Retrieves the metrology file locally or remotely

        Reads in fresh metrology data or grabs it from
        cache if available

        Parameters:
            local (bool):
                If True, does a local system check

        Returns:
            The metrology object data
        '''

        # create the file path
        if self.manga:
            rel_path = u'metrology/manga_simbmap_127.dat'
        else:
            rel_path = u'metrology/lvm_simbmap_1801.dat'
        if local:
            base_path = os.environ['LVMCORE_DIR']
        else:
            base_path = u'https://raw.githubusercontent.com/sdss/lvmcore/master'

        self.simbmap_file = os.path.join(base_path, rel_path)
        return pd.read_table(self.simbmap_file, sep=',', header=0 if self.manga else 21)

    # def _get_plateholes_file(self, plate=None, local=None):
    #     ''' Retrieves the platesholes file locally or remotely

    #     Reads in fresh plateHoles data or grabs it from
    #     cache if available

    #     Parameters:
    #         plate (int):
    #             The plateid to load the plateHoles file for
    #         local (bool):
    #             If True, does a local system check

    #     Returns:
    #         The plateHoles object data
    #     '''

    #     # check for plate
    #     if not plate and not self.plateifu:
    #         raise MarvinError('Must have the plate id to get the correct plateholes file')
    #     elif not plate and self.plateifu:
    #         plate, ifu = self.plateifu.split('-')

    #     # create the file path
    #     pltgrp = '{:04d}XX'.format(int(plate) // 100)
    #     rel_path = u'platedesign/plateholes/{0}/plateHolesSorted-{1:06d}.par'.format(pltgrp, int(plate))
    #     if local:
    #         base_path = os.environ['MANGACORE_DIR']
    #     else:
    #         base_path = u'https://svn.sdss.org/public/repo/manga/mangacore/tags/v1_2_3'

    #     self.platesholes_file = os.path.join(base_path, rel_path)

    #     # read in the Yanny object
    #     global PLATEHOLES
    #     if PLATEHOLES is None:
    #         plateholes_data = self._read_in_yanny(self.platesholes_file, local=local)
    #         PLATEHOLES = plateholes_data['STRUCT1']

    #     return PLATEHOLES

    # @staticmethod
    # def _get_sky_fibers(data, ifu):
    #     """ Return the sky fibers associated with each galaxy.

    #     Parameters:
    #         data (object):
    #             the plateHoles object data
    #         ifu (int):
    #             The ifu design of the target

    #     Returns:
    #         A numpy array of tuples of sky fiber RA, Dec coordinates
    #     """
    #     sky_fibers = data[data['targettype'] == 'SKY']
    #     sky_block = sky_fibers['block'] == int(ifu)
    #     zip_data = list(zip(sky_fibers[sky_block]['target_ra'], sky_fibers[sky_block]['target_dec']))
    #     skies = np.array(zip_data)
    #     return skies

    # def get_sky_coordinates(self, plateifu=None):
    #     ''' Returns the RA, Dec coordinates of the sky fibers

    #     An Nx2 numpy array, with each row the [RA, Dec] coordinate of
    #     the sky fiber.

    #     Parameters:
    #         plateifu (str):
    #             The plateifu of the target
    #     '''

    #     plateifu = plateifu or self.plateifu
    #     if not plateifu:
    #         raise MarvinError('Need the plateifu to get the corresponding sky fibers')
    #     else:
    #         plate, ifu = plateifu.split('-')

    #     data = self._get_plateholes_file(plate=plate)
    #     self.skies = self._get_sky_fibers(data, ifu)

    def get_fiber_coordinates(self, size: int = None):
        """ Returns the RA, Dec coordinates for each fiber.

        Parameters:
            size (int):
                The IFU size.  This extracts only the fibers for the given size

        Returns:
            An Nx3 numpy array, each row containing [fiberid, RA, Dec]
        """

        fibers = self.coord.spherical_offsets_by(
            self.simbmap['raoff']*u.arcsec, self.simbmap['decoff']*u.arcsec)

        # extract only the fibers for the specified IFU size
        size = size or self.size
        fibers = fibers[:size]

        return fibers

    def _calculate_hexagon(self, use_envelope: bool = True, rotated: bool = True):
        """ Calculates the vertices of the bundle hexagon. """

        simbmap = self.simbmap[:self.size].copy()

        # assign the row axes based on the original angle
        row, col = ('raoff', 'decoff') if self.angle != 0 else ('decoff', 'raoff')

        # identify middle, top, bottom row slices
        middleRow = simbmap[simbmap[row].round(2) == 0]
        topRow = simbmap[simbmap[row].round(2).max() - simbmap[row].round(2) == 0]
        bottopRow = simbmap[simbmap[row].round(2).min() - simbmap[row].round(2) == 0]

        if len(topRow) < max(self.rings) or len(bottopRow) < max(self.rings):
            print('Warning, did not select all fibers in top/bottom row' )

        # identify hexagon vertices
        vertice0 = middleRow[middleRow[col] == np.max(middleRow[col])]
        vertice3 = middleRow[middleRow[col] == np.min(middleRow[col])]

        vertice1 = topRow[topRow[col] == np.max(topRow[col])]
        vertice2 = topRow[topRow[col] == np.min(topRow[col])]

        vertice5 = bottopRow[bottopRow[col] == np.max(bottopRow[col])]
        vertice4 = bottopRow[bottopRow[col] == np.min(bottopRow[col])]

        hexagonOff = np.array(
            [[vertice0['raoff'].values[0], vertice0['decoff'].values[0]],
             [vertice1['raoff'].values[0], vertice1['decoff'].values[0]],
             [vertice2['raoff'].values[0], vertice2['decoff'].values[0]],
             [vertice3['raoff'].values[0], vertice3['decoff'].values[0]],
             [vertice4['raoff'].values[0], vertice4['decoff'].values[0]],
             [vertice5['raoff'].values[0], vertice5['decoff'].values[0]]])

        # This array increases the area of the hexagon so that it is an
        # envelope of the bundle.
        if use_envelope:
            halfFibre = self.fiber_diameter.value / 2.
            hexagonExtra = np.array(
                [[halfFibre, 0.0],
                 [halfFibre, halfFibre],
                 [-halfFibre, halfFibre],
                 [-halfFibre, 0.0],
                 [-halfFibre, -halfFibre],
                 [halfFibre, -halfFibre]])

            if self.angle != 0:
                angle_rad = self.angle * np.pi / 180
                rot = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],[-np.sin(angle_rad), np.cos(angle_rad)]])
                hexagonExtra = np.dot(hexagonExtra, rot)

            hexagonOff += hexagonExtra

        r, d = zip(*hexagonOff)
        return self.coord.spherical_offsets_by(r*u.arcsec, d*u.arcsec)

    def _rotate(self, angle: int = 90):
        data = self.simbmap[['raoff', 'decoff']].copy()
        angle_rad = angle * np.pi / 180
        rot = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],[-np.sin(angle_rad), np.cos(angle_rad)]])
        rotated = np.dot(data, rot)
        r, d = zip(*rotated)
        self.simbmap.raoff = r
        self.simbmap.decoff = d
        self.angle += angle
        self.hexagon = self._calculate_hexagon()
        self.sky_region = PolygonSkyRegion(vertices=self.hexagon)


    def plot(self, hdu: fits.ImageHDU = None, wcs: WCS = None, nx: int = 100,
             scale: float = 35.3, sunit: str = 'arcsec', degaxis: bool = True,
             overlay_fibers: bool = False, **kwargs):

        if hdu and not wcs:
            wcs = WCS(hdu.header)
        elif wcs and not hdu:
            ny = nx
            hdu = fits.ImageHDU(np.ones([nx, ny]))
        elif not hdu:
            ny = nx
            hdu = fits.ImageHDU(np.ones([nx, ny]))
            wcs_dict = {
                'CTYPE1': 'RA---TAN', 'CUNIT1': 'deg', 'CD1_1': -2.77778e-4*scale, 'CRPIX1': nx//2 + 1, 'CRVAL1': self.ra, 'NAXIS1': nx,
                'CTYPE2': 'DEC--TAN', 'CUNIT2': 'deg', 'CD2_2': 2.77778e-4*scale, 'CRPIX2': ny//2 + 1, 'CRVAL2': self.dec, 'NAXIS2': ny}
            wcs = WCS(wcs_dict)

        plt.clf()
        if wcs.has_celestial:
            ww = wcs.celestial
        ax = plt.subplot(projection=ww)

        if degaxis:
            ax.coords[0].set_format_unit(u.degree, decimal=True)
            ax.coords[1].set_format_unit(u.degree, decimal=True)

        ax.imshow(hdu.data, **kwargs)

        # plot hexagon
        pp = self.sky_region.to_pixel(ww)
        pp.plot(ax=ax, edgecolor=kwargs.get('edgecolor', 'red'))

        # plot fibers
        if overlay_fibers:
            for fiber in self.simbmap[['raoff', 'decoff']][:self.size].iterrows():
                fiber_coord = self.coord.spherical_offsets_by(fiber[1]['raoff']*u.arcsec, fiber[1]['decoff']*u.arcsec)
                p = SphericalCircle((fiber_coord.ra, fiber_coord.dec), scale / 2 * getattr(u, sunit),
                                    edgecolor='green', facecolor='none',
                                    transform=ax.get_transform('fk5'))
                ax.add_patch(p)

        return ax

    def create_ds9_region(self):
        return self.sky_region.serialize('ds9')

    def print_bundle(self):
        ''' Print the bundle to an Astropy Table '''

        tt = self.fibers.to_table()
        tt.add_column(np.arange(self.size) + 1, index=0, name='fiberid')
        ascii.write(tt, format='fixed_width_two_line',
                    formats={'ra': '{0:.12f}', 'dec': '{0:.12f}'})

    def print_hexagon(self):
        ''' Print the hexagon to an Astropy Table '''

        tt = self.hexagon.to_table()
        ascii.write(tt, format='fixed_width_two_line',
                    formats={'ra': '{0:.12f}', 'dec': '{0:.12f}'})


# plot a series of hex footprints or regions on a blank wcs

# feet=session.scalars(select(drp.RSS.footprint)).all()

# used only to generate a blank wcs
# b=Bundle(ra=ra, dec=dec, ring=7)
# ax=b.plot(nx=100)

# color = iter(plt.cm.rainbow(np.linspace(0, 1, n)))
# for foot in feet:
#     reg = PolygonSkyRegion(vertices=SkyCoord(foot, unit=u.degree))
#     pp = reg.to_pixel(ax.wcs)
#     pp.plot(ax=ax, edgecolor=next(color))