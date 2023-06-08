#import packages used for generating subject set
import matplotlib.pyplot as plt
import gc
import numpy as np
import pandas
import pandas as pd
import time
import os

# Astropy imports
from astropy.wcs import WCS
from astropy.visualization import make_lupton_rgb
from astropy import units as u
from astropy.coordinates import SkyCoord

# Import the Rubin TAP service utilities
from lsst.rsp import get_tap_service, retrieve_query

# Image visualization routines.
import lsst.afw.display as afwDisplay
# The Butler provides programmatic access to LSST data products.
from lsst.daf.butler import Butler
# Geometry package
import lsst.geom as geom
# Object for multi-band exposures
from lsst.afw.image import MultibandExposure

import lsst.daf.butler as dafButler
import lsst.geom
import lsst.afw.display as afwDisplay

# Must explicitly set this to save figures
afwDisplay.setDefaultBackend("matplotlib")

plt.style.use('tableau-colorblind10')

import warnings
from astropy.units import UnitsWarning

pd.set_option('display.max_rows', 20)
warnings.simplefilter("ignore", category=UnitsWarning)

plot_filter_labels = {'u':'u', 'g':'g', 'r':'r', 'i':'i', 'z':'z', 'y':'y'}
plot_filter_colors = {'u': '#56b4e9', 'g': '#008060', 'r': '#ff4000',
                      'i': '#850000', 'z': '#6600cc', 'y': '#000000'}
plot_filter_symbols = {'u': 'o', 'g': '^', 'r': 'v', 'i': 's', 'z': '*', 'y': 'p'}
    
def get_cutout_image(butler, ra_deg, dec_deg, visit, detector, band, cutoutSideLength, datasetType='calexp'):
    """
    Get the cutout image information from butler. 
    This shoudl be followed by make_fig

    Input Parameters
    ----------
    ra : ra of source in degrees
    dec : dec of source in degrees
    visit : visit id
    detector : detector number
    band : band to get cutput for
    cutoutSideLength : size of the cutout
    
    Returns
    ----------
    Cutout image information
    """
    cutoutSize = geom.ExtentI(cutoutSideLength, cutoutSideLength)
    
    radec = geom.SpherePoint(ra_deg,dec_deg, geom.degrees)
    
    dataId = {'visit': visit, 'detector': detector}  
    calexp_wcs = butler.get('calexp.wcs', **dataId)
    
    xy = geom.PointI(calexp_wcs.skyToPixel(radec))
    bbox = geom.BoxI(xy - cutoutSize // 2, cutoutSize)
    parameters = {'bbox': bbox}
    
    cutout_image = butler.get('calexp', parameters=parameters, **dataId)

    return cutout_image

def make_calexp_fig(cutout_image, ra, dec, out_name):
    """
    Create an image.
    should be followed with remove_figure
    
    Parameters
    ----------
    cutout_image : cutout_image from butler.get
    ra : ra of source in degrees
    dec : dec of source in degrees
    out_name : file name where you'd like to save it
    
    Returns
    ----------
    cutout image
    """
#     fig = plt.figure(figsize=(4, 4))
#     afw_display = afwDisplay.Display(frame=fig)
#     afw_display.scale('asinh', 'zscale')
#     afw_display.mtv(cutout_image.image)
    
#     cutout_wcs = cutout_image.getWcs()
#     radec = geom.SpherePoint(ra, dec, geom.degrees)
#     xy = geom.PointI(cutout_wcs.skyToPixel(radec))
    
#     afw_display.dot('x', xy.getX(), xy.getY(), size=1, ctype='orange')
#     plt.gca().axis('off')
#     plt.savefig(out_name)
    
    fig = plt.figure()
    plt.subplot(projection=WCS(cutout_image.getWcs().getFitsMetadata()))
    calexp_extent = (cutout_image.getBBox().beginX, cutout_image.getBBox().endX,
                 cutout_image.getBBox().beginY, cutout_image.getBBox().endY)
    im = plt.imshow(cutout_image.image.array, cmap='gray', vmin=-200.0, vmax=1000,
                extent=calexp_extent, origin='lower')
    plt.colorbar(location='right', anchor=(0, 0.1))
    plt.gca().axis('off')
    # plt.xlabel('Right Ascension')
    # plt.ylabel('Declination')
    plt.savefig(out_name)
    
    return fig


def remove_figure(fig):
    """
    Remove a figure to reduce memory footprint.
    Parameters
    ----------
    fig: matplotlib.figure.Figure
        Figure to be removed.
    Returns
    -------
    None
    """
    # get the axes and clear their images
    for ax in fig.get_axes():
        for im in ax.get_images():
            im.remove()
    fig.clf()       # clear the figure
    plt.close(fig)  # close the figure

    gc.collect()    # call the garbage collector
    
def get_flux(flux_table):
    """
    Create dictionary of light curve.
    This should be follwed by plotlc

    Input Parameters
    ----------
    flux_table : from query_flux
    
    Returns
    ----------
    two dictionaries of days in MJD and flux for each band
    """
    pick = {}
    for filter in plot_filter_labels:
        pick[filter] = (flux_table['band'] == filter)
    mjd_days = {}
    mags = {}
    for filter in plot_filter_labels:
        mjd_days[filter] = np.array(flux_table[pick[filter]]['expMidptMJD']) * u.day
        mags[filter] = np.array(flux_table[pick[filter]]['psfMag'])
        
    return mjd_days, mags
    
def plotlc(bands, days, magnitudes, out_name):
    """
    Create a light curve.

    Input Parameters
    ----------
    days : dictionary for MJD in each band 
    magnitudes : dictionary for flux in each band
    out_name : file name where you'd like to save it
    
    Returns
    ----------
    light curve image
    """
    
    fig = plt.figure(figsize=(10,4))
    for band in bands:
        plt.plot(days[band], magnitudes[band],\
                 plot_filter_symbols[band], ms=4, label=plot_filter_labels[band])
    plt.minorticks_on()
    plt.xlabel('MJD (days)')
    plt.ylabel('magnitude')
    plt.legend('upper right')
    plt.legend()
    plt.savefig(out_name)
    return fig
    
    
def make_figure(exp, out_name):
    """
    Create an image.
    should be followed with remove_figure
    Parameters
    ----------
    exp : calexp from butler.get
    out_name : file name where you'd like to save it
    
    """
    fig = plt.figure(figsize=(10, 8))
    afw_display = afwDisplay.Display(1)
    afw_display.scale('asinh', 'zscale')
    afw_display.mtv(exp.image)
    plt.gca().axis('on')
    plt.savefig(out_name)
    
    return fig

def get_bandtractpatch(ra_deg, dec_deg, skymap):
    """
    get the tract and patch of a source. currently retrieves i band only. 
    Parameters
    ----------
    ra : ra of source in degrees
    dec : dec of source in degrees
    
    """
    spherePoint = lsst.geom.SpherePoint(ra_deg*lsst.geom.degrees, dec_deg*lsst.geom.degrees)
    tract = skymap.findTract(spherePoint)
    patch = tract.findPatch(spherePoint)
    my_tract = tract.tract_id
    my_patch = patch.getSequentialIndex()
    dataId = {'band': 'i', 'tract': my_tract, 'patch': my_patch}
    return dataId

def setup_plotting():

    # Set up some plotting defaults:       
    params = {'axes.labelsize': 20,
              'font.size': 20,
              'legend.fontsize': 14,
              'xtick.major.width': 3,
              'xtick.minor.width': 2,
              'xtick.major.size': 12,
              'xtick.minor.size': 6,
              'xtick.direction': 'in',
              'xtick.top': True,
              'lines.linewidth': 3,
              'axes.linewidth': 3,
              'axes.labelweight': 3,
              'axes.titleweight': 3,
              'ytick.major.width': 3,
              'ytick.minor.width': 2,
              'ytick.major.size': 12,
              'ytick.minor.size': 6,
              'ytick.direction': 'in',
              'ytick.right': True,
              'figure.figsize': [8, 8],
              'figure.facecolor': 'White'
              }

    plt.rcParams.update(params)

    #initializing Tap and Butler
    pandas.set_option('display.max_rows', 20)
    warnings.simplefilter("ignore", category=UnitsWarning)
    
    # Use lsst.afw.display with the matplotlib backend
    afwDisplay.setDefaultBackend('matplotlib')
    
def setup_butler(config, collection):    
    service = get_tap_service()
    assert service is not None
    assert service.baseurl == "https://data.lsst.cloud/api/tap"

    butler = dafButler.Butler(config, collections=collection)
    skymap = butler.get('skyMap')
    
    return service, butler, skymap
