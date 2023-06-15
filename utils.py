#import packages used for generating subject set

import matplotlib.pyplot as plt
import gc
import numpy as np
import pandas
import time

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

    

    #config = 'dp02'
    #collection = '2.2i/runs/DP0.2'
    butler = dafButler.Butler(config, collections=collection)
    skymap = butler.get('skyMap')
    
    return service, butler, skymap

def run_query(number_sources, use_center_coords, use_radius):
    query = (
        "SELECT TOP "
        + str(number_sources)
        + " "
        + "objectId, coord_ra, coord_dec, detect_isPrimary "
        + "g_cModelFlux, r_cModelFlux, r_extendedness, r_inputCount "
        + "FROM dp02_dc2_catalogs.Object "
        + "WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec), "
        + "CIRCLE('ICRS', "
        + use_center_coords
        + ", "
        + use_radius
        + ")) = 1 "
        + "AND detect_isPrimary = 1 "
        + "AND r_extendedness = 1 "
        + "AND scisql_nanojanskyToAbMag(r_cModelFlux) < 18.0 "
        + "ORDER by r_cModelFlux DESC"
    )
    results = service.search(query)
    assert len(results) == number_sources
    return results
