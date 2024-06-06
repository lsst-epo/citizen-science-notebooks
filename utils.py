# import packages used for generating subject set
from astropy.units import UnitsWarning
import matplotlib.pyplot as plt
import gc
import os
import warnings
import pandas as pd
import time
import os

# Astropy imports
from astropy.wcs import WCS
from astropy.visualization import make_lupton_rgb
from astropy import units as u
from astropy.coordinates import SkyCoord

# Import the Rubin TAP service utilities
from lsst.rsp import get_tap_service

# Image visualization routines.
import lsst.afw.display as afwdisplay

# The Butler provides programmatic access to LSST data products.
import lsst.daf.butler as dafbutler
import lsst.geom

# Must explicitly set this to save figures
afwdisplay.setDefaultBackend("matplotlib")

plt.style.use("tableau-colorblind10")

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
    fig = plt.figure(figsize=(4, 4))
    afw_display = afwDisplay.Display(frame=fig)
    afw_display.scale('asinh', 'zscale')
    afw_display.mtv(cutout_image.image)
    
    cutout_wcs = cutout_image.getWcs()
    radec = geom.SpherePoint(ra, dec, geom.degrees)
    xy = geom.PointI(cutout_wcs.skyToPixel(radec))
    
    afw_display.dot('x', xy.getX(), xy.getY(), size=1, ctype='orange')
    plt.gca().axis('off')
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
    fig.clf()  # clear the figure
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
    afw_display = afwdisplay.Display(1)
    afw_display.scale("asinh", "zscale")
    afw_display.mtv(exp.image)
    plt.gca().axis("on")
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
    spherepoint = lsst.geom.SpherePoint(
        ra_deg * lsst.geom.degrees, dec_deg * lsst.geom.degrees
    )
    tract = skymap.findTract(spherepoint)
    patch = tract.findPatch(spherepoint)
    my_tract = tract.tract_id
    my_patch = patch.getSequentialIndex()
    dataid = {"band": "i", "tract": my_tract, "patch": my_patch}
    return dataid


def setup_plotting():
    # Set up some plotting defaults:
    params = {
        "axes.labelsize": 20,
        "font.size": 20,
        "legend.fontsize": 14,
        "xtick.major.width": 3,
        "xtick.minor.width": 2,
        "xtick.major.size": 12,
        "xtick.minor.size": 6,
        "xtick.direction": "in",
        "xtick.top": True,
        "lines.linewidth": 3,
        "axes.linewidth": 3,
        "axes.labelweight": 3,
        "axes.titleweight": 3,
        "ytick.major.width": 3,
        "ytick.minor.width": 2,
        "ytick.major.size": 12,
        "ytick.minor.size": 6,
        "ytick.direction": "in",
        "ytick.right": True,
        "figure.figsize": [8, 8],
        "figure.facecolor": "White",
    }

    plt.rcParams.update(params)

    # initializing Tap and Butler
    pandas.set_option("display.max_rows", 20)
    warnings.simplefilter("ignore", category=UnitsWarning)

    # Use lsst.afw.display with the matplotlib backend
    afwdisplay.setDefaultBackend("matplotlib")


def setup_butler(config, collection):
    service = get_tap_service()
    assert service is not None
    assert service.baseurl == "https://data.lsst.cloud/api/tap"
    butler = dafbutler.Butler(config, collections=collection)
    skymap = butler.get("skyMap")
    return service, butler, skymap


def run_butler_query(service, number_sources, use_center_coords, use_radius):
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


def prep_table(results, skymap):
    results_table = results.to_table().to_pandas()
    results_table["dataId"] = results_table.apply(
        lambda x: get_bandtractpatch(x["coord_ra"], x["coord_dec"], skymap),
        axis=1
    )
    return results_table


def make_manifest_with_images(results_table, butler, batch_dir):
    # In-memory manifest file as an array of dicts
    manifest = []

    # Create directory if it does not already exist
    if os.path.isdir(batch_dir) is False:
        os.mkdir(batch_dir)

    # Loop over results_table, or any other iterable provided by the PI:
    for index, row in results_table.iterrows():
        # Use the Butler to get data for each index, row
        deepCoadd = butler.get("deepCoadd", dataId=row["dataId"])
        filename = "cutout" + str(row["objectId"]) + ".png"
        figout = make_figure(deepCoadd, batch_dir + filename)

        # Create the CSV-file-row-as-dict
        csv_row = {
            "filename": filename,  # required column, do not change the column name
            "objectId": row.objectId,  # required column, do not change the column name
            "objectIdType": "DIRECT",  # required column, do not change the column name
            # Add your desired columns:
            "coord_ra": row.coord_ra,
            "coord_dec": row.coord_dec,
            "g_cModelFlux": row.g_cModelFlux,
            "r_cModelFlux": row.r_cModelFlux,
            "r_extendedness": row.r_extendedness,
            "r_inputCount": row.r_inputCount,
        }
        manifest.append(csv_row)
        remove_figure(figout)

    return manifest


def make_manifest_with_tabular_data(results_table, batch_dir):
    # In-memory manifest file as an array of dicts
    manifest_dict = []

    # Create directory if it does not already exist
    if os.path.isdir(batch_dir) is False:
        os.mkdir(batch_dir)

    # Get field names
    col_names = list(results_table.fieldnames)

    # Loop over results_table, or any other iterable provided by the PI:
    for row in results_table:

        # csv_row = { "sourceId": str(uuid.uuid4()) }
        csv_row = {}

        for col in col_names:
            if col == "objectId":
                csv_row["sourceId"] = row[col]
            else:
                csv_row[col] = row[col]

        manifest_dict.append(csv_row)

    return manifest_dict
