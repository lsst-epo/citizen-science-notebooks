# import packages used for generating subject set
from astropy.units import UnitsWarning
import matplotlib.pyplot as plt
import gc
import os, uuid, itertools
import pandas
import warnings

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

    gc.collect()  # call the garbage collector


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

    # config = 'dp02'
    # collection = '2.2i/runs/DP0.2'
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
    if os.path.isdir(batch_dir) == False:
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
            "sourceId": row.objectId,  # required column, do not change the column name
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
    if os.path.isdir(batch_dir) == False:
        os.mkdir(batch_dir)

    # Get field names
    col_names = list(results_table.fieldnames)

    # Loop over results_table, or any other iterable provided by the PI:
    for row in results_table:

        csv_row = { "sourceId": str(uuid.uuid4()) }

        for col in col_names:

            csv_row[col] = row[col]

        manifest_dict.append(csv_row)

    return manifest_dict
