"""
Based on this notebook
https://github.com/lsst-sitcom/sciunit_wlshear/blob/main/notebooks/access-psf-coadds.ipynb
"""
from lsst.daf.butler import Butler
import numpy as np
import gc

from lsst.afw.geom.ellipses import Quadrupole, SeparableDistortionTraceRadius
from lsst import geom
import lsst.pex.exceptions
import fitsio

# REPO = '/repo/main'
REPO = '/repo/embargo'

# example collections
#    A360 (r and i): u/mgorsuch/a360_cell_coadd/20250513T044026Z
#    A360 (g): u/mgorsuch/a360_cell_coadd_g/20250611T144112Z
# those have tracts
#    [10463, 10464, 10704, 10705]

# LSSTCam
# no cells yet
# LSSTCam/runs/DRP/20250501_20250604/w_2025_23/DM-51284
# has cells
# LSSTCam/runs/DRP/20250501_20250609/w_2025_26/DM-51580


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--nrand', type=int)
    parser.add_argument('--range', nargs=2, type=int, help='range as a slice')
    parser.add_argument('--band', default='i')
    parser.add_argument('--tracts', type=int, nargs='+')
    return parser.parse_args()


def make_field_info(tract, patch):
    st = np.zeros(1, dtype=[('tract', 'i4'), ('patch', 'i4')])
    st['tract'] = tract
    st['patch'] = patch
    return st


def get_field_info(butler, collections, data_kwargs, band, tracts=None):
    """
    Retrieves the unique tract/patch combinations within a specified collection.


    Parameters
    ----------
    butler: Butler object
        The butler for queries
    collection:
        the relevant collection containing cell-based coadds of interest
    data_kwargs: dict
        dictionary of specific instrument and skymap used for butler query
    band: str
        Band string

    Returns
    --------
    field_data: array
        with columns for available tract and patch IDs within collection
    """
    import esutil as eu
    if tracts is not None:
        print('will limit to tracts:', tracts)

    dlist = []
    for ref in butler.registry.queryDatasets(
        # 'deepCoaddCell',
        'deep_coadd_cell_predetection',
        band=band,
        collections=collections,
        instrument=data_kwargs['instrument'],
        skymap=data_kwargs['skymap'],
    ):

        st = make_field_info(
            ref.dataId.get('tract'),
            ref.dataId.get('patch'),
        )
        if tracts is not None and st['tract'][0] not in tracts:
            continue

        dlist.append(st)

    return eu.numpy_util.combine_arrlist(dlist)


def get_cell_count(field_data, butler, collections, data_kwargs, band):
    """
    Retrieve the total number of cells from your input collection

    -- Inputs --

    field_data: array
        containing tract/patch combinations of interest
    butler: Butler object
    collection:
        the relevant collection containing cell-based coadds of interest
    data_kwargs: dictionary
        of specific instrument and skymap used for butler query

    -- Returns --

    cell_count: number of cells with inputs in specified field data

    NOTE: DOES include duplicate cells due to overlap of patches/tracts
    """
    from esutil.pbar import pbar

    cell_count = 0

    for field in pbar(field_data):

        coadd = butler.get(
            # 'deepCoaddCell',
            'deep_coadd_cell_predetection',
            collections=collections,
            instrument=data_kwargs['instrument'],
            skymap=data_kwargs['skymap'],
            tract=field['tract'],
            patch=field['patch'],
            band=band,
        )

        cells = len(list(coadd.cells.keys()))  # get number of non-empty cells
        cell_count += cells
        del coadd
        gc.collect()

    return cell_count


def make_data(n=1, include_cen=False):
    dtype = [
        ('tract', 'i4'),
        ('patch', 'i4'),
        ('visit_count', 'i2'),
        ('x_index', 'i4'),
        ('y_index', 'i4'),
        ('ra', 'f8'),
        ('dec', 'f8'),
        ('T', 'f4'),
        ('e1', 'f4'),
        ('e2', 'f4'),
    ]
    if include_cen:
        dtype += [
            ('Tcen', 'f4'),
            ('e1cen', 'f4'),
            ('e2cen', 'f4'),
        ]
    return np.zeros(n, dtype=dtype)


def get_coadd_psf_e1e2T(coadd_psf, pos):
    shape = coadd_psf.computeShape(pos)

    i_xx, i_yy, i_xy = shape.getIxx(), shape.getIyy(), shape.getIxy()

    q = Quadrupole(i_xx, i_yy, i_xy)
    s = SeparableDistortionTraceRadius(q)

    e1, e2 = s.getE1(), s.getE2()
    T = 2 * shape.getTraceRadius() ** 2
    return e1, e2, T


def get_rand_bbox_xy(rng, bbox):
    xl = bbox.getMinX()
    xh = bbox.getMaxX()
    yl = bbox.getMinY()
    yh = bbox.getMaxY()

    pos = geom.Point2D(
        x=rng.uniform(low=xl, high=xh),
        y=rng.uniform(low=yl, high=yh),
    )
    return pos


def get_cell_data(
    field_data, butler, collections, data_kwargs, band, nrand=None, rng=None,
):
    """
    Iterates through cells in each patch to collect cell PSF infromation

    -- Inputs --

    field_data: array
        containing tract/patch combinations of interest
    butler: Butler object
    collection:
        the relevant collection containing cell-based coadds of interest
    data_kwargs: dictionary
        of specific instrument and skymap used for butler query

    -- Returns --

    data: array containing PSF information for each cell.
    """
    import esutil as eu
    from esutil.pbar import pbar

    print('getting cell count')
    cell_count = get_cell_count(
        field_data, butler, collections, data_kwargs, band=band,
    )

    print('cell num: ', cell_count)

    radec_set = set()

    dlist = []

    print('getting psf stats')
    for i, field in enumerate(pbar(field_data)):
        coadd_cell_patch = butler.get(
            # 'deepCoaddCell',
            'deep_coadd_cell_predetection',
            collections=collections,
            instrument=data_kwargs['instrument'],
            skymap=data_kwargs['skymap'],
            tract=field['tract'],
            patch=field['patch'],
            band=band,
        )
        coadd_psf = butler.get(
            # 'deepCoadd.psf',
            'deep_coadd.psf',
            collections=collections,
            instrument=data_kwargs['instrument'],
            skymap=data_kwargs['skymap'],
            tract=field['tract'],
            patch=field['patch'],
            band=band,
        )
        coadd = butler.get(
            # 'deepCoadd',
            'deep_coadd',
            collections=collections,
            instrument=data_kwargs['instrument'],
            skymap=data_kwargs['skymap'],
            tract=field['tract'],
            patch=field['patch'],
            band=band,
        )
        visits = coadd.getInfo().getCoaddInputs().visits

        wcs = coadd_cell_patch.wcs

        # skips empty cell indices
        cell_list = list(coadd_cell_patch.cells.keys())

        for cell_index in cell_list:

            cell = coadd_cell_patch.cells[cell_index]

            # collect cell center location
            # primarily used for removing duplicates due to patch overlap
            cell_bbox = cell.inner.bbox
            cell_center = cell_bbox.getCenter()
            cell_center_coord = wcs.pixelToSky(cell_center)

            racen = cell_center_coord[0].asDegrees()
            deccen = cell_center_coord[1].asDegrees()
            radec_tuple = (racen, deccen)

            if radec_tuple in radec_set:
                continue
            else:
                radec_set.add(radec_tuple)

            e1cen, e2cen, Tcen = get_coadd_psf_e1e2T(coadd_psf, cell_center)

            if nrand is not None:
                assert rng is not None
                for irand in range(nrand):
                    pos = get_rand_bbox_xy(rng=rng, bbox=cell_bbox)

                    coord = wcs.pixelToSky(pos)
                    try:
                        e1, e2, T = get_coadd_psf_e1e2T(coadd_psf, pos)

                        idata = make_data(include_cen=True)

                        idata['tract'] = field['tract']
                        idata['patch'] = field['patch']
                        idata['visit_count'] = len(
                            visits.subsetContaining(coord)
                        )
                        idata['x_index'] = cell_index.x
                        idata['y_index'] = cell_index.y

                        idata['ra'] = coord[0].asDegrees()
                        idata['dec'] = coord[1].asDegrees()
                        idata['T'] = T
                        idata['e1'] = e1
                        idata['e2'] = e2

                        idata['Tcen'] = Tcen
                        idata['e1cen'] = e1cen
                        idata['e2cen'] = e2cen

                        dlist.append(idata)
                    except lsst.pex.exceptions.Exception:
                        pass
            else:
                idata = make_data()
                idata['tract'] = field['tract']
                idata['patch'] = field['patch']
                idata['visit_count'] = len(
                    visits.subsetContaining(cell_center_coord)
                )
                idata['x_index'] = cell_index.x
                idata['y_index'] = cell_index.y
                idata['ra'] = racen
                idata['dec'] = deccen
                idata['T'] = Tcen
                idata['e1'] = e1cen
                idata['e2'] = e2cen

                dlist.append(idata)

        coadd = 0
        coadd_cell_patch = 0
        if (i % 5) == 0:
            gc.collect()

    gc.collect()

    data = eu.numpy_util.combine_arrlist(dlist)
    return data


def main():
    args = get_args()

    rng = np.random.RandomState()

    comcam_data_id = {
        # 'instrument': 'LSSTComCam',
        'instrument': 'LSSTCam',
        'skymap': 'lsst_cells_v1',
    }

    # cell_collection = 'u/mgorsuch/ComCam_Cells/Rubin_SV_38_7/20250214T210230Z'
    # cell_collection = 'u/mgorsuch/a360_cell_coadd/20250513T044026Z'
    cell_collections = [
        # not sure if this is needed
        'LSSTCam/runs/DRP/20250501_20250609/w_2025_26/DM-51580/20250703T162051Z',
        'LSSTCam/runs/DRP/20250501_20250609/w_2025_26/DM-51580',
        # 'u/mgorsuch/a360_cell_coadd/20250513T044026Z',
        # 'u/mgorsuch/a360_cell_coadd',
    ]
    cell_butler = Butler(REPO, collections=cell_collections)

    field_info = get_field_info(
        cell_butler, cell_collections, comcam_data_id, band=args.band,
        tracts=args.tracts,
    )

    if args.range is not None:
        field_info = field_info[args.range[0]:args.range[1]]

    data = get_cell_data(
        field_info, cell_butler, cell_collections, comcam_data_id,
        band=args.band,
        rng=rng,
        nrand=args.nrand,
    )

    print('writing:', args.output)
    fitsio.write(args.output, data, clobber=True)


main()
