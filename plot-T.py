import numpy as np
import ngmix
import matplotlib.pyplot as mplt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import esutil as eu
from glob import glob


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--diff', action='store_true')
    return parser.parse_args()


def main():
    args = get_args()

    if args.diff:
        flist = glob('psf-data-coadd-psf-nrand10-*.fits.gz')

    else:
        flist = 'psf-data-coadd-psf.fits.gz'

    data = eu.io.read(flist)

    fig, ax = mplt.subplots()
    ax.set(xlabel='RA', ylabel='DEC')

    fwhm = ngmix.moments.T_to_fwhm(data['T']) * 0.2
    if args.diff:
        fwhm_cen = ngmix.moments.T_to_fwhm(data['Tcen']) * 0.2
        fwhm = fwhm - fwhm_cen

    bins = 1000
    counts, xedges, yedges = np.histogram2d(data['ra'], data['dec'], bins=bins)

    # Compute 2D histogram for sum of z-values
    fwhm_sum, _, _ = np.histogram2d(
        data['ra'], data['dec'], bins=bins, weights=fwhm,
    )

    # Compute mean: z_sum / counts (avoid division by zero)
    mean_fwhm = np.divide(
        fwhm_sum,
        counts,
        out=np.zeros_like(fwhm_sum),
        where=counts != 0
    )

    mean_fwhm = mean_fwhm.clip(min=-0.03, max=0.03)

    # Plot the result
    cim = ax.imshow(
        mean_fwhm.T,
        origin='lower',
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        # cmap='viridis',
        cmap='inferno',
        # cmap='gray',
        interpolation='nearest',
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if args.diff:
        name = r'$\Delta$ FWHM [arcsec]'
    else:
        name = 'FWHM [arcsec]'

    fig.colorbar(cim, cax=cax, label=name)

    mplt.show()


main()
