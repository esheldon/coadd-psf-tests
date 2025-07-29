import numpy as np
import ngmix
import matplotlib.pyplot as mplt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import esutil as eu


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--flist', nargs='+', required=True)
    parser.add_argument('--diff', action='store_true')
    parser.add_argument('--bins', type=int, default=100)
    parser.add_argument('--output', required=True)
    parser.add_argument('--show', action='store_true')
    return parser.parse_args()


def main():
    args = get_args()

    data = eu.io.read(args.flist)

    fwhm = ngmix.moments.T_to_fwhm(data['T']) * 0.2
    if args.diff:
        fwhm_cen = ngmix.moments.T_to_fwhm(data['Tcen']) * 0.2
        fwhm = fwhm - fwhm_cen

    fwhm_stats = eu.stat.print_stats(fwhm)
    fm = f'{fwhm_stats["mean"]:.2g}'
    fs = f'{fwhm_stats["std"]:.2g}'
    ftext = (
        r'$\mu$ = ' + fm
        + r' $\sigma$ = ' + fs
    )

    counts, xedges, yedges = np.histogram2d(data['ra'], data['dec'], bins=args.bins)

    # Compute 2D histogram for sum of z-values
    fwhm_sum, _, _ = np.histogram2d(
        data['ra'], data['dec'], bins=args.bins, weights=fwhm,
    )

    # Compute mean: z_sum / counts (avoid division by zero)
    mean_fwhm = np.divide(
        fwhm_sum,
        counts,
        out=np.zeros_like(fwhm_sum),
        # where=counts != 0
    )

    # if args.diff:
    #     mean_fwhm = mean_fwhm.clip(min=-0.03, max=0.03)

    radiff = data['ra'].max() - data['ra'].min()
    decdiff = data['dec'].max() - data['dec'].min()

    if radiff > decdiff:
        width = 10
        height = width * decdiff / radiff
    else:
        height = 10
        width = height * radiff / decdiff

    fig, ax = mplt.subplots(figsize=(width, height))

    ax.set(xlabel='RA', ylabel='DEC')
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

    ax.text(
        0.1, 0.95,
        ftext,
        color='black',
        transform=ax.transAxes
    )
    ax.text(
        0.1+0.001, 0.95+0.001,
        ftext,
        color='white',
        transform=ax.transAxes
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if args.diff:
        name = r'$\Delta$ FWHM [arcsec]'
    else:
        name = 'FWHM [arcsec]'

    fig.colorbar(cim, cax=cax, label=name)

    if args.show:
        mplt.show()

    print('writing:', args.output)
    fig.savefig(args.output, dpi=150)


main()
