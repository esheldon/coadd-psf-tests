import numpy as np
import matplotlib.pyplot as mplt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import esutil as eu


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--flist', nargs='+', required=True)
    parser.add_argument('--bins', type=int, default=100)
    parser.add_argument('--output', required=True)
    parser.add_argument('--show', action='store_true')
    return parser.parse_args()


def main():
    args = get_args()

    if len(args.flist) == 1:
        print('reading:', args.flist)

    data = eu.io.read(args.flist)

    nepoch = data['visit_count']

    nepoch_stats = eu.stat.print_stats(nepoch)
    fm = f'{nepoch_stats["mean"]:.2g}'
    fs = f'{nepoch_stats["std"]:.2g}'
    ftext = (
        r'$\mu$ = ' + fm
        + r' $\sigma$ = ' + fs
    )

    counts, xedges, yedges = np.histogram2d(data['ra'], data['dec'], bins=args.bins)

    # Compute 2D histogram for sum of z-values
    nepoch_sum, _, _ = np.histogram2d(
        data['ra'], data['dec'], bins=args.bins, weights=nepoch,
    )

    # Compute mean: z_sum / counts (avoid division by zero)
    mean_nepoch = np.divide(
        nepoch_sum,
        counts,
        out=np.zeros_like(nepoch_sum),
        # where=counts != 0
    )

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
        mean_nepoch.T,
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
    name = r'N$_{\mathrm{epoch}}$'

    fig.colorbar(cim, cax=cax, label=name)

    fig.tight_layout()
    if args.show:
        mplt.show()

    print('writing:', args.output)
    fig.savefig(args.output, dpi=150)


main()
