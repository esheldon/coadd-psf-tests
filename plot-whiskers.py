import numpy as np
import matplotlib.pyplot as mplt
import esutil as eu
from glob import glob


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--diff', action='store_true')
    parser.add_argument('--flist', nargs='+', required=True)
    parser.add_argument('--frac', type=float)
    return parser.parse_args()


def main():
    args = get_args()

    # if args.diff:
    #     flist = glob('psf-data-coadd-psf-nrand10-*.fits.gz')
    #
    # else:
    #     flist = 'psf-data-coadd-psf.fits.gz'

    data = eu.io.read(args.flist)

    fig, ax = mplt.subplots()
    ax.set(xlabel='RA', ylabel='DEC')

    if args.frac is not None:
        rng = np.random.RandomState()
        rind = rng.choice(data.size, size=int(args.frac * data.size))
        data = data[rind]

    if args.diff:
        e1 = data['e1'] - data['e1cen']
        e2 = data['e2'] - data['e2cen']
    else:
        e1 = data['e1']
        e2 = data['e2']

    u, v = eu.plotting.polar2whisker(
        e1=e1,
        e2=e2,
    )

    print('ra stats')
    eu.stat.print_stats(data['ra'])
    print('dec stats')
    eu.stat.print_stats(data['dec'])

    print('e1 stats')
    e1stats = eu.stat.print_stats(e1)
    print('e2 stats')
    e2stats = eu.stat.print_stats(e2)
    print('u stats')
    eu.stat.print_stats(u)
    print('v stats')
    eu.stat.print_stats(v)

    # return

    scale = 0.2
    eu.plotting.mwhiskers(
        plt=ax,
        xin=data['ra'],
        yin=data['dec'],
        uin=u,
        vin=v,
        scale=scale,
        color='black',
    )

    e1m = f'{e1stats["mean"]:.2g}'
    e1s = f'{e1stats["std"]:.2g}'
    e1text = (
        r'$\langle$ e$_1$ $\rangle$ = ' + e1m
        + r' $\sigma$ = ' + e1s
    )
    e2m = f'{e2stats["mean"]:.2g}'
    e2s = f'{e2stats["std"]:.2g}'
    e2text = (
        r'$\langle$ e$_2$ $\rangle$ = ' + e2m
        + r' $\sigma$ = ' + e2s
    )

    ax.text(
        0.1, 0.95,
        e1text,
        color='blue',
        transform=ax.transAxes
    )
    ax.text(
        0.65, 0.95,
        e2text,
        color='blue',
        transform=ax.transAxes
    )

    # val = 0.02
    # scaled_val = val * scale
    #
    # # transformed_val = scaled_val *
    # ax.plot(
    #     [0.1, 0.15],
    #     [0.9, 0.9],
    #     color='red',
    #     transform=ax.transAxes,
    # )
    # ax.text(
    #     0.1, 0.9,
    #     'hello',
    #     transform=ax.transAxes
    # )

    # fig.savefig('radec.png')
    mplt.show()


main()
