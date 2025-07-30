def plot_rho_stats(args, data, logger=None, **kwargs):
    import numpy as np

    from matplotlib.figure import Figure
    fig = Figure(figsize=(10.5, 5))
    # In matplotlib 2.0, this will be
    # axs = fig.subplots(ncols=2)
    axs = [fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)]
    axs = np.array(axs, dtype=object)

    axs[0].set_xlim(0.5 * args.min_sep, 1.5 * args.max_sep)
    axs[1].set_xlim(0.5 * args.min_sep, 1.5 * args.max_sep)

    gcolor = 'lightgrey'
    for ax in axs:
        ax.grid(which='minor', color=gcolor, linestyle=':', linewidth=0.5)
        ax.grid(True, color=gcolor)

    axs[0].set_xlabel(r'$\theta$ (arcmin)')
    axs[0].set_ylabel(r'$\rho(\theta)$')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log', nonpositive='clip')

    axs[1].set_xlabel(r'$\theta$ (arcmin)')
    # axs[1].set_ylabel(r'$\rho(\theta)$')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log', nonpositive='clip')

    # Left plot is rho1,3,4
    rho1 = _plot_single(axs[0], data['rho1'], 'blue', 'o')
    rho3 = _plot_single(axs[0], data['rho3'], 'green', 's', 0.1)
    rho4 = _plot_single(axs[0], data['rho4'], 'red', '^', 0.2)

    axs[0].legend([rho1, rho3, rho4],
                  [r'$\rho_1(\theta)$',
                   r'$\rho_3(\theta)$',
                   r'$\rho_4(\theta)$'],
                  loc='upper right', fontsize=12)

    # Right plot is rho2,5
    rho2 = _plot_single(axs[1], data['rho2'], 'blue', 'o')
    rho5 = _plot_single(axs[1], data['rho5'], 'green', 's', 0.1)

    axs[1].legend([rho2, rho5],
                  [r'$\rho_2(\theta)$', r'$\rho_5(\theta)$'],
                  loc='upper right', fontsize=12)

    axs[0].set_ylim(args.ymin, args.ymax)
    axs[1].set_ylim(args.ymin, args.ymax)
    fig.tight_layout()
    return fig, axs


def _plot_single(ax, rho, color, marker, offset=0.):
    import numpy as np

    # Add a single rho stat to the plot.
    meanr = rho.meanr * (1. + rho.bin_size * offset)
    xip = rho.xip
    sig = np.sqrt(rho.varxip)
    ax.plot(meanr, xip, color=color)
    ax.plot(meanr, -xip, color=color, ls=':')
    ax.errorbar(
        meanr[xip > 0], xip[xip > 0], yerr=sig[xip > 0], color=color, ls='',
        marker=marker
    )
    ax.errorbar(
        meanr[xip < 0], -xip[xip < 0], yerr=sig[xip < 0], color=color, ls='',
        marker=marker, fillstyle='none', mfc='white'
    )
    return ax.errorbar(-meanr, xip, yerr=sig, color=color, marker=marker)


def get_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--front', required=True,
                        help='front for file names')

    parser.add_argument('--plot-file', required=True)

    parser.add_argument(
        '--min-sep', type=float, default=0.5,
        help='minimum separation, default 0.5 arcmin'
    )
    parser.add_argument(
        '--max-sep', type=float, default=50,
        help='maximum separation, default 50 arcmin'
    )

    parser.add_argument(
        '--ymin', type=float, default=1.0e-12,
        help='ymin for plot, default 1.0e-10'
    )
    parser.add_argument(
        '--ymax', type=float, default=1.e-6,
        help='ymax for plot, default 1,0e-6'
    )

    return parser.parse_args()


def read_rho_stats(front):
    import treecorr

    data = {}
    for key in ['rho1', 'rho2', 'rho3', 'rho4', 'rho5']:
        fname = front + f'-{key}.fits'
        print('reading:', fname)
        data[key] = treecorr.GGCorrelation.from_file(fname)

    return data


def main():
    args = get_args()

    data = read_rho_stats(args.front)

    fig, ax = plot_rho_stats(args, data, logger=None)
    print('writing:', args.plot_file)
    fig.savefig(args.plot_file, dpi=150)


main()
