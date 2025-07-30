#!/usr/bin/env python


def calc_rho_stats(
    ra, dec,
    T,
    g1, g2,
    weights,
    Tcen,
    g1cen, g2cen,
    nproc=1,
    npatch=None,
    min_sep=0.5,
    max_sep=50,
    bin_size=0.2,
    file_name=None,
    model_properties=None,
    low_mem=False,
    **kwargs,
):
    import treecorr

    # treecorr.set_max_omp_threads(1)

    tckwargs = kwargs
    tckwargs['min_sep'] = min_sep
    tckwargs['max_sep'] = max_sep
    tckwargs['bin_size'] = bin_size
    # tckwargs['var_method'] = 'jackknife'
    tckwargs['num_threads'] = nproc
    # tckwargs['verbose'] = 2  # show some progress
    tckwargs['var_method'] = 'bootstrap'

    if 'sep_units' not in tckwargs:
        tckwargs['sep_units'] = 'arcmin'

    # Set this to true if there is a problem and we need to skip plots.
    # skip = False

    # get the shapes
    print(f'using {ra.size} stars')

    dT = T - Tcen
    dg1 = g1 - g1cen
    dg2 = g2 - g2cen

    # make the treecorr catalogs
    print("creating Treecorr Catalogs")

    print('    cat_g')
    cat_g = treecorr.Catalog(
        ra=ra, dec=dec,
        ra_units='deg', dec_units='deg',
        g1=g1, g2=g2,
        w=weights,
        npatch=npatch,
    )
    print('    cat_dg')
    cat_dg = treecorr.Catalog(
        ra=ra, dec=dec,
        ra_units='deg', dec_units='deg',
        g1=dg1, g2=dg2,
        w=weights,
        patch_centers=cat_g.patch_centers,
    )
    print('    cat_gdTT')
    cat_gdTT = treecorr.Catalog(
        ra=ra, dec=dec,
        ra_units='deg', dec_units='deg',
        g1=g1 * dT / T, g2=g2 * dT / T,
        w=weights,
        patch_centers=cat_g.patch_centers,
    )

    # setup and run the correlations
    print("doing rho stats")

    # save the rho objects
    data = {}
    print('    rho1')
    data['rho1'] = treecorr.GGCorrelation(tckwargs)
    data['rho1'].process(cat_dg, low_mem=low_mem)

    print('    rho2')
    data['rho2'] = treecorr.GGCorrelation(tckwargs)
    data['rho2'].process(cat_g, cat_dg, low_mem=low_mem)

    print('    rho3')
    data['rho3'] = treecorr.GGCorrelation(tckwargs)
    data['rho3'].process(cat_gdTT, low_mem=low_mem)

    print('    rho4')
    data['rho4'] = treecorr.GGCorrelation(tckwargs)
    data['rho4'].process(cat_dg, cat_gdTT, low_mem=low_mem)

    print('    rho5')
    data['rho5'] = treecorr.GGCorrelation(tckwargs)
    data['rho5'].process(cat_g, cat_gdTT, low_mem=low_mem)
    # treecorr.set_max_omp_threads(None)

    return data


def read_data(args):
    import esutil as eu
    import numpy as np

    columns = [
        'ra', 'dec',
        'T',
        'e1',
        'e2',
        'Tcen',
        'e1cen',
        'e2cen',
        'visit_count',
    ]

    data = eu.io.read(args.flist, columns=columns)

    if args.min_nepoch is not None:
        print('trimming to visit_count > ', args.min_nepoch)
        w, = np.where(data['visit_count'] > args.min_nepoch)
        print(f'kept {w.size}/{data.size}')
        data = data[w]

    return data


def extract_data(st):
    import ngmix

    g1, g2 = ngmix.shape.e1e2_to_g1g2(st['e1'], st['e2'])
    g1cen, g2cen = ngmix.shape.e1e2_to_g1g2(st['e1cen'], st['e2cen'])

    T = st['T']
    Tcen = st['Tcen']

    return g1, g2, T, g1cen, g2cen, Tcen


def main():
    args = get_args()

    print('reading data')
    st = read_data(args)

    g1, g2, T, g1cen, g2cen, Tcen = extract_data(st)

    data = calc_rho_stats(
        ra=st['ra'],
        dec=st['dec'],
        g1=g1,
        g2=g2,
        weights=None,
        T=T,
        g1cen=g1cen,
        g2cen=g2cen,
        Tcen=Tcen,
        npatch=args.npatch,
        min_sep=args.min_sep,
        max_sep=args.max_sep,
        nproc=args.nproc,
        low_mem=args.low_mem,
    )

    for key in data:
        fname = args.front + f'-{key}.fits'
        print(f'writing: {fname}')
        data[key].write(fname)


def get_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--front', required=True,
                        help='front for file names')

    parser.add_argument('flist', nargs='+')
    parser.add_argument('--min-nepoch', type=int)

    parser.add_argument(
        '--low-mem', action='store_true',
        help='use low_mem mode',
    )

    parser.add_argument(
        '--nstar-min',
        type=int, default=50,
        help='only use images with at least this many stars, default 50'
    )
    parser.add_argument(
        '--nproc', type=int, default=1,
        help='number of processes to use, default 1'
    )
    parser.add_argument(
        '--npatch',
        type=int,
        default=40,
        help='number of patches for bootstrap, default 40'
    )

    parser.add_argument(
        '--min-sep', type=float, default=0.5,
        help='minimum separation, default 0.5 arcmin'
    )
    parser.add_argument(
        '--max-sep', type=float, default=50,
        help='maximum separation, default 50 arcmin'
    )

    parser.add_argument(
        '--seeing-min', type=float, default=0,
        help='min seeing, default 0'
    )

    return parser.parse_args()


main()
