
import argparse
from compute_cn_pcf import compute_cn_pcf
def arg_parse(args=None):
    """
    Parse command line arguments
    Returns:
    """
    description = ''
    parser = argparse.ArgumentParser(
        prog='hatchet_pcf compute-cn',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        '-i',
        '--bbfile',
        type=str,
        required=True,
        help='Prefix path to bb file (required)',
    )

    parser.add_argument(
        '-m',
        '--minK',
        type=int,
        required=False,
        default=None,
        help='Minimum number of cluster to infer. Default (None: Number of chromosomal arms analyzed)',
    )

    parser.add_argument(
        '-M',
        '--maxK',
        type=int,
        required=False,
        default=None,
        help='Maximum number of clones. Default (None: Number of genes used as input)',
    )
    
    parser.add_argument(
        '-v',
        '--covar',
        type=str,
        required=False,
        default='diag',
        help='Form of covariance matrix: spherical, diag, full, or tied (default: diag)',
    )

    parser.add_argument(
        '-z',
        '--decode_alg',
        type=str,
        required=False,
        default='map',
        help='Decoding algorithm to use: map or viterbi (default: map)'
    )

    parser.add_argument(
        '-t',
        '--tmat',
        type=str,
        required=False,
        default='diag',
        help='Form of transition matrix to infer: fixed, diag (1-parameter), or full (default: diag)',
    )


    parser.add_argument(
        '-T',
        '--tau',
        type=float,
        required=False,
        default=0.000001,
        help='Off-diagonal value for initializing transition matrix (Default: 0.000001)'
    )

    parser.add_argument(
        '-r',
        '--restarts',
        type=int,
        required=False,
        default=10,
        help="Number of restarts performed by the clustering to choose the best (default: 10)",
    )


    parser.add_argument(
        '-y',
        '--subset',
        required=False,
        default=None,
        type=str,
        nargs='+',
        help=(
            'List of sample names to use as a subset of those included in binning' '(default: none, run on all samples)'
        ),
    )


    parser.add_argument(
        '-p',
        '--preserve_arm',
        required=False,
        default=False,
        type=bool,
        help='Preserve chromosome arms in local clustering (default: none)'
    )


    parser.add_argument(
        '-x',
        '--runningdir',
        type=str,
        required=False,
        default="./",
        help='Running directory (default: ./)',
    )
   
    parser.add_argument(
        '-f',
        '--ampdel',
        action='store_true',
        default=False,
        required=False,
        help=
            'Remove amp-del assumption where each mutated allele of every segment can be either amplified or '
            'deleted in all tumor clones w.r.t. base (2 for diploid and 4 for tetraploid) (default: use assumption)'
        ,
    )

    parser.add_argument(
        '-d',
        '--cnstates',
        type=int,
        required=False,
        default=-1,
        help='Maximum number of distinct copy-number states for each segment (default: None, no limit)',
    )

    parser.add_argument(
        '-tS',
        type=float,
        required=False,
        default=0.75,
        help='The minimum proportion of samples within clonal cluster to call cluster as clonal (default: 0.75)',
    )

    parser.add_argument(
        '-tR',
        type=float,
        required=False,
        default=0.08,
        help='RDR tolerance for clonal clusters (default: 0.08)',
    )

    parser.add_argument(
        '-tB',
        type=float,
        required=False,
        default=0.04,
        help='BAF tolerance for clonal clusters (default: 0.04)',
    )

    parser.add_argument(
        '-tD',
        type=float,
        required=False,
        default=0.1,
        help=(
            'Maximum BAF shift for neutral cluster used to automatically infer the diploid/tetraploid cluster '
            '(default: 0.1)'
        ),
    )

    parser.add_argument(
        '-minC',
        type=int,
        required=False,
        default=2,
        help='Minimum number of clines. default: 2',
    )

    parser.add_argument(
        '-maxC',
        type=int,
        required=False,
        default=8,
        help='Maximum number of clones. default: 8',
    )


    parser.add_argument(
        '-solver',
        type=str,
        required=False,
        default="cpp",
        help='Solver type to use. default cpp',
    )

    parser.add_argument(
        '-solvermode',
        type=str,
        required=False,
        default="cd",
        help='Solver type to use. default cd',
    )

    parser.add_argument(
        '-cnmax',
        type=int,
        required=False,
        default=None,
        help='Maximum number of copy number. default None',
    )


    parser.add_argument(
        '-l',
        '--limitinc',
        type=float,
        required=False,
        default=0.6,
        help=(
            'Upper bound to the relative increase of objective function. When there are significant small CNAs, '
            'their effect on the objective function may be confounded by only larger events, use this value to '
            'limit the relative increase of OBJ so that fitting small CNAs is more considered (default: None)'
        ),
    )

    parser.add_argument(
        '-g',
        '--ghostprop',
        type=float,
        required=False,
        default=0.3,
        help=(
            'Increasing proportion used to compute the value of the first ghost point added in the solution '
            'selection (default: 0.3)'
        ),
    )

    parser.add_argument(
        '-u',
        '--mu',
        type=float,
        required=False,
        default=0.03,
        help='Minimum clone proporion in each sample (default: 0.03)',
    )

    parser.add_argument(
        '--maxiterations',
        type=int,
        required=False,
        default=10,
        help='Maximum number of iterations composed of C-step/U-step for each seed (default: 10)',
    )

    parser.add_argument(
        '-s',
        '--seeds',
        type=int,
        required=False,
        default=400,
        help='Number of seeds for coordinate-descent method (default: 400)',
    )

    parser.add_argument(
        '-j',
        '--jobs',
        type=int,
        required=False,
        default=1,
        help='Number of parallel jobs (default: 1)',
    )

    parser.add_argument(
        '-rs',
        '--randomseed',
        type=int,
        required=False,
        default=None,
        help='Random seed (default: None)',
    )
    parser.add_argument(
        '-tl',
        '--timelimit',
        type=int,
        required=False,
        default=None,
        help='Time limit for each ILP run (default: None)',
    )
    parser.add_argument(
        '-ml',
        '--memlimit',
        type=int,
        required=False,
        default=None,
        help='Memory limit for each ILP run (default: None)',
    )

    args = parser.parse_args(args)
    return args

def main():
    args=vars(arg_parse())
    compute_cn_pcf(
        bbfile=args["bbfile"],
        minK=args["minK"],
        maxK=args["maxK"],
        covar=args["covar"],
        decode_alg=args["decode_alg"],
        tmat=args["tmat"],
        tau=args["tau"],
        restarts=args["restarts"],
        preserve_arm=args["preserve_arm"],
        subset=args["subset"],
        tD=args["tD"],
        tB=args["tB"],
        tR=args["tR"],
        tS=args["tS"],
        minC=args["minC"], 
        maxC=args["maxC"], 
        solver=args["solver"],
        solve_mode=args["solvermode"],
        d=args["cnstates"],
        cn_max=args["cnmax"],
        mu=args["mu"],
        dT=args["tD"],
        g=args["ghostprop"],
        limit=args["limitnc"],
        ampdel=args["ampdel"],
        n_seed=args["seeds"],
        n_worker=args["jobs"],
        random_seed=args["randomseed"],
        timelimit=args["timelimit"],
        rundir=args["runningdir"]
    )
if __name__ == '__main__':
    main()