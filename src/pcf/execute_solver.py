import pandas as pd
import os
from solve_pcf import solve_pcf

def execute_solver(
    outprefix,
    n,
    df,
    clonal=None,
    solver='gurobi',
    solve_mode='cd',
    d=-1,
    cn_max=None,
    mu=0.1,
    dT=0.1,
    ampdel=FALSE,
    n_seed=10,
    n_worker=6,
    random_seed=None,
    max_iters=None,
    timelimit=None
):

    bbc_out_file = outprefix + '.bbc.ucn.tsv'


    obj, cA, cB, u, cluster_ids, sample_ids = solve_pcf(
        solver=solver,
        clonal=clonal,
        df=df,
        n=n,
        solve_mode=solve_mode,
        d=d,
        cn_max=cn_max,
        mu=mu,
        dT=dT,
        ampdel=ampdel,
        n_seed=n_seed,
        n_worker=n_worker,
        random_seed=random_seed,
        max_iters=max_iters,
        timelimit=timelimit
    )

    segmentation(
        cA,
        cB,
        u,
        cluster_ids,
        sample_ids,
        df=df,
        bbc_out_file=bbc_out_file
    )

    return obj, cA, cB, u, cluster_ids, sample_ids


def segmentation(
    cA,
    cB,
    u,
    cluster_ids,
    sample_ids,
    bbc_file,
    bbc_out_file=None
):
    df = pd.read_csv(bbc_file, sep='\t')
 
    n_clone = len(cA[0])
    cA = pd.DataFrame(cA, index=cluster_ids, columns=range(n_clone))
    cB = pd.DataFrame(cB, index=cluster_ids, columns=range(n_clone))
    u = pd.DataFrame(u, index=range(n_clone), columns=sample_ids)
    # Make (n_sample, n_clone) in shape; easier to merge later
    u = u.T

    # copy-numbers represented as <CN_A>|<CN_B> strings
    cN = cA.astype(str) + '|' + cB.astype(str)
    cN.columns = ['cn_normal'] + [f'cn_clone{i}' for i in range(1, n_clone)]

    # Merge in copy-number + proportion information to our original Dataframe
    df = df.merge(cN, left_on='CLUSTER', right_index=True)
    u.columns = ['u_normal'] + [f'u_clone{i}' for i in range(1, n_clone)]
    df = df.merge(u, left_on='SAMPLE', right_index=True)

    # Sorting the values by start/end position critical for merging contiguous
    # segments with identical copy numbers later on
    df = df.sort_values(['#CHR', 'START', 'END', 'SAMPLE'])
    df = df.reset_index(drop=True)

    # last 2*n_clone columns names = [cn_normal, u_normal, cn_clone1, u_clone1, cn_clone2, ...]
    extra_columns = [col for sublist in zip(cN.columns, u.columns) for col in sublist]
    all_columns = df.columns.values[: -2 * n_clone].tolist() + extra_columns

    if bbc_out_file is not None:
        # rearrange columns for easy comparison to legacy files
        df = df[all_columns]
        df.to_csv(bbc_out_file, sep='\t', index=False)


def runningTetraploid(
    outprefix,
    df,
    minC=2, 
    maxC=12, 
    clonal,
    scale,
    solver='gurobi',
    solve_mode='cd',
    d=-1,
    cn_max=None,
    mu=0.1,
    dT=0.1,
    ampdel=FALSE,
    n_seed=10,
    n_worker=6,
    random_seed=None,
    max_iters=None,
    timelimit=None,
    rundir=None
):

    cn = '{}:{}:{}'.format(scale[0], clonal[scale[0]][0], clonal[scale[0]][1])
    cn += ',{}:{}:{}'.format(scale[1], clonal[scale[1]][0], clonal[scale[1]][1])
    if len(clonal) > 2:
        cn = cn + ',' + ','.join(['{}:{}:{}'.format(s, clonal[s][0], clonal[s][1]) for s in clonal if s not in scale])

    obj={}
    for n in range(minC, maxC + 1):
        outprefix = os.path.join(rundir, 'results.tetraploid.n{}'.format(n))
        
        obj[n]=execute_solver(
                    outprefix=outprefix,
                    df=df,
                    clonal=cn,
                    solver=solver,
                    solve_mode=solve_mode,
                    d=d,
                    cn_max=cn_max,
                    mu=mu,
                    dT=dT,
                    ampdel=ampdel,
                    n_seed=n_seed,
                    n_worker=n_worker,
                    random_seed=random_seed,
                    max_iters=max_iters,
                    timelimit=timelimit
        )

    return obj
