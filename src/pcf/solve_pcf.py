
import os
import numpy as np
import pandas as pd
from pyomo import environ as pe
from collections import OrderedDict

from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit
from hatchet.utils.solve.cd import CoordinateDescent, CoordinateDescentSplit
from hatchet.utils.solve.utils import parse_clonal, scale_rdr
from hatchet import config
import hatchet.utils.Supporting as sp

def solver_available(solver=None):
    solver = solver or config.compute_cn.solver
    if solver == 'cpp':
        return os.getenv('GRB_LICENSE_FILE') is not None
    elif solver == 'gurobipy':
        return pe.SolverFactory('gurobi', solver_io='python').available(exception_flag=False)
    return pe.SolverFactory(solver).available(exception_flag=False)


def solve_pcf(
    n,
    df,
    clonal=None,
    solver='gurobi',
    solve_mode='cd',
    d=-1,
    cn_max=None,
    mu=0.01,
    dT=0.1,
    ampdel=True,
    n_seed=400,
    n_worker=8,
    random_seed=None,
    max_iters=None,
    timelimit=None
):

    assert solve_mode in ('ilp', 'cd', 'both'), 'Unrecognized solve_mode'
    assert solver_available(solver), f'Solver {solver} not available or not licensed'

    if max_iters is None:
        max_iters = 10

    df = df.sort_values(['ID', 'SAMPLE'])

    # sanity-check
    sample_ids = np.sort(df['SAMPLE'].unique())
    for _cluster_id, _df in df.groupby('ID'):
        _sample_ids = _df['SAMPLE'].values
        if not np.all(_sample_ids == sample_ids):
            raise ValueError(f'Sample IDs for cluster {_cluster_id} do not match {sample_ids}')

    rdr = df.pivot(index='ID', columns='SAMPLE', values='RDR')
    baf = df.pivot(index='ID', columns='SAMPLE', values='BAF')

    bins = {}  # cluster_id => no. of bins
    for cluster_id, _df in df.groupby('ID'):
        if len(_df['BINS'].unique()) != 1:
            raise ValueError(f'Bin sizes for cluster {cluster_id} across tumor samples are not identical!')
        bins[cluster_id] = _df.iloc[0]['BINS']
    bins = pd.Series(bins)

    weights = 100 * bins / sum(bins)

    if clonal is None:
        _candidate_cluster_ids = np.all(rdr > 0.5 - dT, axis=1)
        if not np.any(_candidate_cluster_ids):
            raise RuntimeError(f'Unable to determine cluster with diploid RDR threshold {dT}')
        dipoid_cluster_id = (_candidate_cluster_ids * weights).idxmax()
        copy_numbers = {dipoid_cluster_id: (1, 1)}
    else:
        copy_numbers = parse_clonal(clonal)

    gamma = scale_rdr(rdr, copy_numbers)

    
 

    rdr = rdr * gamma
    f_a = rdr * baf
    f_b = rdr - f_a

    gamma=gamma.to_frame()

    seg=seg.set_index(['ID', 'SAMPLE'])

    seg=seg.merge(gamma, left_index=True, right_index=True)

    seg['CN']=seg['RDR']*seg["GAMMA"]
    seg['CN.int']=seg['CN'].round()
    seg['fA']=seg['CN']*seg["BAF"]
    seg['fA.int']=seg['fA'].round()
    seg['fB']=seg['CN']-seg["fA"]
    seg['fB.int']=seg['fB'].round()

    if cn_max==None:
        cn_max=np.ceil(seg['CN'])


    sp.log(
        msg='Maximum copy-numner = \n' + str(cn_max) + '\n',
        level='INFO',
    )

    if solve_mode == 'ilp':
        ilp = ILPSubset(
            n,
            cn_max,
            d=d,
            mu=mu,
            ampdel=ampdel,
            copy_numbers=copy_numbers,
            f_a=f_a,
            f_b=f_b,
            w=weights,
        )
        ilp.create_model(pprint=True)
        return ilp.run(solver_type=solver, timelimit=timelimit)
    elif solve_mode == 'cd':
        cd = CoordinateDescent(
            f_a=f_a,
            f_b=f_b,
            n=n,
            mu=mu,
            d=d,
            cn_max=cn_max,
            w=weights,
            ampdel=ampdel,
            cn=copy_numbers,
        )
        return cd.run(
            solver_type=solver,
            max_iters=max_iters,
            n_seed=n_seed,
            j=n_worker,
            random_seed=random_seed,
            timelimit=timelimit,
        )
    else:
        cd = CoordinateDescent(
            f_a=f_a,
            f_b=f_b,
            n=n,
            mu=mu,
            d=d,
            cn_max=cn_max,
            w=weights,
            ampdel=ampdel,
            cn=copy_numbers,
        )
        _, cA, cB, _, _, _ = cd.run(
            solver_type=solver,
            max_iters=max_iters,
            n_seed=n_seed,
            j=n_worker,
            random_seed=random_seed,
            timelimit=timelimit,
        )

        ilp = ILPSubset(
            n,
            cn_max,
            d=d,
            mu=mu,
            ampdel=ampdel,
            copy_numbers=copy_numbers,
            f_a=f_a,
            f_b=f_b,
            w=weights,
        )
        ilp.create_model()
        ilp.hot_start(cA, cB)
        return ilp.run(solver_type=solver, timelimit=timelimit)

