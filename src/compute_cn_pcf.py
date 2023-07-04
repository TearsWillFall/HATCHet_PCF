
import os
import sys
import argparse
import shutil
import subprocess
import shlex
from collections import Counter

import hatchet
from hatchet import config, __version__
from hatchet.utils.Supporting import ensure, log, error
from hatchet.utils.solve import solve
from hatchet.utils.solve.utils import segmentation
from cluster_bins_pcf import cluster_bins
from find_cluster_pcf import findNeutralCluster,findClonalCluster
from execute_solver_pcf import runningDiploid,runningTetraploid
from select_pcf import select
import hatchet.utils.Supporting as sp

def compute_cn_pcf(
    bbfile,
    minK=None,
    maxK=None,
    covar="diag",
    decode_alg="map",
    tmat="diag",
    tau=0.000001,
    restarts=10,
    preserve_arm=True,
    subset=None,
    tD=0.1,
    tB=0.04,
    tR=0.08,
    tS=0.75,
    minC=2, 
    maxC=12, 
    solver='gurobi',
    solve_mode='cd',
    d=-1,
    cn_max=None,
    mu=0.1,
    dT=0.1,
    g=0.3,
    limit=0.6,
    ampdel=False,
    n_seed=10,
    n_worker=6,
    random_seed=None,
    max_iters=None,
    timelimit=None,
    rundir=""
):

    seg=cluster_bins(
        bbfile=bbfile,
        minK=minK,
        maxK=maxK,
        covar=covar,
        decode_alg=decode_alg,
        tmat=tmat,
        tau=tau,
        preserve_arm=preserve_arm,
        restarts=restarts,
        subset=subset,
        rundir=rundir
    )
    
    neutral = findNeutralCluster(
        seg=seg,
        tD=tD
    )

    sp.log(
        msg='Found neutral solutions = \n' + str(neutral) + '\n',
        level='INFO',
    )

    clonal, scale = findClonalCluster(
            seg=seg,
            neutral=neutral,
            tB=tB,
            tR=tR,
            tS=tS
        )


    sp.log(
        msg='Found clonal solutions = \n' + str(clonal) + '\n',
        level='INFO',
    )
    sp.log(
        msg='Found scales = \n' + str(scale) + '\n',
        level='INFO',
    )

    
    diploid=runningDiploid(
        df=seg,
        neutral=neutral,
        minC=minC, 
        maxC=maxC, 
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
        timelimit=timelimit,
        rundir=rundir)



    sp.log(
        msg='Diploid objective values = \n' + str(diploid) + '\n',
        level='INFO',
    )


    tetraploid=runningTetraploid(
        df=seg,
        clonal=clonal,
        scale=scale,
        minC=minC, 
        maxC=maxC, 
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
        timelimit=timelimit,
        rundir=rundir
    )



    sp.log(
        msg='Tetraploid objective values = \n' + str(tetraploid) + '\n',
        level='INFO',
    )

    dbest,tbest=select(
        diploid=diploid,
        tetraploid=tetraploid,
        g=g,
        limit=limit,
        rundir=rundir
    )
    return dbest,tbest


    















