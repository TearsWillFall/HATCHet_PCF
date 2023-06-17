
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
from find_cluster_pcf import findNeutralCluster,findClonalClusters


def comput_cn_pcf (
    bbfile,
    minK=2,
    maxK=30,
    covar="diag",
    decode_alg="viterbi",
    tmat="diag",
    tau=10e-6,
    restarts=10,
    subset=None,
    rundir=None,
    tD=0.1,
    tB=0.04,
    tR=0.08,
    tS=0.75
):

    seg=cluster_bins(
        bbfile=bbfile,
        minK=minK,
        maxK=maxK,
        covar=covar,
        decode_alg=decode_alg,
        tmat=tmat,
        tau=tau,
        restarts=restarts,
        subset=subset,
        rundir=rundir
    )
    
    neutral = findNeutralCluster(
        seg=seg,
        tD=tD
    )

    clonal, scale = findClonalClusters(
            seg=seg,
            neutral=neutral,
            tB=tB,
            tR=tR,
            tS=tS
        )
















