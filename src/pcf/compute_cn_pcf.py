
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
from solve_pcf import solve_pcf
import find_cluster_pcf


def comput_cn_pcf (
    args
):

    seg=cluster_bins(
        bbfile=args['bbfile'],
        args=args
    )
    
     neutral = findNeutralCluster(
        seg=seg,
        tD=args['tD']
    )

    clonal, scale = findClonalClusters(
            seg=seg,
            neutral=neutral,
            tB=args['tB'],
            tR=args['tR'],
            tR=args['tS']
        )
















