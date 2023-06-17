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



def read_bb(bbfile, subset=None):
    """
    Constructs arrays to represent the bin in each chromosome or arm.
    If bbfile was binned around chromosome arm, then uses chromosome arms.
    Otherwise, uses chromosomes.

    Returns:
        botht: list of np.ndarrays of size (n_bins, n_tracks)
            where n_tracks = n_samples * 2
        bb: table read from input bbfile
        sample_labels: order in which samples are represented in each array in botht
        chr_lables: order in which chromosomes or arms are represented in botht

    each array contains
    1 track per sample for a single chromosome arm.
    """

    bb = pd.read_table(bbfile)
    bb["BAF"]=bb["BETA"]/(1+bb["BETA"])
    bb["SIZE"]=bb["END"]-bb["START"]
    bb["CHR_ARM"]=bb["CHR"]+"_"+bb["ARM"]
    if subset is not None:
        bb = bb[bb.SAMPLE.isin(subset)]

    tracks = []

    sample_labels = []
    populated_labels = False

    chr_labels = []
    for ch, df0 in bb.groupby('CHR_ARM'):
        df0 = df0.sort_values('START')
        arrs=[]
        for sample, df in df0.groupby('SAMPLE'):
            df = df.sort_values('SAMPLE')
            if not populated_labels:
                sample_labels.append(sample)
              
            arrs.append(df.BAF.to_numpy())
            arrs.append(df.RDR.to_numpy())
    
            if len(arrs) > 0:
                tracks.append(np.array(arrs))
                chr_labels.append(str(ch))
        

        populated_labels = True

    return (
        tracks,
        bb.sort_values(by=['CHR', 'START', 'SAMPLE']),
        sample_labels,
        chr_labels,
    )





