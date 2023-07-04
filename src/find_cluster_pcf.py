import pandas as pd
import numpy as np
import os
from collections import Counter
from hatchet.utils.Supporting import ensure, log, error
from util_pcf import argmax,argmin


def findNeutralCluster(seg,tD=0.1):
    selected={}
    for id, df in seg.groupby('ID'):
        if all(abs(0.5-df["BAF"])<=tD):
            selected[id]=df["SIZE"].unique()
    if len(selected) == 0:
        raise ValueError(error('None potential neutral cluster found with given parameters!'))
    neutral = max(selected, key=(lambda k: selected[k]))

    return neutral



def findClonalCluster(
    seg,
    neutral,
    tB=0.04,
    tR=0.08,
    tS=0.75
):  
    seg=seg.sort_values(by=['ID','SAMPLE','SIZE'])
    samples=seg.SAMPLE.unique()
    ids=seg.ID.unique()


    topbaf={}
    heigh={}

    topbaf = {p: seg[(seg.SAMPLE==p)].BAF.min() for p in samples}
    heigh = {p: (abs(seg[(seg.ID==neutral)&(seg.SAMPLE==p)].BAF.iloc[0] - topbaf[p]) / 4.0) > tB for p in samples}


    location = {}
    level = {}

    for id in ids:
        locp = []
        levp = {}
        for p in samples:
            t = topbaf[p]
            n = seg[(seg.ID==neutral)&(seg.SAMPLE==p)].BAF.iloc[0]
            x = seg[(seg.ID==id)&(seg.SAMPLE==p)].BAF.iloc[0]

            d = {
                'top': abs(t - x),
                'mid': abs(((t + n) / 2.0) - x),
                'midbot': abs(((t + 3 * n) / 4.0) - x),
                'bot': abs(n - x),
            }
            

            if d['top'] < d['bot'] and abs(d['top'] - d['bot']) > tB:
                levp[p] = ['top']
            elif d['bot'] < d['top'] and abs(d['top'] - d['bot']) > tB:
                levp[p] = ['bot']
            else:
                levp[p] = ['top', 'bot']

            c = argmin(d)
            locp.append(c)
            if c != 'top' and d['top'] <= tB:
                locp.append('top')
            if c != 'mid' and d['mid'] <= tB:
                locp.append('mid')
            if c != 'midbot' and d['midbot'] <= tB:
                locp.append('midbot')
            if c != 'bot' and d['bot'] <= tB:
                locp.append('bot')
        count = Counter([levv for p in samples for levv in levp[p]])
        count = sorted(count.keys(), key=(lambda x: count[x]), reverse=True)

        for lev in count:
            if False not in set(lev in levp[p] for p in samples):
                level[id] = lev
                locp = Counter(locp)
                loc = argmax(locp)
                location[id] = loc
                break


    clusters=seg[seg.ID!=neutral]
    clusters=clusters[clusters.ID.isin(location.keys())]
    cluster_ids=clusters.ID.unique()
    allclonal = [
    (2, 0),
    (2, 1),
    (3, 2),
    (4, 2),
    (1, 0),
    (3, 0),
    (3, 1),
    (4, 0),
    (4, 1),
    (5, 0),
    ]
    found_pattern = []
    best_pattern = {}
    best_scale = ()
    best_value = 0

    for id in cluster_ids:
        rightpos = sum((seg[(seg.ID==id)&(seg.SAMPLE==p)].RDR.iloc[0] - seg[(seg.ID==neutral)&(seg.SAMPLE==p)].RDR.iloc[0]) > tR for p in samples)
        leftpos = sum((seg[(seg.ID==id)&(seg.SAMPLE==p)].RDR.iloc[0] - seg[(seg.ID==neutral)&(seg.SAMPLE==p)].RDR.iloc[0]) < -tR for p in samples)

        eqbaf = lambda p: abs(seg[(seg.ID==id)&(seg.SAMPLE==p)].BAF.iloc[0] - seg[(seg.ID==neutral)&(seg.SAMPLE==p)].BAF.iloc[0]) <= tB
        eqrdr = lambda p: abs(seg[(seg.ID==id)&(seg.SAMPLE==p)].RDR.iloc[0] - seg[(seg.ID==neutral)&(seg.SAMPLE==p)].RDR.iloc[0]) <= tR

        if True in set(eqbaf(p) and heigh[p] for p in samples):
            continue
        if True in set(eqbaf(p) and eqrdr(p) and heigh[p] for p in samples):
            continue

        if rightpos/len(samples) >tS:
            if location[id] == 'bot' or location[id] == 'midbot':
                options = [(3, 2), (4, 2)]
            elif location[id] == 'mid' or location[id] == 'top':
                options = [(4, 2), (3, 2)]
            else:
                assert False
            # if sum(abs(fseg[cluster][p]['baf'] - fseg[neutral][p]['baf']) <= tB for p in samples) == len(samples):
            #     options = [(4, 4)] + options
        elif leftpos/len(samples) >tS:
            if level[id] == 'top' and location[id] == 'top':
                options = [(2, 0)]
            elif level[id] == 'top':
                options = [(2, 0), (2, 1)]
            elif level[id] == 'bot':
                options = [(2, 1)]
            else:
                assert False
            # if sum(abs(fseg[cluster][p]['baf'] - fseg[neutral][p]['baf']) <= tB for p in samples) == len(samples):
            #     options = [(1, 1), (0, 0)] + options
        else:
            options = []
            continue
   
        regPurity = lambda v: 1.0 if 1.0 <= v <= 1.05 else (0.0 if -0.05 <= v <= 0.0 else v)
        calcPurity = lambda d, c, r: regPurity(float(2 * d - 2 * r) / float(2 * r + 2 * d - c * d))
        calcScalingFactor = lambda p, d: float(2.0 + 2.0 * p) / float(d)
        calcFraction = lambda p, cn: float(2.0 * (1.0 - p) + sum(cn) * p)
        calcRDR = lambda p, cn, s: calcFraction(p, cn) / float(s)
        calcBAF = lambda p, cn: float(1.0 * (1.0 - p) + min(cn) * p) / calcFraction(p, cn)

        for opt in options:
            purity = {p: calcPurity(seg[(seg.ID==neutral)&(seg.SAMPLE==p)].RDR.iloc[0], sum(opt), seg[(seg.ID==id)&(seg.SAMPLE==p)].RDR.iloc[0]) for p in samples}
            if False in set(0.0 <= purity[p] <= 1.0 for p in samples):
                continue
            scaling = {p: calcScalingFactor(purity[p], seg[(seg.ID==neutral)&(seg.SAMPLE==p)].RDR.iloc[0]) for p in samples}
            if False in set(scaling[p] >= 0.0 for p in samples):
                continue
            curr_pattern = {}
            curr_pattern[neutral] = (2, 2)
            curr_pattern[id] = opt
            curr_scale = (neutral, id)
            curr_value = seg[(seg.ID==neutral)].SIZE.unique() + seg[(seg.ID==id)].SIZE.unique()
            for cn in [a for a in allclonal if a != opt]:
                estRDR = {p: calcRDR(purity[p], cn, scaling[p]) for p in samples}
                estBAF = {p: calcBAF(purity[p], cn) for p in samples}
                candidates=[]
                for idx in ids:
                    checkRDR=0
                    checkBAF=0
                    for p in samples:
                        checkRDR+=(estRDR[p]-seg[(seg.ID==idx)&(seg.SAMPLE==p)].RDR.iloc[0])<= tR
                        checkBAF+=(estBAF[p]-seg[(seg.ID==idx)&(seg.SAMPLE==p)].BAF.iloc[0])<= tB
                    if idx not in curr_pattern and checkRDR==len(samples) and checkBAF==len(samples):
                        candidates.append(idx)
                if len(candidates) > 0:
                    choice = max(candidates, key=(lambda i: seg[seg.ID==i].SIZE.unique()))
                    curr_pattern[choice] = cn
                    curr_value += seg[seg.ID==choice].SIZE.unique()

            if curr_pattern not in found_pattern:
                found_pattern.append(curr_pattern)

            if curr_value > best_value:
                best_pattern = curr_pattern
                best_value = curr_value
                best_scale = curr_scale

    return best_pattern,best_scale