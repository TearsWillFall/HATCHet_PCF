import pandas as pd
import numpy as np
import os
from hatchet.utils.Supporting import ensure, log, error


def findNeutralCluster(seg,tD):
    selected={}
    for id, df in seg.groupby('ID'):
        if all(abs(0.5-df["BAF"])<=tD):
            selected[id]=df["SIZE"].unique()
    if len(selected) == 0:
        raise ValueError(error('None potential neutral cluster found with given parameters!'))
    neutral = max(selected, key=(lambda k: selected[k]))

    return neutral



def findClonalClusters(
    seg,
    neutral,
    tB,
    tR,
    tS
):  
    topbaf={}
    heigh={}
    for sample, df in seg.groupby('SAMPLE'):
        topbaf[sample]=df.BAF.min()
        heigh[sample]=(abs(df[df.ID==neutral].BAF-topbaf[sample])/4.0)>tB
  
    location = {}
    level = {}

    
    for id, df in seg.groupby('ID'):
        locp = []
        levp = {}
        samples=[]
        for sample, df in seg.groupby('SAMPLE'):
            t = topbaf[sample]
            n = seg[(seg.ID==neutral)&(seg.SAMPLE==sample)].BAF
            x = df[(seg.ID==id)&(seg.SAMPLE==sample)].BAF
            d = {
                'top': abs(t - x),
                'mid': abs(((t + n) / 2.0) - x),
                'midbot': abs(((t + 3 * n) / 4.0) - x),
                'bot': abs(n - x),
            }

            if d['top'] < d['bot'] and abs(d['top'] - d['bot']) > tB:
                levp[sample] = ['top']
            elif d['bot'] < d['top'] and abs(d['top'] - d['bot']) > tB:
                levp[sample] = ['bot']
            else:
                levp[sample] = ['top', 'bot']

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
    clusters=clusters[~clusters.ID.isin(location)]

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

    samples=df.groupby('SAMPLE')
    ids=df.groupby('ID')

    for id, df in clusters.groupby('ID'):
        
        rightpos =sum(df.RDR-seg[seg.ID==neutral].RDR >tR)
        leftpos =sum(df.RDR-seg[seg.ID==neutral].RDR <-tR)


        if rightpos >tS:
            if location[id] == 'bot' or location[id] == 'midbot':
                options = [(3, 2), (4, 2)]
            elif location[id] == 'mid' or location[id] == 'top':
                options = [(4, 2), (3, 2)]
            else:
                assert False
            # if sum(abs(fseg[cluster][p]['baf'] - fseg[neutral][p]['baf']) <= tB for p in samples) == len(samples):
            #     options = [(4, 4)] + options
        elif leftpos == len(samples):
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
            purity = {p: calcPurity(seg[(seg.ID==neutral)&(seg.ID==p)].RDR, sum(opt), df[df.SAMPLE==p].RDR) for p in samples}
            if False in set(0.0 <= purity[p] <= 1.0 for p in samples):
                continue
            scaling = {p: calcScalingFactor(purity[p], seg[(seg.ID==neutral)&(seg.ID==p)].RDR) for p in samples}
            if False in set(scaling[p] >= 0.0 for p in samples):
                continue
            curr_pattern = {}
            curr_pattern[neutral] = (2, 2)
            curr_pattern[id] = opt
            curr_scale = (neutral, id)
            curr_value = seg[(seg.ID==neutral)].SIZE.unique() + df.SIZE.unique()
            for cn in [a for a in allclonal if a != opt]:
                estRDR = {p: calcRDR(purity[p], cn, scaling[p]) for p in samples}
                estBAF = {p: calcBAF(purity[p], cn) for p in samples}
                checkRDR = lambda r, i: set(p for p in samples if abs(r[p] - seg[(seg.ID==i)&(seg.SAMPLE==p)].RDR) <= tR) == samples
                checkBAF = lambda b, i: set(p for p in samples if abs(b[p] -  seg[(seg.ID==i)&(seg.SAMPLE==p)].BAF) <= tB) == samples
                candidates = [
                    idx
                    for idx in ids
                    if idx not in curr_pattern and checkRDR(estRDR, idx) and checkBAF(estBAF, idx)
                ]
                if len(candidates) > 0:
                    choices=seg[seg.ID.isin(candidates)]
                    choice=choices[choices.SIZE.idxmax()].ID.unique()
                    curr_pattern[choice] = cn
                    curr_value += choices.max()

            if curr_pattern not in found_pattern:
                found_pattern.append(curr_pattern)

            if curr_value > best_value:
                best_pattern = curr_pattern
                best_value = curr_value
                best_scale = curr_scale

    return best_pattern, best_scale
