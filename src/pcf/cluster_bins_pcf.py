from collections import Counter
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from scipy.special import logsumexp
from scipy.spatial.distance import pdist, squareform
from hmmlearn import hmm

from hatchet.utils.ArgParsing import parse_cluster_bins_args
import hatchet.utils.Supporting as sp

import read_pcf


def cluster_bins(
    bbfile,
    args
):
    tracks, bb, sample_labels, chr_labels =read_pcf(bbfile,args['subset'])

    bb.to_csv(args['outbbc'], index=False, sep='\t')

    (best_score, best_model, best_labels, best_K, results) =hmm_model_select(
            tracks,
            minK=args['minK'],
            maxK=args['maxK'],
            covar=args['covar'],
            decode_alg=args['decoding'],
            tmat=args['transmat'],
            tau=args['tau'],
            restarts=args['restarts'],
        )

    best_labels = reindex(best_labels)
    bb['CLUSTER'] = np.repeat(best_labels, len(sample_labels))
    seg = form_seg(bb)
    seg.to_csv(args['outsegments'], index=False, sep='\t')
    return seg



def hmm_model_select(tracks, minK=20, maxK=50, tau=10e-6, tmat='diag', decode_alg='viterbi', covar='diag', restarts=10):
    assert tmat in ['fixed', 'diag', 'free']
    assert decode_alg in ['map', 'viterbi']

    # format input
    tracks = [a for a in tracks if a.shape[0] > 0 and a.shape[1] > 0]
    if len(tracks) > 1:
        X = np.concatenate(tracks, axis=1).T
        lengths = [a.shape[1] for a in tracks]
    else:
        X = tracks[0].T
        lengths = [tracks[0].shape[1]]

    best_K = 0
    best_score = -1.01   # below minimum silhouette score value
    best_model = None
    best_labels = None

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    C = squareform(pdist(X_scaled))

    rs = {}
    for K in range(minK, maxK + 1):
        # print(K, datetime.now())

        my_best_ll = -1 * np.inf
        my_best_labels = None
        my_best_model = None
        for s in range(restarts):
            # construct initial transition matrix
            A = make_transmat(1 - tau, K)
            assert np.all(A > 0), (
                'Found 0 or negative elements in transition matrix.'
                'This is likely a numerical precision issue -- try increasing tau.',
                A,
            )
            assert np.allclose(np.sum(A, axis=1), 1), ('Not all rows in transition matrix sum to 1.', A)

            if tmat == 'fixed':
                model = hmm.GaussianHMM(
                    n_components=K,
                    init_params='mc',
                    params='smc',
                    covariance_type=covar,
                    random_state=s,
                )
            elif tmat == 'free':
                model = hmm.GaussianHMM(
                    n_components=K,
                    init_params='mc',
                    params='smct',
                    covariance_type=covar,
                    random_state=s,
                )
            else:
                model = DiagGHMM(
                    n_components=K,
                    init_params='mc',
                    params='smct',
                    covariance_type=covar,
                    random_state=s,
                )

            model.startprob_ = np.ones(K) / K
            model.transmat_ = A
            model.fit(X, lengths)

            prob, labels = model.decode(X, lengths, algorithm=decode_alg)
            if prob > my_best_ll:
                my_best_labels = labels
                my_best_ll = prob
                my_best_model = model

        score = silhouette_score(C, my_best_labels, metric='precomputed')

        rs[K] = my_best_ll, score, my_best_labels
        if score > best_score:
            best_score = score
            best_model = my_best_model
            best_labels = my_best_labels
            best_K = K

    return best_score, best_model, best_labels, best_K, rs





def make_transmat(diag, K):
    offdiag = (1 - diag) / (K - 1)
    transmat_ = np.diag([diag - offdiag] * K)
    transmat_ += offdiag
    return transmat_


def reindex(labels):
    """
    Given a list of labels, reindex them as integers from 1 to n_labels
    Also orders them in nonincreasing order of prevalence
    """
    old2new = {}
    j = 1
    for i, _ in Counter(labels).most_common():
        old2new[i] = j
        j += 1
    old2newf = lambda x: old2new[x]

    return [old2newf(a) for a in labels]


def form_seg(bbc):
    segments = []
    for key, df_ in bbc.groupby('CLUSTER'):
        nbins = []
        size=[]
        chr = []
        genes = []
        rdr = []
        beta=[]
        baf = []
        samples=[]

        for sample, df in df_.groupby('SAMPLE'):
            size.append(df.SIZE.sum())
            nbins.append(len(df))
            chr.append(df.CHR_ARM.unique())
            genes.append(df.GENE.unique())
            rdr.append(df.RDR.mean())
            baf.append(df.BAF.mean())
            beta.append(df.BETA.mean())
            samples.append(sample)
        keys = [key] * len(baf)
        [segments.append(t) for t in zip(keys,samples,chr,genes,nbins,size, rdr, baf, beta)]

    seg = pd.DataFrame(
        segments,
        columns=[
            'ID',
            'SAMPLE',
            'CHR',
            'GENES',
            'BINS',
            'SIZE',
            'RDR',
            'BAF',
            'BETA'
        ],
    )
    return seg
