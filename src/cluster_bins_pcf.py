from collections import Counter
import numpy as np
import os
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from scipy.special import logsumexp
from scipy.spatial.distance import pdist, squareform
from hmmlearn import hmm

from hatchet.utils.ArgParsing import parse_cluster_bins_args
import hatchet.utils.Supporting as sp

from read_pcf import read_bb

def cluster_bins(
    bbfile,
    minK=2,
    maxK=30,
    covar="diag",
    decode_alg="viterbi",
    tmat="diag",
    tau=10e-6,
    restarts=10,
    subset=None,
    rundir=""
):
    tracks, bb, sample_labels, chr_labels =read_bb(bbfile,subset=subset)

    outbb= os.path.join(rundir, 'output.bbc')

    bb.to_csv(outbb,index=False, sep=',')

    if minK==None:
        minK=len(chr_labels)
    elif maxK==None:
        maxK=len(bb.GENE.unique())

    (best_score, best_model, best_labels, best_K, results) =hmm_model_select(
            tracks,
            minK=minK,
            maxK=maxK,
            covar=covar,
            decode_alg=decode_alg,
            tmat=tmat,
            tau=tau,
            restarts=restarts,
        )

    best_labels = reindex(best_labels)
    bb['CLUSTER'] = np.repeat(best_labels, len(sample_labels))
    seg = form_seg(bb)
    outseg= os.path.join(rundir, 'output.seg')
    seg.to_csv(outseg, index=False, sep='\t')

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


class DiagGHMM(hmm.GaussianHMM):
    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice):
        super()._accumulate_sufficient_statistics(stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice)

        if 't' in self.params:
            # for each ij, recover sum_t xi_ij from the inferred transition matrix
            bothlattice = fwdlattice + bwdlattice
            loggamma = (bothlattice.T - logsumexp(bothlattice, axis=1)).T

            # denominator for each ij is the sum of gammas over i
            denoms = np.sum(np.exp(loggamma), axis=0)
            # transpose to perform row-wise multiplication
            stats['denoms'] = denoms

    def _do_mstep(self, stats):
        super()._do_mstep(stats)
        if 't' in self.params:

            denoms = stats['denoms']
            x = (self.transmat_.T * denoms).T

            # numerator is the sum of ii elements
            num = np.sum(np.diag(x))
            # denominator is the sum of all elements
            denom = np.sum(x)

            # (this is the same as sum_i gamma_i)
            # assert np.isclose(denom, np.sum(denoms))

            stats['diag'] = num / denom
            # print(num.shape)
            # print(denom.shape)

            self.transmat_ = self.form_transition_matrix(stats['diag'])

    def form_transition_matrix(self, diag):
        tol = 1e-10
        diag = np.clip(diag, tol, 1 - tol)

        offdiag = (1 - diag) / (self.n_components - 1)
        transmat_ = np.diag([diag - offdiag] * self.n_components)
        transmat_ += offdiag
        # assert np.all(transmat_ > 0), (diag, offdiag, transmat_)
        return transmat_



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
