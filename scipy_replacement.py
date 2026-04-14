"""
Replacement for scipy used by PhyloWGS when scipy is not available (e.g. Python 2.7 on Apple Silicon).
Uses only numpy and the standard library (math).
"""
from __future__ import division
import math
import numpy as np

# ---------------------------------------------------------------------------
# scipy.special.gammaln  ->  math.lgamma (log of gamma function)
# ---------------------------------------------------------------------------
def gammaln(x):
    """Log of the absolute value of the gamma function. Handles arrays via numpy."""
    try:
        iter(x)
    except TypeError:
        return float(math.lgamma(float(x)))
    x = np.asarray(x, dtype=float)
    if x.shape == ():
        return float(math.lgamma(float(x)))
    out = np.vectorize(math.lgamma)(x)
    return out

# ---------------------------------------------------------------------------
# scipy.misc.logsumexp  ->  stable log-sum-exp
# ---------------------------------------------------------------------------
def logsumexp(X, axis=None):
    """Stable log-sum-exp. X can be array; axis as in numpy."""
    X = np.asarray(X)
    maxes = np.max(X, axis=axis)
    return np.log(np.sum(np.exp(X - maxes), axis=axis)) + maxes

# ---------------------------------------------------------------------------
# scipy.stats.mstats.gmean  ->  geometric mean = exp(mean(log(x)))
# ---------------------------------------------------------------------------
def gmean(a, axis=0):
    """Geometric mean along axis. Supports masked arrays (masked values ignored)."""
    try:
        import numpy.ma as ma
        if isinstance(a, ma.MaskedArray):
            a = ma.array(a, dtype=float)
            return np.exp(ma.mean(ma.log(a), axis=axis))
    except (ImportError, AttributeError):
        pass
    a = np.asarray(a, dtype=float)
    return np.exp(np.mean(np.log(a), axis=axis))

# ---------------------------------------------------------------------------
# scipy.stats.gaussian_kde replacement: simple Gaussian KDE
# ---------------------------------------------------------------------------
class gaussian_kde(object):
    """Simple Gaussian kernel density estimator. Replaces scipy.stats.gaussian_kde for 1D/2D."""
    def __init__(self, dataset):
        self.dataset = np.atleast_2d(np.asarray(dataset))
        if self.dataset.shape[0] > self.dataset.shape[1]:
            self.dataset = self.dataset.T
        self.n, self.d = self.dataset.shape
        # Silverman's rule of thumb for bandwidth
        sigmas = np.std(self.dataset, axis=1, ddof=1)
        sigmas = np.where(sigmas > 1e-10, sigmas, 1.0)
        self.bw = sigmas * (self.n * (self.d + 2) / 4.0) ** (-1.0 / (self.d + 4))

    def __call__(self, points):
        points = np.atleast_2d(np.asarray(points))
        if points.shape[0] != self.d:
            points = points.T
        m = points.shape[1]
        density = np.zeros(m)
        for i in range(self.n):
            diff = (points - self.dataset[:, i:i + 1]) / (self.bw[:, np.newaxis] + 1e-300)
            density += np.exp(-0.5 * np.sum(diff ** 2, axis=0)) / (np.prod(self.bw) * (2 * np.pi) ** (self.d / 2.0))
        density /= self.n
        return density

# ---------------------------------------------------------------------------
# scipy.misc.comb  (n choose k)  ->  exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1))
# ---------------------------------------------------------------------------
def comb(n, k, exact=False):
    """Number of combinations. exact=True not supported; uses float."""
    n, k = int(n), int(k)
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    return int(round(np.exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))))
