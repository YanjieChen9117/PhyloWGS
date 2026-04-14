from pwgsresults.index_calculator import IndexCalculator
import json
import gzip
import zipfile
import numpy as np
try:
    import scipy.stats
    _has_scipy = True
except ImportError:
    _has_scipy = False
    import sys
    import os
    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _root not in sys.path:
        sys.path.insert(0, _root)

np.seterr(invalid='raise')

def _gaussian_kde(dataset):
    if _has_scipy:
        return scipy.stats.gaussian_kde(dataset)
    from scipy_replacement import gaussian_kde
    return gaussian_kde(dataset)

def calc_tree_densities(summaries):
  tidxs = sorted(summaries.keys())
  _extract = lambda idxname: np.array([summaries[tidx][idxname] for tidx in tidxs])

  epsilon = 0.0001
  indices = {I: _extract(I + '_index') for I in ('linearity', 'branching', 'clustering')}
  X = indices['clustering']
  # Epsilon prevents division by zero in case of single-node trees.
  Y = indices['branching'] / (indices['branching'] + indices['linearity'] + epsilon)

  # Must be (# dimensions, # data points)
  XY = np.vstack((X, Y))
  # Must conver to Python list so it can be serialized to JSON.

  try:
    density = list(_gaussian_kde(XY)(XY))
  except (np.linalg.LinAlgError, FloatingPointError, AttributeError):
    try:
      density = list(_gaussian_kde(X)(X))
    except (np.linalg.LinAlgError, FloatingPointError, AttributeError):
      density = np.zeros(len(X))

  return dict(zip(tidxs, density))

class JsonWriter(object):
  def __init__(self, dataset_name):
    self._dataset_name = dataset_name

  def write_mutlist(self, mutlist, mutlist_outfn):
    with gzip.GzipFile(mutlist_outfn, 'w') as mutf:
      mutlist['dataset_name'] = self._dataset_name
      json.dump(mutlist, mutf)

  def write_summaries(self, summaries, params, summaries_outfn):
    for summary in summaries.values():
      calculator = IndexCalculator(summary)
      summary['linearity_index'] = calculator.calc_linearity_index()
      summary['branching_index'] = calculator.calc_branching_index()
      summary['clustering_index'] = calculator.calc_clustering_index()

    to_dump = {
      'dataset_name': self._dataset_name,
      'params': params,
      'trees': summaries,
      'tree_densities': calc_tree_densities(summaries),
    }
    with gzip.GzipFile(summaries_outfn, 'w') as summf:
      json.dump(to_dump, summf)

  def write_mutass(self, mutass, mutass_outfn):
    with zipfile.ZipFile(mutass_outfn, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as muts_file:
      for tree_idx, tree_mutass in mutass.items():
        to_dump = {
          'mut_assignments': tree_mutass,
          'dataset_name': self._dataset_name
        }
        muts_file.writestr('%s.json' % tree_idx, json.dumps(to_dump))

