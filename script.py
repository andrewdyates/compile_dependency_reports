#!/usr/bin/python
"""Given a set of dependency matrices for a microarray, generate a suite of analytics."""
import json
import numpy as np
import numpy.ma as ma
from py_symmetric_matrix import *
from enriched import *
from util import *
import errno, os
import sys

TOP_PLOTS = 500
TOP_ENRICHMENTS = [1000, 10000, 100000, 500000]


def main(json_input=None, outdir=None, pina_file=None):
  assert (json_input and outdir)

  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

  D = json.load(open(json_input))
  R = {} # results dictionary

  # create list of enrichment sets
  # think of some more elegant way to do this?
  enriched = {'None': None}
  if pina_file:
    enriched['PINA'] = PINAEnriched(open(pina_file))


  varlist = []
  Q = np.genfromtxt(name_iter(open(D['data_file']), varlist), usemask=True, delimiter='\t')
  print "Loaded %d variables of %d entries." % (np.size(Q,0), np.size(Q,1))
  print

  # For each dependency:
  for name, d in D['dependencies'].items():
    print "Loading %s..." % name
    print "values: ", d
    M = np.load(d['values_file'])
    B = np.load(d['bool_file'])
    A = np.load(d['argsorted_file'])
    # report missing values, get size
    n = np.size(M,0)
    try:
      assert n == np.size(B,0) == np.size(A,0)
    except AssertionError:
      print  "Fail: n[%d] == np.size(B,0)[%d] == np.size(A,0)[%d]" % (n, np.size(B,0), np.size(A,0))
      print "Skipping %s..." % name
      continue
    n_missing = np.size(np.where(B == 0))
    # create masked array 
    Q = ma.array(data=M, mask=B)
    # save results in results dict R
    R[name] = {
      'n': n,
      'n_missing': n_missing,
      'Q': Q,
      'A': A,
    }
    print "Loaded %s." % name
    print R[name]
    print
    

  # for each enriched set (including no enrichment)
  for enrich_name, E in enriched.items():

    # skip enrichment for now
    if E is not None:
      continue
    
    for name, d in R:
      print name
      stats = {
        'mean': d['Q'].mean(),
        'std': d['Q'].std(),
        'mode': d['Q'].mode(),
        'median': d['Q'].median(),
        }
      print stats
      d.update(stats)

    #   generate histogram
    #   if not no enrichment, compute enrichment
    #     output stats
    #     plot stats
    #   print top 500 pairs
    
    # for all dependencies:
    # plot histograms
    # plot logscale histogram
    # plot enrichments to 
      # 1k
      # 10k
      # 100k
      # 500k
    # plot 
    
    # for all pairs of dependencies:
      # compute pcc
      # compute dcor
      # compute spearman
      # scatterplot

if __name__ == "__main__":
  main(dict(**[s.split('=') for s in sys.argv[1]]))
