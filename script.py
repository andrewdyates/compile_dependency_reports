#!/usr/bin/python
"""Given a set of dependency matrices for a microarray, generate a suite of analytics.

SAMPLE USE:

  python script.py json_input="gse2034.json" outdir="test" pina_file="/nfs/01/osu6683/PINA_gene_folder_enrichment/Homo-sapiens-20110628.txt"
"""
import json
import numpy as np
import numpy.ma as ma
from py_symmetric_matrix import *
from enriched import *
from util import *
import errno, os
import sys
import itertools
from scipy.stats import mstats

import matplotlib 
matplotlib.use('agg') # required for use on OSC servers
import matplotlib.pyplot as plt

TOP_PLOTS = 500
TOP_ENRICHMENTS = [1000, 10000, 100000, 500000]


def main(json_input=None, outdir=None, pina_file=None):
  assert (json_input and outdir is not None)

  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise
  print "Created output directory %s" % os.path.abspath(outdir)

  D = json.load(open(json_input))
  R = {} # results dictionary

  # create list of enrichment sets
  enriched = {'All': None}
  if pina_file:
    enriched['PINA'] = PINAEnriched(open(pina_file))
    print "Loaded %d pairs from %d variables from %s." % \
        (len(enriched['PINA'].pairs), len(enriched['PINA'].genes), pina_file)


  varlist = []
  Q = np.genfromtxt(name_iter(open(D['data_file']), varlist), usemask=True, delimiter='\t')
  print "Loaded %d variables of %d entries from %s." % (np.size(Q,0), np.size(Q,1), D['data_file'])
  print

  # For each dependency:
  for name, d in D['dependencies'].items():
    print "Loading %s..." % name
    print "values: ", d
    M = np.load(os.path.join(d['dir'], d['values_file']))
    B = np.load(os.path.join(d['dir'], d['bool_file']))
    A = np.load(os.path.join(d['dir'], d['argsorted_file']))
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
    Q = ma.array(data=M, mask=~B)
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

    if E is not None:
      # Generate enrichment mask
      E_Mask = np.zeros(np.size(Q), dtype=np.bool)
      set_list = E.indices(varlist)
      np.put(E_Mask, set_list, True)
    else:
      E_Mask = None
      # Mask by enrichment
    
    for dep_name, d in R.items():

      M = d['Q']
      if E_Mask is not None:
        M = M[E_Mask] # cache this?
        
      print enrich_name, dep_name
      stats = {
        'size': ma.count(M), # is this correct given masked values?
        'mean': M.mean(),
        'std': M.std(),
        'median': ma.median(M),
        }
      print stats
      d.update(stats)

      # generate histograms
      name_tuple = (enrich_name, dep_name, D['gse_id'])
      print "Plotting histogram for %s..." % "_".join(name_tuple)
      plt.clf(); plt.cla()
      plt.title("%s %s %s" % name_tuple)
      plt.hist(M.compressed(), bins=1000, histtype='step', normed=True)
      plot_path = os.path.join(outdir, "hist"+"_".join(name_tuple)+".png")
      print "Saving figure as %s..." % plot_path
      plt.savefig(plot_path, dpi=600)
      # log histogram
      print "Plotting log histogram for %s..." % "_".join(name_tuple)
      plt.clf(); plt.cla()
      plt.title("Logscale %s %s %s" % name_tuple)
      plt.hist(M.compressed(), bins=1000, histtype='step', normed=True, log=True)
      plot_path = os.path.join(outdir, "loghist"+"_".join(name_tuple)+".png")
      print "Saving figure as %s..." % plot_path
      plt.savefig(plot_path, dpi=600)


    # compute all-pairs 
    for x, y in itertools.combinations(R.keys(), 2):

      name_tuple = (x, y, enrich_name, D['gse_id'])
      print "Scatter Plotting %s versus %s for %s from %s" % (name_tuple)
      X_Q = R[x]['Q']
      Y_Q = R[y]['Q']
      if E_Mask is not None:
        X_Q = X_Q[E_Mask]
        Y_Q = Y_Q[E_Mask]
        
      pcc = mstats.pearsonr(X_Q, Y_Q)  
      print "PCC of %s and %s:" % (x,y), pcc
      
      plt.clf(); plt.cla()
      plt.title(" ".join(name_tuple))
      plt.xlabel(x); plt.ylabel(y)
      plot_path = os.path.join(outdir, "scatter"+"_".join(name_tuple)+".png")
      plt.plot(X_Q, Y_Q, 'b.')
      print "Saving figure as %s..." % plot_path
      plt.savefig(plot_path, dpi=600)
      

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
  main(**dict([s.split('=') for s in sys.argv[1:]]))
