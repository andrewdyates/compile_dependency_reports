#!/usr/bin/python
"""Plot histogram of numpy array.

python $HOME/compile_dependency_reports/enrichment_hist.py json_input=$HOME/gse2034/gse2034.json outdir=$HOME pina_file=$HOME/PINA_gene_folder_enrichment/Homo-sapiens-20110628.txt measures=mic
"""
import sys
import random
import json
import numpy as np
import numpy.ma as ma
from py_symmetric_matrix import *
from enriched import *
from util import *
import errno, os

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

r_base = importr('base')
grdevices = importr('grDevices')
stats = importr('stats')
graphics = importr('graphics')

def main(json_input=None, outdir=None, pina_file=None, measures=None):
  assert json_input and outdir and pina_file
  if type(measures) == str:
    measures = measures.split(',')

  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

  # Load JSON file
  D = json.load(open(json_input))

  varlist = []
  Q = np.genfromtxt(name_iter(open(D['data_file']), varlist), usemask=True, delimiter='\t')

  # Load Enrichments and create mask
  Enrich = PINAEnriched(open(pina_file))
  RandEnrich = RandomPINAEnriched(open(pina_file))
  
  E_Mask = np.zeros(D['n_pairs'], dtype=np.bool)
  np.put(E_Mask, Enrich.indices(varlist), True)
  R_Mask = np.zeros(D['n_pairs'], dtype=np.bool)
  np.put(R_Mask, RandEnrich.indices(varlist), True)
  print "Loaded %d pairs from %d variables from %s." % \
    (len(Enrich.pairs), len(Enrich.genes), pina_file)

  # Per dependency, compare histograms
  for name, d in D['dependencies'].items():
    if measures is not None:
      if name not in measures:
        continue
        
    print "Loading %s..." % name
    print "values: ", d
    M = np.load(os.path.join(d['dir'], d['values_file']))
    B = np.load(os.path.join(d['dir'], d['bool_file']))
    n = np.size(M,0)
    n_missing = n - np.count_nonzero(B)
    print "%s missing %d pairs of %d..." % (name, n_missing, n)
    assert n == D['n_pairs']

    # B: 1 is valid. E_Mask: 1 is selected. ~(B&E_mask)
    #M_E = ma.array(data=M, mask=~(B&E_Mask))
    #M_R = ma.array(data=M, mask=~(B&R_Mask))
    M_R = np.compress((B&R_Mask), M)
    M_E = np.compress((B&E_Mask), M)
    # perhaps not the most efficient...
    a=random.sample(xrange(np.count_nonzero(B)), len(M_E))
    M_All = np.take(np.compress(B, M), a)
    print "Enriched size: %d. Randomized PINA size: %d. From All Random size: %d" % \
        (len(M_E), len(M_R), len(M_All))
    
    # Draw plot
    plot_pdfname = os.path.join(outdir, "%s_PINA_enriched_hist.pdf" % name)
    print "Plotting %s..." % plot_pdfname
    grdevices.pdf(plot_pdfname)
    graphics.plot(stats.density(M_E), main="%s PINA enriched Histogram" % name, xlab=name, col="blue")
    graphics.lines(stats.density(M_All), col="grey")
    graphics.lines(stats.density(M_R), col="black")
    graphics.rug(M_E, col="blue")
    graphics.rug(M_R, col="black")
    grdevices.dev_off()

    
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
  
  
