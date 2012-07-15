#!/usr/bin/python
"""Generate scatterplots of top k pairs."""
from py_symmetric_matrix import *
from matplotlib import pyplot as pp
import numpy as np
import json
import os, sys, errno


def generate_top_k_scatters(M, D, A, varlist, dep_name, study_id, k=500):
  for i in xrange(1,k+1):
    n = len(varlist)
    idx = A[-i]
    score = D[idx]
    xi, yi = inv_sym_idx(idx, n)
    x, y = varlist[xi], varlist[yi]
    pp.clf(); pp.cla();
    pp.title("%d %s %.3f" % (i, dep_name, score))
    pp.xlabel(x); pp.ylabel(y)
    pp.plot(M[xi], M[yi], 'b.')
    filename = "%d_%s_%.3f_%s.png" % (i, dep_name, score, study_id)
    pp.save("%d_%s_%.3f_%s.png" % (i, dep_name, score, study_id))
    print filename

def main(jsonfile=None, k=500, outdir=""):
  J = json.load(open(jsonfile))
  M = np.ma.load(J["data_matrix"])
  varlist = [s.partition('\t')[0] for s in open(J["data_matrix"]) if s[0] not in ['#', '\n']]
  for dep_name, d in J['dependencies'].items():
    plot_dir = os.path.join(outdir, dep_name)
    make_dir(plot_dir)
    D = np.load(d["values_file"])
    try:
      A = np.load(d["argsorted_file"])
    except KeyError:
      print "no argsorted file. Sorting %s..." % d["values_file"]
      A = D.argsort()
    generate_top_k_scatters(M, D, A, varlist, k=k, dep_name=dep_name, study_id=J['gse_id'])
    
def make_dir(outdir):
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

    
if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)
