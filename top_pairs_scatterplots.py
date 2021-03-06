#!/usr/bin/python
"""Generate scatterplots of top k pairs.

EXAMPLE USE:

python top_pairs_scatterplots.py jsonfile=$HOME/gse2034/gse2034_best.json k=500 outdir=$HOME/gse2034/top_scatters
"""
from py_symmetric_matrix.cpy import *
from matplotlib import pyplot as pp
import numpy as np
import json
import os, sys, errno


def generate_top_k_scatters(M, D, A, varlist, dep_name, study_id, k=500, plot_dir=None):
  fp_log = open(os.path.join(plot_dir, "%s_%s_%d.log.txt" % (study_id, dep_name, k)), "w")
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
    filename = "%d_%s_%s_%s_%.3f_%s.png" % (i, x, y, dep_name, score, study_id)
    print filename
    filename = filename.replace('///','-').replace(' ', "_")
    if len(filename) > 100:
      filename = filename[:100] + ".png"
    pp.savefig(os.path.join(plot_dir, filename))
    fp_log.write("%s\t%s\t%d\t%.5f\n" % (x, y, i, score))
  fp_log.close()


def generate_top_k_named_scatters(M, D, A, varlist, dep_name, study_id, k=500, plot_dir=None, genelist=None):
  fp_log = open(os.path.join(plot_dir, "%s_%s_%d.log.txt" % (study_id, dep_name, k)), "w")
  i,rank = 0,0
  geneset = set(genelist)
  n = len(varlist)
  assert len(set(varlist) & geneset) == len(geneset)
  while i < k:
    if rank % 1000000 == 0:
      print "@ rank %d..." % rank
    if rank >= n*(n-1)//2:
      print "rank exhausted. Breaking..."
      break
    idx = A[-rank]
    score = D[idx]
    xi, yi = inv_sym_idx(idx, n)
    x, y = varlist[xi], varlist[yi]
    if x in geneset or y in geneset:
      print x,y
    rank += 1
    if (x not in geneset) and (y not in geneset):
      continue
    i += 1
    pp.clf(); pp.cla();
    pp.title("%d rank:%d %s %.3f" % (i, rank, dep_name, score))
    pp.xlabel(x); pp.ylabel(y)
    pp.plot(M[xi], M[yi], 'b.')
    filename = "%d_%d_%s_%s_%s_%.3f_%s.png" % (i, rank, x, y, dep_name, score, study_id)
    print filename
    filename = filename.replace('///','-').replace(' ', "_")
    if len(filename) > 100:
      filename = filename[:100] + ".png"
    pp.savefig(os.path.join(plot_dir, filename))
    fp_log.write("%s\t%s\t%d\t%.5f\n" % (x, y, i, score))
  fp_log.close()


def main(jsonfile=None, k=500, outdir="", genelist=None):
  k = int(k)
  J = json.load(open(jsonfile))
  if genelist is not None:
    if type(genelist) == str:
      genelist = genelist.split(',')
      
  M = np.ma.load(J["data_matrix"])
  varlist = [s.partition('\t')[0] for s in open(J["data_file"]) if s[0] not in ['#', '\n']]
  assert len(varlist) == J['n_vars'], "%d != %d" % (len(varlist), J['n_vars'])
  for dep_name, d in J['dependencies'].items():
    print "Plotting for %s..." % dep_name
    plot_dir = os.path.join(outdir, dep_name)
    make_dir(plot_dir)
    D = np.load(os.path.join(d["dir"], d["values_file"]))
    assert np.size(D) == J["n_pairs"]
    try:
      A = np.load(os.path.join(d["dir"], d["argsorted_file"]))
    except KeyError:
      print "no argsorted file. Sorting %s..." % d["values_file"]
      A = D.argsort()
      assert np.size(A) == J["n_pairs"]
    if not genelist:
      generate_top_k_scatters(M, D, A, varlist, k=k, dep_name=dep_name, study_id=J['gse_id'], plot_dir=plot_dir)
    else:
      print "Filtering for genes in %s..." % genelist
      generate_top_k_named_scatters(M, D, A, varlist, k=k, dep_name=dep_name, study_id=J['gse_id'], plot_dir=plot_dir, genelist=genelist)
    
def make_dir(outdir):
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

    
if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)
