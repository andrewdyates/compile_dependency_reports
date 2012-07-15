#!/usr/bin/python
"""Plot histogram of numpy array.

python test_histogram.py np_matrix=~/gse2034/gse2034_dcor.values.npy plot_pdfname=~/test.pdf
"""
import sys
import random
import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

r_base = importr('base')
grdevices = importr('grDevices')
stats = importr('stats')
graphics = importr('graphics')

def main(np_matrix=None, plot_pdfname=None, k=70000, title="GSE2034 dCOR", xlab="dCOR"):
  assert np_matrix and plot_pdfname
  M = np.load(np_matrix)
  S = np.array(random.sample(M, k))
  grdevices.pdf(plot_pdfname)
  graphics.plot(stats.density(S), main=title, xlab=xlab)
  graphics.rug(S)
  grdevices.dev_off()

if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
  
  
