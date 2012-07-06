#!/usr/bin/python
"""Given a set of dependency matrices for a microarray, generate a suite of analytics."""
import json
from numpy import *
from py_symmetric_matrix import *

def main(json_input=None, out_dir=None):
  
