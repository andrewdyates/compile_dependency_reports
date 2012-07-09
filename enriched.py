from util import *
from py_symmetric_matrix import *
import random


class Enriched(object):

  def __init__(self, pairs):
    """Initialize.
    Args:
      pairs [*(str, str)] iterable of pairs of names.
    """
    # TODO: implement pair set as adjancency list
    self.genes = set()
    self.pairs = []
    self.pairs_hash = set()
    
    for x, y in pairs:
      x, y = clean(x), clean(y)
      x, y = sorted((x,y))
      if not (x and y) or x == y:
        continue
      q = self._hash(x,y)
      if not q in self.pairs_hash:
        self.genes.add(x)
        self.genes.add(y)
        self.pairs.append((x,y))
        self.pairs_hash.add(q)

  def _hash(self, x,y):
    return "%s,%s" % (x,y)
      
  def exists(self, x, y):
    x, y = sorted((clean(x),clean(y)))
    if not (x and y):
      return False
    if not (x in self.genes and y in self.genes):
      return False
    return self._hash(x,y) in self.pairs_hash

  def indices(self, varlist):
    """Return indices of self.pairs for sym dependency matrix indexed by varlist.
    """
    n = len(varlist)
    a = []
    cleaned_varlist = [clean(v) for v in varlist]
    var_d = dict([(v, i) for i, v in enumerate(cleaned_varlist)])
    for x, y in self.pairs:
      try:
        i, j = var_d[x], var_d[y]
      except KeyError:
        continue
      idx = sym_idx(i,j,n)
      a.append(idx)
    return a
    


class PINAEnriched(Enriched):

  RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")
  NAME_COL_A = 2
  NAME_COL_B = 3
  
  def _pair_gen(self, fp):
    for line in fp:
      row = line[:-1].split('\t')
      s_a, s_b = (row[self.NAME_COL_A], row[self.NAME_COL_B])
      m_a, m_b = self.RX_GENE_NAME.match(s_a), self.RX_GENE_NAME.match(s_b)
      if not (m_b and m_b):
        continue
      yield (m_a.group(1), m_b.group(1))
      
  def __init__(self, fp):
    """Initialize from new open fp to PINA formatted file."""
    super(PINAEnriched, self).__init__(pairs=self._pair_gen(fp))


    
class RandomPINAEnriched(PINAEnriched):
  """Generate randomized enrichment set from a PINA set.

  Rather hacky...
  """
  def __init__(self, fp):
    super(RandomPINAEnriched, self).__init__(fp)
    # replace self.pairs and self.pairs_hash with randomly generated sets of equal length
    n_pairs = len(self.pairs_hash)
    n = len(self.genes)
    gene_list = list(self.genes)
    self.pairs = []
    self.pairs_hash = set()
    
    idxs = random.sample(xrange(int(n*(n-1)/2)), n_pairs)
    for i in idxs:
      x, y = inv_sym_idx(i, n)
      gene_x, gene_y = sorted((gene_list[x], gene_list[y]))
      self.pairs.append((gene_x, gene_y))
      self.pairs_hash.add("%s,%s"%(gene_x, gene_y))
