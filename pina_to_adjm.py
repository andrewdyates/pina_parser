#!/usr/bin/python
"""Convert edge list file to adj matrix.

WARNING: adjm is a large file. Move it out of the source code respository on creation.
"""
import numpy as np
import matrix_io as mio
import sys

FNAME = "pina_compiled_may22_2013.tab"
FNAME_OUT = "pina_compiled_may22_2013_adjm.tab"

def main():
  fp = open(FNAME)
  fp.next()
  ppi = {}
  genes = set()
  for line in fp:
    row = line.strip('\n\r').split('\t')
    genes.add(row[0]); genes.add(row[1])
    ppi.setdefault(row[0],set()).add(row[1])

  rownames = sorted(genes)
  colnames = rownames
  row_idx = dict(((s,i) for i,s in enumerate(rownames)))
  A = np.zeros((len(rownames), len(colnames)))

  for gene_i in ppi:
    for gene_j in ppi[gene_i]:
      i = row_idx[gene_i]
      j = row_idx[gene_j]
      A[i,j] = 1; A[j,i] = 1

  print >>sys.stderr, "# interactions x2: %d" % np.sum(A)
  print >>sys.stderr, "# (genes x genes) size:", A.shape

  mio.save(A, fp=FNAME_OUT, row_ids=['"%s"'%s for s in rownames], col_ids=['"%s"'%s for s in colnames], fmt="%d")

if __name__ == "__main__":
  main()
