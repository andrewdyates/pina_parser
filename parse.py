#!/usr/bin/python
"""Parse PINA file into useful object using HUGO name corrections.

Sample line:
uniprotkb:Q14141	uniprotkb:Q96B97	uniprotkb:SEPT6(gene name)	uniprotkb:SH3KBP1(gene name)	-	-	MI:0004(affinity chromatography technology)	-	pubmed:19531213	taxid:9606(Homo sapiens)	taxid:9606(Homo sapiens)	MI:0915(physical association)	(biogrid)	BIOGRID_185461	-	prey	bait	go:GO:0005525|go:GO:0043679|go:GO:0007049|go:GO:0032154|go:GO:0000777|go:GO:0000910|go:GO:0005737|go:GO:0030496|go:GO:0031105|go:GO:0005819|go:GO:0008021|go:GO:0019048	go:GO:0006915|go:GO:0016477|go:GO:0007267|go:GO:0030659|go:GO:0005856|go:GO:0007010|go:GO:0005829|go:GO:0006897|go:GO:0007173|go:GO:0005925|go:GO:0042059|go:GO:0043005|go:GO:0005886|go:GO:0008360|go:GO:0045202	unspecified:32644
"""
import sys
import re
import hugo_gene_symbols

FNAME = "Homo sapiens-20121210.txt"
EXPECTED_HEADER = '"ID(s) interactor A"	"ID(s) interactor B"	"Alt. ID(s) interactor A"	"Alt. ID(s) interactor B"	"Alias(es) interactor A"	"Alias(es) interactor B"	"Interaction detection method(s)"	"Publication 1st author(s)"	"Publication Identifier(s)"	"Taxid interactor A"	"Taxid interactor B"	"Interaction type(s)"	"Source database(s)"	"Interaction identifier(s)"	"Confidence value(s)"	"Experimental role(s) interactor A"	"Experimental role(s) interactor B"	"Properties interactor A"	"Properties interactor B"	"HostOrganism(s)"'

RX_INTERACT = re.compile('[^(]+\(([^)]+)\)')
RX_UNIPROT = re.compile('uniprotkb:(.+)')
RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")
FNAME_OUT = "pina_compiled_may22_2013.tab"
HUMAN = "taxid:9606(Homo sapiens)"

def write_ppi(ppi):
  fp = open(FNAME_OUT, "w")
  print >>fp, "Gene A\tGene B\tInteraction Methods"
  for pair in sorted(ppi):
    a,b = pair.split(';')
    print >>fp, "%s\t%s\t%s" % (a,b,"|".join(ppi[pair]))
  fp.close()

def get_interact(s):
  r = set()
  for ss in s.split('|'):
    m = RX_INTERACT.match(s)
    if m:
      r.add(m.group(1))
  return r

def get_sym(H, prot_s, sym_s, alias_s, allow_miss=True):
  """Return official HUGO gene symbol."""
  prot_m = RX_UNIPROT.match(prot_s)
  if prot_m:
    prot = prot_m.group(1)
  else:
    prot=None
    
  sym_m = RX_GENE_NAME.match(sym_s)
  if sym_m:
    sym = sym_m.group(1)
  else:
    sym = None
    
  if alias_s in ("-","None"):
    alias = None
  else:
    alias = alias_s
  
  if prot in H.uniprot:
    return H.uniprot[prot]
  elif sym in H.official:
    return sym
  elif alias in H.official:
    return alias
  else:
    s = H.find_sym(sym)
    if s:
      return s
    else:
      ss = H.find_sym(alias)
      if ss:
        return ss
      else:
        # go ahead and return uniprot gene symbol that is not official HUGO gene symbol
        if allow_miss and sym not in ("-","None"):
          return sym
        else:
          return None

def main():
  fp = open(FNAME)
  header = fp.next().strip('\n\r')
  n_not_human = 0
  problems = {}
  ppi = {}
  syms = set()
  
  assert header == EXPECTED_HEADER
  H = hugo_gene_symbols.load()
  # protein ID to gene id
  for i, line in enumerate(fp):
    linenum = i+2
    row = line.strip('\n\r').split('\t')
    prot1, prot2 = row[0], row[1]
    sym1, sym2 = row[2], row[3]
    alias1, alias2 = row[4], row[5]
    tax1, tax2 = row[9], row[10]
    if tax1 != tax2 or tax1 != HUMAN:
      n_not_human += 1
      continue
    
    gene_a = get_sym(H, prot1, sym1, alias1)
    gene_b = get_sym(H, prot2, sym2, alias2)
    interact = get_interact(row[6])
    
    if gene_a: syms.add(gene_a)
    if gene_b: syms.add(gene_b)
    if not gene_a:
      problems.setdefault("%s;%s;%s"%(prot1, sym1, alias1),set()).add(linenum)
    elif not gene_b:
      problems.setdefault("%s;%s;%s"%(prot2, sym2, alias2),set()).add(linenum)
    else:
      key = ";".join(sorted((gene_a,gene_b)))
      ppi.setdefault(key,set()).update(interact)

  print "STATS"
  print "Lines read:", linenum-1
  print "#PPI:", len(ppi)
  print "Lines not human:", n_not_human
  print "#Unique Syms:", len(syms)
  print "#Problem Syms:", len(problems)
  print "#Problem PPI:", sum(len(s) for s in problems.values())
  print "Problem syms"
  for i, v in enumerate(problems.items()):
    print v

  print "Writing output file to %s..." % FNAME_OUT
  write_ppi(ppi)
  

if __name__ == "__main__":
  main()
