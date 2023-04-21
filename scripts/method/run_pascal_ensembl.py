import argparse
from PascalX import genescorer

parser = argparse.ArgumentParser(description='Arguments to Pascal scoring.')
parser.add_argument("-r", "--refpanel", help = "Reference panel path")
parser.add_argument("-g", "--gwas", help = "GWAS summary statistics path")
parser.add_argument("-a", "--annotation", help = "Genome annotation file")
parser.add_argument("-o", "--output", help = "Output file path")

args = parser.parse_args()

Scorer = genescorer.chi2sum()
Scorer.load_refpanel(args.refpanel, keepfile=None, parallel=4)
print("Reference panel loaded.")

Scorer.load_genome(args.annotation, useNAgenes=True, csymb=0)
print("Genome annotation loaded.")

Scorer.load_GWAS(args.gwas,rscol=0,pcol=6,header=True)
print("GWAS dataset loaded.")

print("Starting scoring...")
R = Scorer.score_all(parallel=4)
Scorer.save_scores(args.output)