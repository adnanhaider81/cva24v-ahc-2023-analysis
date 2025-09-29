#!/usr/bin/env python3
import argparse, re, sys
ap = argparse.ArgumentParser(description='Select a CV-A24v reference from BLAST hits or fallback')
ap.add_argument('--blast_tsv', required=True)
ap.add_argument('--fallback', required=True)
ap.add_argument('--out_acc', required=True)
a = ap.parse_args()
acc = None
for line in open(a.blast_tsv):
    if not line.strip() or line.startswith('contig\t'):
        continue
    parts = line.strip().split('\t')
    sciname = parts[7].lower()
    sacc = parts[1]
    if 'coxsackievirus a24' in sciname or 'variant a24' in sciname:
        acc = sacc
        break
if acc is None:
    acc = a.fallback
open(a.out_acc, 'w').write(acc + '\n')
print('Selected', acc)
