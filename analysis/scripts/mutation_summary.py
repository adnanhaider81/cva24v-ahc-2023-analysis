#!/usr/bin/env python3
import argparse, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
ap = argparse.ArgumentParser(description='Summarize AA changes vs prototype and 2005 strains for ORF and VP1')
ap.add_argument('--consensus', required=True, help='combined consensus fasta in ref coords')
ap.add_argument('--prototype', required=True, help='prototype fasta D90457.1')
ap.add_argument('--pk2005', required=True, help='multi fasta with AB365074-AB365078')
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()

def translate(seq):
    seq = seq.replace('N','A')
    return str(Seq(seq).translate(to_stop=False))

ref = next(SeqIO.parse(a.prototype, 'fasta'))
samples = list(SeqIO.parse(a.consensus, 'fasta'))
pk2005 = list(SeqIO.parse(a.pk2005, 'fasta'))

with open(a.out_tsv, 'w') as out:
    out.write('sample\tpos\trefAA\tsampleAA\tpk2005_variants\n')
    refAA = translate(str(ref.seq))
    for s in samples:
        samAA = translate(str(s.seq))
        length = min(len(refAA), len(samAA))
        for i in range(length):
            if samAA[i] != refAA[i]:
                # presence in any 2005 variant
                pkset = set()
                for r in pk2005:
                    aa = translate(str(r.seq))[i] if i < len(translate(str(r.seq))) else '-'
                    if aa != refAA[i]:
                        pkset.add(aa)
                out.write(f'{s.id}\t{i+1}\t{refAA[i]}\t{samAA[i]}\t{",".join(sorted(pkset)) if pkset else ""}\n')
print('Wrote', a.out_tsv)
