#!/usr/bin/env python3
import argparse, os, re
from Bio import SeqIO
ap = argparse.ArgumentParser(description='Extract VP1 region from consensus FASTA using VP1 coordinates from a GenBank flat file of the reference')
ap.add_argument('--ref_gb', required=True, help='GenBank flat file with VP1 annotated')
ap.add_argument('--consensus_fa', required=True, help='FASTA of consensus sequences in reference coordinates')
ap.add_argument('--out_vp1_fa', required=True)
a = ap.parse_args()

# Find VP1 coordinates
gb = next(SeqIO.parse(a.ref_gb, 'genbank'))
vp1_locs = []
for feat in gb.features:
    if feat.type == 'CDS':
        prod = str(feat.qualifiers.get('product', [''])[0]).lower()
        gene = str(feat.qualifiers.get('gene', [''])[0]).lower()
        if 'vp1' in prod or gene == 'vp1':
            vp1_locs.append(feat.location)
if not vp1_locs:
    raise SystemExit('No VP1 feature found in reference GenBank')
loc = vp1_locs[0]
start = int(loc.start) + 1
end = int(loc.end)
with open(a.out_vp1_fa, 'w') as out:
    for rec in SeqIO.parse(a.consensus_fa, 'fasta'):
        seq = str(rec.seq)
        vp1 = seq[start-1:end]
        out.write(f'>{rec.id}_VP1\n{vp1}\n')
print(f'Wrote {a.out_vp1_fa} using VP1 {start}-{end}')
