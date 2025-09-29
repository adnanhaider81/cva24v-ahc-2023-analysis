#!/usr/bin/env python3
import argparse, subprocess, csv, sys
ap = argparse.ArgumentParser(description='BLAST contigs against nt and extract top hits')
ap.add_argument('--contigs', required=True)
ap.add_argument('--out_tsv', required=True)
ap.add_argument('--remote', action='store_true')
ap.add_argument('--db', default='nt')
ap.add_argument('--max_hits', type=int, default=5)
a = ap.parse_args()
cmd = ['blastn', '-query', a.contigs, '-outfmt', '6 qseqid sacc pident length evalue bitscore staxids sscinames sskingdoms']
cmd += ['-remote', '-db', 'nt'] if a.remote else ['-db', a.db]
res = subprocess.run(cmd, capture_output=True, text=True, check=True)
lines = [l for l in res.stdout.strip().splitlines() if l.strip()]
rows = {}
for l in lines:
    parts = l.split('\t')
    qid = parts[0]
    rows.setdefault(qid, []).append(parts)
with open(a.out_tsv, 'w') as out:
    out.write('contig\tacc\tpident\tlength\tevalue\tbitscore\ttaxid\tsciname\tskingdom\n')
    for qid, hits in rows.items():
        for h in hits[:a.max_hits]:
            out.write('\t'.join([qid] + h[1:]) + '\n')
print('Wrote', a.out_tsv)
