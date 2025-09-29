#!/usr/bin/env python3
import argparse, csv, re, yaml
ap = argparse.ArgumentParser(description='Merge Kraken2 and Kaiju reports and flag likely targets')
ap.add_argument('--k2_reads', required=False)
ap.add_argument('--k2_contigs', required=False)
ap.add_argument('--kaiju', required=False)
ap.add_argument('--out_tsv', required=True)
ap.add_argument('--out_targets', required=True)
a = ap.parse_args()

def parse_kraken_report(path):
    rows = []
    if not path: return rows
    with open(path) as f:
        for line in f:
            if not line.strip(): 
                continue
            # Kraken2 report: pct  reads  taxa... name
            parts = re.split(r'\s{2,}', line.strip())
            try:
                pct = float(parts[0])
            except Exception:
                continue
            name = parts[-1]
            rows.append(('kraken2', name, pct))
    return rows

def parse_kaiju_table(path):
    rows = []
    if not path: return rows
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        hdr = next(reader, None)
        for r in reader:
            if len(r) < 4:
                continue
            name = r[1]
            try:
                pct = float(r[3])
            except Exception:
                pct = 0.0
            rows.append(('kaiju', name, pct))
    return rows

rows = []
rows += parse_kraken_report(a.k2_reads)
rows += parse_kraken_report(a.k2_contigs)
rows += parse_kaiju_table(a.kaiju)

# Sum pct by name
from collections import defaultdict
agg = defaultdict(lambda: [0.0, set()])
for src, name, pct in rows:
    agg[name][0] += pct
    agg[name][1].add(src)

ordered = sorted(agg.items(), key=lambda x: x[1][0], reverse=True)

with open(a.out_tsv, 'w') as out:
    out.write('name\tcombined_pct\tsources\n')
    for name, (pct, srcs) in ordered:
        out.write(f'{name}\t{pct:.3f}\t{",".join(sorted(srcs))}\n')

# Flag targets
targets = []
names_norm = [n.lower() for n,_ in ordered]
def present_any(subs):
    return any(any(s in n for s in subs) for n in names_norm)

if present_any(['coxsackievirus a24', 'cv-a24v', 'enterovirus a24']):
    targets.append('cva24v')
if present_any(['enterovirus 70', 'ev70']):
    targets.append('ev70')
if present_any(['adenovirus']):
    targets.append('adenovirus')

yaml.safe_dump({'selected_targets': targets}, open(a.out_targets,'w'))
print('Wrote', a.out_tsv, 'and', a.out_targets)
