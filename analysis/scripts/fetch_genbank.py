#!/usr/bin/env python3
import argparse, time, os
from Bio import Entrez

ap = argparse.ArgumentParser(description='Fetch GenBank FASTA by accession list')
ap.add_argument('--email', required=False, help='Entrez email, or set env NCBI_EMAIL')
ap.add_argument('--api_key', required=False, help='NCBI API key, or set env NCBI_API_KEY')
ap.add_argument('--acc', required=True, help='Text file with one accession per line')
ap.add_argument('--out_fasta', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email:
    raise SystemExit('Set --email or env NCBI_EMAIL')
Entrez.email = email
api_key = a.api_key or os.getenv('NCBI_API_KEY')
if api_key:
    Entrez.api_key = api_key

accs = [x.strip() for x in open(a.acc) if x.strip() and not x.startswith('#')]
with open(a.out_fasta, 'w') as out:
    for acc in accs:
        h = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
        out.write(h.read())
        h.close()
        time.sleep(0.34)
print('Wrote', a.out_fasta)
