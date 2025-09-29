#!/usr/bin/env python3
import argparse, os
from Bio import Entrez, SeqIO
ap = argparse.ArgumentParser(description='Fetch single GenBank FASTA and optional GenBank flat file')
ap.add_argument('--acc', required=True)
ap.add_argument('--out_fasta', required=True)
ap.add_argument('--out_gb', required=False)
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
a = ap.parse_args()
email = a.email or os.getenv('NCBI_EMAIL')
if not email:
    raise SystemExit('Set --email or env NCBI_EMAIL')
Entrez.email = email
api_key = a.api_key or os.getenv('NCBI_API_KEY')
if api_key:
    Entrez.api_key = api_key
r = Entrez.efetch(db='nucleotide', id=a.acc, rettype='fasta', retmode='text')
open(a.out_fasta, 'w').write(r.read()); r.close()
if a.out_gb:
    r = Entrez.efetch(db='nucleotide', id=a.acc, rettype='gb', retmode='text')
    open(a.out_gb, 'w').write(r.read()); r.close()
print('Fetched', a.acc)
