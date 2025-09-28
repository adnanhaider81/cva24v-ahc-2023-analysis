#!/usr/bin/env python3
import argparse, os
from Bio import Entrez, SeqIO

def fetch_record(acc, email, api_key=None):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(handle, 'genbank')
    handle.close()
    return rec

def get_orf_coords(rec):
    for feat in rec.features:
        if feat.type == 'CDS':
            if 'product' in feat.qualifiers and any('polyprotein' in p.lower() for p in feat.qualifiers['product']):
                start = int(feat.location.start)
                end = int(feat.location.end)
                strand = int(feat.location.strand or 1)
                return start, end, strand
    for feat in rec.features:
        if feat.type == 'CDS':
            start = int(feat.location.start)
            end = int(feat.location.end)
            strand = int(feat.location.strand or 1)
            return start, end, strand
    raise SystemExit('No CDS found in reference')

def main():
    ap = argparse.ArgumentParser(description='Translate ORF and compare AA to prototype reference')
    ap.add_argument('--email', required=False, help='Entrez email, or set NCBI_EMAIL')
    ap.add_argument('--api_key', required=False, help='NCBI API key, or set NCBI_API_KEY')
    ap.add_argument('--ref_acc', default='D90457.1')
    ap.add_argument('--genomes', required=True)
    ap.add_argument('--out_tsv', required=True)
    args = ap.parse_args()

    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        raise SystemExit('Set --email or env NCBI_EMAIL')
    api_key = args.api_key or os.getenv('NCBI_API_KEY')

    ref = fetch_record(args.ref_acc, email, api_key)
    rs, re, rstrand = get_orf_coords(ref)
    rnt = ref.seq[rs:re] if rstrand == 1 else ref.seq[rs:re].reverse_complement()
    raa = rnt.translate()

    with open(args.out_tsv, 'w') as out:
        out.write('sample	pos	ref_aa	alt_aa
')
        for rec in SeqIO.parse(args.genomes, 'fasta'):
            qaa = rec.seq.translate()
            L = min(len(raa), len(qaa))
            for i in range(L):
                if raa[i] != qaa[i]:
                    out.write(f'{rec.id}	{i+1}	{raa[i]}	{qaa[i]}
')

if __name__ == '__main__':
    main()
