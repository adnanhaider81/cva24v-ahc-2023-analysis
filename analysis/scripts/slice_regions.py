#!/usr/bin/env python3
import argparse, os
from Bio import Entrez, SeqIO
from Bio import pairwise2
from pathlib import Path

REGIONS = {
  '5UTR': ('five_prime_UTR',),
  '2A': ('2A',),
  '2B': ('2B',),
  '2C': ('2C',),
  '3A': ('3A',),
  '3B': ('3B',),
  '3C': ('3C',),
  '3D': ('3D',)
}

def fetch_record(acc, email, api_key=None):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(handle, 'genbank')
    handle.close()
    return rec

def find_region(rec, names):
    for feat in rec.features:
        if feat.type in ['CDS', 'gene', 'misc_feature']:
            q = []
            for key in ['gene','product','note']:
                if key in feat.qualifiers:
                    q += [x.lower() for x in feat.qualifiers[key]]
            if any(any(n.lower() in s for s in q) for n in names):
                return int(feat.location.start), int(feat.location.end), int(feat.location.strand or 1)
    return None

def slice_and_write(ref_seq, q_seq, s, e, strand, out_path, qid):
    aln = pairwise2.align.globalms(str(ref_seq), str(q_seq), 2, -1, -5, -1, one_alignment_only=True)[0]
    ref_aln, q_aln = aln.seqA, aln.seqB
    ref_i = 0; q_i = 0; q_s = None; q_e = None
    for ra, qa in zip(ref_aln, q_aln):
        if ra != '-':
            if ref_i == s and q_s is None:
                q_s = q_i
            ref_i += 1
            if ref_i == e:
                q_e = q_i
                break
        if qa != '-':
            q_i += 1
    if q_s is None or q_e is None:
        return False
    frag = q_seq[q_s:q_e]
    if strand == -1:
        frag = frag.reverse_complement()
    with open(out_path, 'w') as out:
        out.write(f">{qid}
{frag}
")
    return True

def main():
    ap = argparse.ArgumentParser(description='Slice named genome regions from query genomes by mapping to annotated reference')
    ap.add_argument('--email', required=False, help='Entrez email, or set NCBI_EMAIL')
    ap.add_argument('--api_key', required=False, help='NCBI API key, or set NCBI_API_KEY')
    ap.add_argument('--ref_acc', default='D90457.1')
    ap.add_argument('--genomes', required=True)
    ap.add_argument('--out_dir', required=True)
    args = ap.parse_args()

    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        raise SystemExit('Set --email or env NCBI_EMAIL')
    api_key = args.api_key or os.getenv('NCBI_API_KEY')

    ref = fetch_record(args.ref_acc, email, api_key)
    ref_seq = ref.seq

    coords = {}
    for key, names in REGIONS.items():
        hit = find_region(ref, names)
        if hit:
            coords[key] = hit

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    for rec in SeqIO.parse(args.genomes, 'fasta'):
        for key, (s, e, strand) in coords.items():
            out_path = Path(args.out_dir)/f"{rec.id}_{key}.fasta"
            slice_and_write(ref_seq, rec.seq, s, e, strand, out_path, f"{rec.id}_{key}")

if __name__ == '__main__':
    main()
