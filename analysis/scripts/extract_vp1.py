#!/usr/bin/env python3
import argparse, os
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

def fetch_record(acc, email, api_key=None):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(handle, 'genbank')
    handle.close()
    return rec

def find_vp1_coords(rec):
    for feat in rec.features:
        if feat.type == 'CDS' and 'gene' in feat.qualifiers:
            if any(g.lower()=='vp1' for g in feat.qualifiers['gene']):
                start = int(feat.location.start)
                end = int(feat.location.end)
                strand = int(feat.location.strand or 1)
                return start, end, strand
    raise SystemExit('VP1 not found in reference annotation')

def map_coords(ref_seq, q_seq, start, end):
    from Bio import pairwise2
    aln = pairwise2.align.globalms(str(ref_seq).upper(), str(q_seq).upper(), 2, -1, -5, -1, one_alignment_only=True)[0]
    ref_aln, q_aln = aln.seqA, aln.seqB
    ref_i = 0
    q_i = 0
    q_start = None
    q_end = None
    for ra, qa in zip(ref_aln, q_aln):
        if ra != '-':
            if ref_i == start and q_start is None:
                q_start = q_i
            ref_i += 1
            if ref_i == end:
                q_end = q_i
                break
        if qa != '-':
            q_i += 1
    if q_start is None or q_end is None:
        raise SystemExit('Failed to map VP1 coords to query')
    return q_start, q_end

def main():
    ap = argparse.ArgumentParser(description='Extract VP1 region from genomes by mapping to annotated reference')
    ap.add_argument('--email', required=False, help='Entrez email, or set NCBI_EMAIL')
    ap.add_argument('--api_key', required=False, help='NCBI API key, or set NCBI_API_KEY')
    ap.add_argument('--ref_acc', required=True, help='GenBank accession for annotated reference, e.g. D90457.1')
    ap.add_argument('--genomes', required=True, help='FASTA with one or more genomes')
    ap.add_argument('--out_fasta', required=True)
    args = ap.parse_args()

    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        raise SystemExit('Set --email or env NCBI_EMAIL')
    api_key = args.api_key or os.getenv('NCBI_API_KEY')

    ref = fetch_record(args.ref_acc, email, api_key)
    start, end, strand = find_vp1_coords(ref)
    ref_seq = ref.seq
    with open(args.out_fasta, 'w') as out:
        for rec in SeqIO.parse(args.genomes, 'fasta'):
            qs, qe = map_coords(ref_seq, rec.seq, start, end)
            frag = rec.seq[qs:qe]
            if strand == -1:
                frag = frag.reverse_complement()
            SeqIO.write(SeqIO.SeqRecord(frag, id=f'{rec.id}_VP1', description=''), out, 'fasta')

if __name__ == '__main__':
    main()
