"""Reverse-complement every record in a FASTA file."""
import sys
from Bio import SeqIO

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

out_seqs = list()
for record in SeqIO.parse(in_fasta, "fasta"):
    record.seq = record.seq.reverse_complement()
    out_seqs.append(record)

SeqIO.write(out_seqs, out_fasta, "fasta")
