#!/usr/bin/env python3

"""
Rotate an MSA so that any leading gaps in the reference sequence
are moved to the end of the alignment, then remove all gaps.

This script ONLY handles gaps at the START of the reference.
"""

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def count_start_gaps(seq):
    """Return number of leading gaps in a sequence."""
    s = str(seq)
    i = 0
    while i < len(s) and s[i] == "-":
        i += 1
    return i


def rotate_alignment_left(alignment, shift):
    """Rotate alignment columns left by shift."""
    aln_len = alignment.get_alignment_length()
    shift = shift % aln_len

    rotated = []
    for rec in alignment:
        new_seq = rec.seq[shift:] + rec.seq[:shift]
        rotated.append(SeqRecord(new_seq, id=rec.id, description=""))

    return rotated


def remove_all_gaps(records):
    """Remove all gaps from all sequences."""
    return [
        SeqRecord(
            Seq(str(rec.seq).replace("-", "")),
            id=rec.id,
            description=""
        )
        for rec in records
    ]


def main(input_fasta, reference_id, output_fasta):
    alignment = AlignIO.read(input_fasta, "fasta")

    ref = next((r for r in alignment if r.id == reference_id), None)
    if ref is None:
        raise ValueError(f"Reference ID '{reference_id}' not found")

    start_gaps = count_start_gaps(ref.seq)

    records = alignment
    if start_gaps > 0:
        records = rotate_alignment_left(records, start_gaps)

    records = remove_all_gaps(records)

    SeqIO.write(records, output_fasta, "fasta")
    print(f"Output written to {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rotate MSA based on leading gaps in reference sequence"
    )
    parser.add_argument("input_fasta")
    parser.add_argument("reference_id")
    parser.add_argument("output_fasta")

    args = parser.parse_args()
    main(args.input_fasta, args.reference_id, args.output_fasta)
