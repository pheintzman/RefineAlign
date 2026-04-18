#!/usr/bin/env python3

"""
Remove columns from a multiple sequence alignment where all sequences
have a gap ('-') or 'N'/'n', EXCEPT when the reference sequence has an
'N'/'n' at that position (to preserve reference coordinates).

Optionally remove any position that has a gap ('-') in the reference sequence.
"""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def remove_all_gap_and_N_columns(
    input_fasta,
    reference_id,
    output_fasta,
    remove_reference_gaps=False
):
    alignment = AlignIO.read(input_fasta, "fasta")
    seq_length = alignment.get_alignment_length()

    ref_record = next((r for r in alignment if r.id == reference_id), None)
    if ref_record is None:
        raise ValueError(f"Reference ID '{reference_id}' not found in alignment")

    # ------------------------------------------------------------------
    # Step 1: Reference-aware removal of all-gap / all-N columns
    # ------------------------------------------------------------------
    keep_indices = []

    for col_idx in range(seq_length):
        ref_base = ref_record.seq[col_idx]
        column = alignment[:, col_idx]

        # Always keep columns where reference has N/n
        if ref_base in ('N', 'n'):
            keep_indices.append(col_idx)
            continue

        # Otherwise, keep if column contains any real base
        if any(base not in ('-', 'N', 'n') for base in column):
            keep_indices.append(col_idx)

    cleaned_records = []
    for rec in alignment:
        new_seq = ''.join(rec.seq[i] for i in keep_indices)
        cleaned_records.append(
            SeqRecord(Seq(new_seq), id=rec.id, description="")
        )

    cleaned_alignment = MultipleSeqAlignment(cleaned_records)

    # ------------------------------------------------------------------
    # Step 2 (optional): Remove positions where reference has gaps
    # ------------------------------------------------------------------
    if remove_reference_gaps:
        ref_clean = next(r for r in cleaned_alignment if r.id == reference_id)
        ref_seq = str(ref_clean.seq)

        final_keep_indices = [i for i, base in enumerate(ref_seq) if base != "-"]

        if len(final_keep_indices) < len(ref_seq):
            trimmed_records = []
            for rec in cleaned_alignment:
                trimmed_seq = ''.join(rec.seq[i] for i in final_keep_indices)
                trimmed_records.append(
                    SeqRecord(Seq(trimmed_seq), id=rec.id, description="")
                )
            cleaned_alignment = MultipleSeqAlignment(trimmed_records)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    AlignIO.write(cleaned_alignment, output_fasta, "fasta")

    print(f"Cleaned alignment written to {output_fasta}")
    print(f"Original length: {seq_length}")
    print(f"Length after column filtering: {len(keep_indices)}")
    if remove_reference_gaps:
        print(f"Final length after removing reference gaps: {cleaned_alignment.get_alignment_length()}")
    else:
        print(f"Final length: {cleaned_alignment.get_alignment_length()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Remove all-gap/all-N columns from an alignment, preserving "
            "reference N positions, with optional removal of positions where "
            "the reference sequence has gaps."
        )
    )
    parser.add_argument("input_fasta")
    parser.add_argument("reference_id")
    parser.add_argument("output_fasta")
    parser.add_argument(
        "--remove-reference-gaps",
        action="store_true",
        help="Remove any column where the reference sequence has a '-'"
    )

    args = parser.parse_args()

    remove_all_gap_and_N_columns(
        args.input_fasta,
        args.reference_id,
        args.output_fasta,
        remove_reference_gaps=args.remove_reference_gaps
    )
