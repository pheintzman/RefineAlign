#!/usr/bin/env python3

import argparse
import sys
import shutil
from pathlib import Path


def read_fasta(filename):
    """Read a FASTA file and return headers and sequences."""

    headers = []
    sequences = []
    current_sequence = []

    with open(filename, 'r') as infile:
        for line in infile:
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                if headers:
                    sequences.append(''.join(current_sequence))

                headers.append(line)
                current_sequence = []

            else:
                if not headers:
                    raise ValueError("FASTA sequence found before first header")

                current_sequence.append(line)

    if headers:
        sequences.append(''.join(current_sequence))

    if not headers:
        raise ValueError("No sequences found in FASTA file")

    return headers, sequences


def parse_regions(region_string):
    """Parse comma-separated, 1-based inclusive coordinate ranges."""

    regions = []

    for region in region_string.split(','):
        try:
            start, end = map(int, region.strip().split('-'))
        except ValueError:
            raise ValueError(
                "Invalid region '{}'. Expected format: start-end".format(region)
            )

        if start < 1 or end < start:
            raise ValueError(
                "Invalid coordinates '{}'. Coordinates must be positive "
                "and end must be >= start.".format(region)
            )

        regions.append((start, end))

    return regions


def mask_regions(sequences, reference_index, regions):
    """Mask specified reference-coordinate regions in all non-reference sequences."""

    reference = sequences[reference_index]

    # Check that all sequences have the same alignment length
    alignment_lengths = set(len(seq) for seq in sequences)

    if len(alignment_lengths) != 1:
        raise ValueError("Sequences have different alignment lengths")

    # Identify alignment columns corresponding to reference coordinates
    problematic_columns = set()

    reference_position = 0

    for column, base in enumerate(reference):

        # Gaps do not count towards reference coordinates
        if base == '-':
            continue

        reference_position += 1

        # Check whether this reference position belongs to a problematic region
        for start, end in regions:
            if start <= reference_position <= end:
                problematic_columns.add(column)
                break

    # Check that all specified coordinates exist in the reference
    reference_length = reference_position

    for start, end in regions:
        if end > reference_length:
            raise ValueError(
                "Region {}-{} exceeds reference length ({})".format(
                    start, end, reference_length
                )
            )

    # Mask problematic columns in all sequences except the reference
    masked_sequences = []

    for i, seq in enumerate(sequences):

        # Preserve the reference sequence unchanged
        if i == reference_index:
            masked_sequences.append(seq)
            continue

        seq = list(seq)

        for column in problematic_columns:

            # Preserve alignment gaps; convert all other characters to N
            if seq[column] != '-':
                seq[column] = 'N'

        masked_sequences.append(''.join(seq))

    return masked_sequences, len(problematic_columns)


def main():

    parser = argparse.ArgumentParser(
        description="Mask problematic alignment regions using reference coordinates."
    )

    parser.add_argument(
        "fastafile",
        help="Input aligned FASTA file"
    )

    parser.add_argument(
        "reference_sequence",
        help="ID of the reference sequence (FASTA header, without >)"
    )

    parser.add_argument(
        "--problematic_regions",
        default="no_problematic_regions",
        help="Comma-separated reference-coordinate ranges (e.g. 15-25,30-34), "
             "or no_problematic_regions to skip masking."
    )

    args = parser.parse_args()

    input_path = Path(args.fastafile)
    output_path = input_path.with_name(
        input_path.stem + ".problematic_regions_removed" + input_path.suffix
    )

    if args.problematic_regions == "no_problematic_regions":
        print("No problematic regions specified; copying alignment unchanged.",
              file=sys.stderr)
        shutil.copyfile(input_path, output_path)
        print(f"Output written to: {output_path}", file=sys.stderr)
        sys.exit(0)
    
    try:

        # Read alignment
        headers, sequences = read_fasta(args.fastafile)

        # Identify reference sequence
        reference_ids = [
            header[1:].split()[0] for header in headers
        ]

        if args.reference_sequence not in reference_ids:
            raise ValueError(
                "Reference sequence '{}' not found in FASTA file".format(
                    args.reference_sequence
                )
            )

        reference_index = reference_ids.index(args.reference_sequence)

        # Parse problematic regions
        regions = parse_regions(args.problematic_regions)

        # Mask problematic regions
        masked_sequences, masked_columns = mask_regions(
            sequences,
            reference_index,
            regions
        )

        # Generate output filename
        input_path = Path(args.fastafile)

        output_path = input_path.with_name(
            input_path.stem + ".problematic_regions_removed" + input_path.suffix
        )

        # Write output FASTA
        with open(output_path, 'w') as outfile:

            for header, seq in zip(headers, masked_sequences):
                outfile.write(header + '\n')
                outfile.write(seq + '\n')

        print("Reference sequence:", args.reference_sequence, file=sys.stderr)
        print("Problematic regions:", args.problematic_regions, file=sys.stderr)
        print("Masked alignment columns:", masked_columns, file=sys.stderr)
        print("Output file:", output_path, file=sys.stderr)

    except (ValueError, OSError) as error:
        print("ERROR:", error, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
