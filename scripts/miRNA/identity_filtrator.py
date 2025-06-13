"""Script to filter miRNA BLAST results by identity threshold and extract unique miRNA names.

This script processes BLAST results in format 6, filters hits based on identity threshold,
and generates two output files: filtered BLAST results and a list of unique miRNA identifiers.
"""

import argparse
import logging
import os
import sys


def setup_logging(log_file=None):
    """Configure logging for the script.

    Args:
        log_file (str, optional): Path to log file. If None, logs to stderr.
    """
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    if log_file:
        logging.basicConfig(filename=log_file, level=logging.INFO, format=log_format)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format, stream=sys.stderr)


def filter_miRNA_hits(
    input_file, output_filtered, output_unique, identity_threshold=90.0, id_column=1
):
    """Filter BLAST results by identity and extract unique miRNA names.

    Args:
        input_file (str): Path to input BLAST results file
        output_filtered (str): Path to output file for filtered results
        output_unique (str): Path to output file for unique identifiers
        identity_threshold (float): Minimum percent identity threshold (default: 90.0)
        id_column (int): Column number containing miRNA IDs (default: 1)
    """
    try:
        with open(input_file, "r") as infile, open(
            output_filtered, "w"
        ) as out_filtered, open(output_unique, "w") as out_unique:

            unique_names = set()
            invalid_lines = 0

            for line_num, line in enumerate(infile, 1):
                columns = line.strip().split("\t")
                if len(columns) > 2:
                    try:
                        identity = float(columns[2])
                        if identity >= identity_threshold:
                            out_filtered.write(line)
                            if len(columns) > id_column:
                                unique_names.add(columns[id_column])
                    except ValueError:
                        invalid_lines += 1
                        logging.warning(
                            f"Invalid identity value in line {line_num}: '{columns[2]}'. "
                            f"Expected a float number."
                        )
                        continue

            if invalid_lines > 0:
                logging.info(f"Total invalid lines skipped: {invalid_lines}")

            # Write unique identifiers
            for name in sorted(unique_names):
                out_unique.write(f"{name}\n")

        logging.info(f"Successfully processed {input_file}:")
        logging.info(f"  - Filtered hits saved to {output_filtered}")
        logging.info(f"  - Unique IDs saved to {output_unique}")

    except Exception as e:
        logging.error(f"Error processing file {input_file}: {str(e)}")
        sys.exit(1)


def main():
    """Parse command line arguments and run the miRNA filtering process."""
    parser = argparse.ArgumentParser(
        description="Filter miRNA hits by identity threshold and extract unique miRNA names",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input BLAST results file (format 6)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Base name for output files (default: derived from input)",
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=90.0,
        help="Minimum percent identity threshold",
    )
    parser.add_argument(
        "--id-column",
        type=int,
        default=1,
        help="Column number containing miRNA IDs (0-based)",
    )
    parser.add_argument(
        "--no-filtered", action="store_true", help="Skip saving filtered BLAST results"
    )
    parser.add_argument(
        "--log", default=None, help="Path to log file (default: log to stderr)"
    )

    args = parser.parse_args()

    # Setup logging
    setup_logging(args.log)

    # Determine output filenames
    base_name = args.output if args.output else os.path.splitext(args.input)[0]
    filtered_file = f"{base_name}_filtered.txt" if not args.no_filtered else os.devnull
    unique_file = f"{base_name}_unique_list.txt"

    # Run filtering
    filter_miRNA_hits(
        input_file=args.input,
        output_filtered=filtered_file,
        output_unique=unique_file,
        identity_threshold=args.identity,
        id_column=args.id_column,
    )


if __name__ == "__main__":
    main()
