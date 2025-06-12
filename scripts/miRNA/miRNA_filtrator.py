"""Script to filter miRNA candidates based on free energy (dG) and sequence length.

This script processes RNAfold output files, extracts miRNA candidates, and filters them
based on minimum free energy threshold and sequence length requirements.
"""

import argparse
import logging
import re
import sys
from typing import List, Optional, Tuple


def setup_logging(log_file: Optional[str] = None) -> None:
    """Configure logging for the script.

    Args:
        log_file: Path to log file. If None, logs to stderr.
    """
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    if log_file:
        logging.basicConfig(filename=log_file, level=logging.INFO, format=log_format)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format, stream=sys.stderr)


def parse_rnafold_output(filename: str) -> List[Tuple[str, str, float]]:
    """Parse RNAfold output file and extract miRNA candidates with their free energy values.

    Args:
        filename: Path to the RNAfold output file

    Returns:
        List of tuples containing (header, sequence, dG)

    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If the file format is invalid
    """
    try:
        with open(filename) as f:
            data = f.read().split(">")[1:]
    except FileNotFoundError:
        logging.error(f"Input file not found: {filename}")
        raise

    candidates: List[Tuple[str, str, float]] = []
    invalid_records = 0

    for record in data:
        lines = record.strip().split("\n")
        if len(lines) < 2:
            invalid_records += 1
            continue

        header: str = lines[0]
        seq: str = lines[1]

        # Find dG in line like ".... ( -3.10)"
        dG_match = re.search(r"\(([-\s]\d+\.\d+)\)", lines[-1])
        if not dG_match:
            logging.warning(f"No dG value found in line: {lines[-1]}")
            invalid_records += 1
            continue

        try:
            dG: float = float(dG_match.group(1).strip())
            candidates.append((header, seq, dG))
        except ValueError as e:
            logging.warning(f"Invalid dG value in line: {lines[-1]}. Error: {str(e)}")
            invalid_records += 1

    if invalid_records > 0:
        logging.info(f"Total invalid records skipped: {invalid_records}")

    return candidates


def main() -> None:
    """Parse command line arguments and filter miRNA candidates."""
    parser = argparse.ArgumentParser(
        description="Filter miRNA candidates based on dG and sequence length."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input file (RNAfold output)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output file (filtered candidates)"
    )
    parser.add_argument("--log", help="Path to log file (default: log to stderr)")
    args = parser.parse_args()

    # Setup logging
    setup_logging(args.log)

    try:
        candidates = parse_rnafold_output(args.input)

        with open(args.output, "w") as out:
            filtered_count = 0
            for header, seq, dG in candidates:
                # Filter candidates with dG < -20 and length >= 60
                if dG < -20 and len(seq) >= 60:
                    out.write(f">{header}\n{seq}\n")
                    filtered_count += 1

        logging.info(f"Successfully processed {args.input}:")
        logging.info(f"  - Total candidates: {len(candidates)}")
        logging.info(f"  - Filtered candidates: {filtered_count}")
        logging.info(f"  - Output saved to: {args.output}")

    except Exception as e:
        logging.error(f"Error processing files: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
