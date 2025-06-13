"""Script to filter miRNA target predictions from psRNATarget results.

This script processes psRNATarget output files and filters miRNA-target pairs based on
various criteria such as expectation value, multiplicity, and target description.
"""

import pandas as pd
import argparse
import logging
import sys
from typing import Optional


def setup_logging(log_file: Optional[str] = None) -> None:
    """Configure logging for the script.
    
    Args:
        log_file: Path to log file. If None, logs to stderr.
    """
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format=log_format
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            stream=sys.stderr
        )


def filter_mirna_targets(
    input_file: str,
    output_file: str,
    max_expectation: float = 1.0,
    inhibition: str = 'Cleavage',
    min_multiplicity: int = 2,
    max_upe: Optional[float] = None,
    target_desc_keyword: Optional[str] = None
) -> None:
    """Filter miRNA target predictions based on specified parameters.
    
    Args:
        input_file: Path to input file
        output_file: Path for saving filtered data
        max_expectation: Maximum expectation value (default: 1.0)
        inhibition: Type of inhibition (default: 'Cleavage')
        min_multiplicity: Minimum multiplicity value (default: 2)
        max_upe: Maximum UPE value (None = no filtering)
        target_desc_keyword: Keyword in gene description (None = no filtering)
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If input data is invalid
    """
    try:
        # Load data
        df = pd.read_csv(input_file, sep='\t', comment='#')
        
        # Basic filtering
        df_filtered = df[
            (df['Expectation'] <= max_expectation) &
            (df['Inhibition'] == inhibition) &
            (df['Multiplicity'] >= min_multiplicity)
        ].copy()
        
        # Additional UPE filtering
        if max_upe is not None:
            df_filtered = df_filtered[df_filtered['UPE$'] <= max_upe]
        
        # Filter by keyword in gene description
        if target_desc_keyword:
            df_filtered = df_filtered[
                df_filtered['Target_Desc.'].str.contains(
                    target_desc_keyword,
                    case=False,
                    na=False
                )
            ]
        
        # Remove duplicates
        df_filtered = df_filtered.drop_duplicates(
            subset=['miRNA_Acc.', 'Target_Acc.']
        )
        
        # Save results
        df_filtered.to_csv(output_file, sep='\t', index=False)
        logging.info(
            f"Results saved to {output_file}. "
            f"Filtered rows: {len(df_filtered)}"
        )
        
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}")
        raise
    except Exception as e:
        logging.error(f"Error processing file: {str(e)}")
        raise


def main() -> None:
    """Parse command line arguments and run miRNA target filtering."""
    parser = argparse.ArgumentParser(
        description='Filter miRNA target predictions from psRNATarget results.'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file (psRNATarget format)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file'
    )
    parser.add_argument(
        '--max_expectation',
        type=float,
        default=1.0,
        help='Maximum expectation value (default: 1.0)'
    )
    parser.add_argument(
        '--min_multiplicity',
        type=int,
        default=2,
        help='Minimum multiplicity value (default: 2)'
    )
    parser.add_argument(
        '--max_upe',
        type=float,
        help='Maximum UPE value (default: no filtering)'
    )
    parser.add_argument(
        '--target_desc',
        help='Keyword in gene description (e.g., "transcription factor")'
    )
    parser.add_argument(
        '--log',
        help='Path to log file (default: log to stderr)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log)
    
    try:
        filter_mirna_targets(
            input_file=args.input,
            output_file=args.output,
            max_expectation=args.max_expectation,
            min_multiplicity=args.min_multiplicity,
            max_upe=args.max_upe,
            target_desc_keyword=args.target_desc
        )
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
