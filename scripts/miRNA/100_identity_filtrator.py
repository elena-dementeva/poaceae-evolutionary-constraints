import argparse
import os
import sys

def filter_miRNA_hits(input_file, output_filtered, output_unique, identity_threshold=90.0, id_column=1):
    """
    Фильтрует результаты BLAST по идентичности и извлекает уникальные имена miRNA
    :param input_file: входной файл с результатами BLAST
    :param output_filtered: файл для отфильтрованных результатов
    :param output_unique: файл для уникальных идентификаторов
    :param identity_threshold: порог идентичности (по умолчанию 90%)
    :param id_column: номер колонки с идентификаторами (по умолчанию 1)
    """
    try:
        with open(input_file, 'r') as infile, \
             open(output_filtered, 'w') as out_filtered, \
             open(output_unique, 'w') as out_unique:
            
            unique_names = set()
            
            for line in infile:
                columns = line.strip().split('\t')
                if len(columns) > 2:
                    try:
                        identity = float(columns[2])
                        if identity >= identity_threshold:
                            out_filtered.write(line)
                            if len(columns) > id_column:
                                unique_names.add(columns[id_column])
                    except ValueError:
                        continue  # Пропускаем строки с некорректными значениями
            
            # Записываем уникальные идентификаторы
            for name in sorted(unique_names):
                out_unique.write(f"{name}\n")
                
        print(f"Successfully processed {input_file}:")
        print(f"  - Filtered hits saved to {output_filtered}")
        print(f"  - Unique IDs saved to {output_unique}")
        
    except Exception as e:
        print(f"Error processing file {input_file}: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Filter miRNA hits by identity threshold and extract unique miRNA names',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, 
                       help='Input BLAST results file (format 6)')
    parser.add_argument('-o', '--output', default=None,
                       help='Base name for output files (default: derived from input)')
    parser.add_argument('--identity', type=float, default=90.0,
                       help='Minimum percent identity threshold')
    parser.add_argument('--id-column', type=int, default=1,
                       help='Column number containing miRNA IDs (0-based)')
    parser.add_argument('--no-filtered', action='store_true',
                       help='Skip saving filtered BLAST results')
    
    args = parser.parse_args()
    
    # Определяем имена выходных файлов
    base_name = args.output if args.output else os.path.splitext(args.input)[0]
    filtered_file = f"{base_name}_filtered.txt" if not args.no_filtered else os.devnull
    unique_file = f"{base_name}_unique_list.txt"
    
    # Запускаем фильтрацию
    filter_miRNA_hits(
        input_file=args.input,
        output_filtered=filtered_file,
        output_unique=unique_file,
        identity_threshold=args.identity,
        id_column=args.id_column
    )

if __name__ == '__main__':
    main()