import pandas as pd
import argparse

def filter_mirna_targets(input_file, output_file, max_expectation=1.0, inhibition='Cleavage',
                         min_multiplicity=2, max_upe=None, target_desc_keyword=None):
    """
    Фильтрует результаты предсказания мишеней miRNA по заданным параметрам.
    
    Параметры:
        input_file (str): Путь к входному файлу.
        output_file (str): Путь для сохранения отфильтрованных данных.
        max_expectation (float): Максимальное значение Expectation (по умолчанию 1.0).
        inhibition (str): Тип ингибирования (по умолчанию 'Cleavage').
        min_multiplicity (int): Минимальное значение Multiplicity (по умолчанию 2).
        max_upe (float): Максимальное значение UPE (None — не фильтровать).
        target_desc_keyword (str): Ключевое слово в описании гена (None — не фильтровать).
    """
    # Загрузка данных
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    # Базовая фильтрация (по умолчанию)
    df_filtered = df[
        (df['Expectation'] <= max_expectation) &
        (df['Inhibition'] == inhibition) &
        (df['Multiplicity'] >= min_multiplicity)
    ].copy()  # Добавлен .copy() для избежания предупреждений
    
    # Дополнительная фильтрация по UPE
    if max_upe is not None:
        df_filtered = df_filtered[df_filtered['UPE$'] <= max_upe]
    
    # Фильтрация по ключевому слову в описании гена
    if target_desc_keyword:
        df_filtered = df_filtered[
            df_filtered['Target_Desc.'].str.contains(target_desc_keyword, case=False, na=False)
        ]
    
    # Удаление дубликатов
    df_filtered = df_filtered.drop_duplicates(subset=['miRNA_Acc.', 'Target_Acc.'])
    
    # Сохранение результатов
    df_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"Результаты сохранены в {output_file}. Отфильтровано строк: {len(df_filtered)}")

def main():
    parser = argparse.ArgumentParser(description='Фильтрация результатов предсказания мишеней miRNA.')
    parser.add_argument('-i', '--input', required=True, help='Входной файл (формат psRNATarget)')
    parser.add_argument('-o', '--output', required=True, help='Выходной файл')
    parser.add_argument('--max_expectation', type=float, default=1.0, 
                        help='Максимальное значение Expectation (по умолчанию: 1.0)')
    # parser.add_argument('--inhibition', default='Cleavage', 
    #                     help='Тип ингибирования (по умолчанию: Cleavage)')
    parser.add_argument('--min_multiplicity', type=int, default=2, 
                        help='Минимальное значение Multiplicity (по умолчанию: 2)')
    parser.add_argument('--max_upe', type=float, 
                        help='Максимальное значение UPE (по умолчанию: не фильтровать)')
    parser.add_argument('--target_desc', 
                        help='Ключевое слово в описании гена (например, "transcription factor")')
    
    args = parser.parse_args()
    
    filter_mirna_targets(
        input_file=args.input,
        output_file=args.output,
        max_expectation=args.max_expectation,
        #inhibition=args.inhibition,
        min_multiplicity=args.min_multiplicity,
        max_upe=args.max_upe,
        target_desc_keyword=args.target_desc
    )

if __name__ == '__main__':
    main()