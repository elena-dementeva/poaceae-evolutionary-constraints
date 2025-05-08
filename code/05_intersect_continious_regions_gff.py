#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

def parse_gff(gff_file):
    """Парсит GFF файл и возвращает DataFrame с features."""
    features = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature_type, start, end, score, strand, frame, attributes = fields
            
            # Парсинг атрибутов
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key.strip()] = value.strip()
            
            # Извлечение ID признака (если есть)
            feature_id = attr_dict.get('ID', '')
            
            features.append({
                'chromosome': chrom,
                'start': int(start),
                'end': int(end),
                'feature_type': feature_type,
                'strand': strand,
                'feature_id': feature_id,  # Добавляем ID признака
                'attributes': attr_dict
            })
    
    return pd.DataFrame(features)

def find_overlaps(constrained_df, features_df):
    """Находит пересечения между constrained регионами и features."""
    overlaps = []
    
    # Группируем features по хромосомам для ускорения поиска
    features_by_chrom = features_df.groupby('chromosome')
    
    for idx, region in constrained_df.iterrows():
        chrom = region['Chromosome']
        reg_start = region['Start']
        reg_end = region['End']
        
        if chrom in features_by_chrom.groups:
            chrom_features = features_by_chrom.get_group(chrom)
            
            # Находим пересечения
            mask = ((chrom_features['start'] <= reg_end) & 
                    (chrom_features['end'] >= reg_start))
            
            overlapping_features = chrom_features[mask]
            
            for _, feature in overlapping_features.iterrows():
                # Вычисляем координаты пересечения
                overlap_start = max(reg_start, feature['start'])
                overlap_end = min(reg_end, feature['end'])
                overlap_length = overlap_end - overlap_start + 1
                
                overlaps.append({
                    'Chromosome': chrom,
                    'Constrained_Start': reg_start,
                    'Constrained_End': reg_end,
                    'Constrained_Length': region['Length'],
                    'Status': region['Status'],
                    'Feature_Type': feature['feature_type'],
                    'Feature_ID': feature['feature_id'],  # Добавляем ID признака
                    'Feature_Start': feature['start'],
                    'Feature_End': feature['end'],
                    'Strand': feature['strand'],
                    'Overlap_Start': overlap_start,
                    'Overlap_End': overlap_end,
                    'Overlap_Length': overlap_length,
                    'Overlap_Percent': (overlap_length / region['Length']) * 100
                })
    
    return pd.DataFrame(overlaps)

def generate_summary(overlaps_df, constrained_df):
    """Генерирует сводную статистику."""
    summary = {}
    
    # Общая статистика
    summary['total_constrained_regions'] = len(constrained_df)
    
    # Подсчитываем уникальные constrained регионы, которые имеют хотя бы одно перекрытие
    unique_regions = overlaps_df[['Chromosome', 'Constrained_Start', 'Constrained_End']].drop_duplicates()
    summary['regions_with_overlap'] = len(unique_regions)
    
    summary['total_overlaps'] = len(overlaps_df)
    
    # Создаем словари для подсчета уникальных регионов по типам features
    unique_regions_by_feature = {}
    contained_regions_by_feature = {}
    
    # Группируем данные по типам features и подсчитываем уникальные регионы
    feature_types = overlaps_df['Feature_Type'].unique()
    
    for feature_type in feature_types:
        # Фильтруем по типу feature
        feature_overlaps = overlaps_df[overlaps_df['Feature_Type'] == feature_type]
        
        # Подсчитываем уникальные constrained регионы для этого типа feature
        unique_regions_for_feature = feature_overlaps[['Chromosome', 'Constrained_Start', 'Constrained_End']].drop_duplicates()
        unique_regions_by_feature[feature_type] = len(unique_regions_for_feature)
        
        # Подсчитываем регионы, которые полностью содержатся внутри features
        # (Overlap_Length == Constrained_Length означает полное включение)
        contained_overlaps = feature_overlaps[feature_overlaps['Overlap_Length'] == feature_overlaps['Constrained_Length']]
        contained_regions = contained_overlaps[['Chromosome', 'Constrained_Start', 'Constrained_End']].drop_duplicates()
        contained_regions_by_feature[feature_type] = len(contained_regions)
    
    summary['unique_regions_by_feature'] = unique_regions_by_feature
    summary['contained_regions_by_feature'] = contained_regions_by_feature
    
    # Средний процент перекрытия
    if len(overlaps_df) > 0:
        summary['mean_overlap_percent'] = overlaps_df['Overlap_Percent'].mean()
    else:
        summary['mean_overlap_percent'] = 0
    
    return summary

def plot_results(overlaps_df, output_prefix):
    """Создает визуализации результатов."""
    # График распределения по типам features
    plt.figure(figsize=(10, 6))
    feature_counts = overlaps_df['Feature_Type'].value_counts()
    sns.barplot(x=feature_counts.index, y=feature_counts.values)
    plt.xticks(rotation=45, ha='right')
    plt.title('Distribution of Overlaps by Feature Type')
    plt.xlabel('Feature Type')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_feature_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # График процента перекрытия
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Feature_Type', y='Overlap_Percent', data=overlaps_df)
    plt.xticks(rotation=45, ha='right')
    plt.title('Overlap Percentage by Feature Type')
    plt.xlabel('Feature Type')
    plt.ylabel('Overlap Percentage')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_overlap_percentage.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Find overlaps between constrained regions and GFF features')
    parser.add_argument('input_table', help='Input table with constrained regions')
    parser.add_argument('gff_file', help='GFF file with genomic features')
    parser.add_argument('--output-prefix', default='overlap_results', help='Prefix for output files')
    
    args = parser.parse_args()
    
    # Чтение входной таблицы
    print("Reading input files...")
    constrained_df = pd.read_csv(args.input_table, sep='\t')
    features_df = parse_gff(args.gff_file)
    
    # Поиск пересечений
    print("Finding overlaps...")
    overlaps_df = find_overlaps(constrained_df, features_df)
    
    # Сохранение результатов
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"{args.output_prefix}_{timestamp}.tsv"
    overlaps_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to: {output_file}")
    
    # Генерация сводной статистики
    summary = generate_summary(overlaps_df, constrained_df)
    summary_file = f"{args.output_prefix}_summary_{timestamp}.txt"
    
    with open(summary_file, 'w') as f:
        f.write("=== Overlap Analysis Summary ===\n\n")
        f.write(f"Total constrained regions: {summary['total_constrained_regions']}\n")
        f.write(f"Unique constrained regions with overlaps: {summary['regions_with_overlap']}\n")
        f.write(f"Total overlap relationships found: {summary['total_overlaps']}\n")
        f.write(f"Mean overlap percentage: {summary['mean_overlap_percent']:.2f}%\n\n")
        
        f.write("Feature type distribution (unique constrained regions):\n")
        for feature_type, count in summary['unique_regions_by_feature'].items():
            percentage = (count / summary['regions_with_overlap']) * 100 if summary['regions_with_overlap'] > 0 else 0
            contained_count = summary['contained_regions_by_feature'][feature_type]
            contained_percentage = (contained_count / count) * 100 if count > 0 else 0
            
            f.write(f"  {feature_type}:\n")
            f.write(f"    - At least partially overlapping: {count} ({percentage:.1f}% of regions with overlap)\n")
            f.write(f"    - Completely contained within feature: {contained_count} ({contained_percentage:.1f}% of {feature_type} overlaps)\n")
        
    print(f"Summary saved to: {summary_file}")
    
    # Создание визуализаций
    if len(overlaps_df) > 0:
        print("Creating visualizations...")
        plot_results(overlaps_df, args.output_prefix)
        print(f"Plots saved as: {args.output_prefix}_*.png")
    else:
        print("No overlaps found, skipping visualizations")
    
    # Вывод первых нескольких результатов на экран
    if len(overlaps_df) > 0:
        print("\nFirst 10 overlaps:")
        print(overlaps_df[['Chromosome', 'Constrained_Start', 'Constrained_End', 
                          'Feature_Type', 'Feature_ID', 'Overlap_Start', 'Overlap_End', 
                          'Overlap_Length', 'Overlap_Percent']].head(10))

if __name__ == "__main__":
    main()