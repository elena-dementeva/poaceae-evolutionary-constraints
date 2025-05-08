import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from collections import defaultdict
from tqdm import tqdm

def parse_gff3(gff_file, feature_type=None):
    """Parse GFF3 file and extract features"""
    features = []
    feature_set = set()  # Для проверки дубликатов
    id_counter = {}  # Для генерации уникальных ID
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom = parts[0]
            feature = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            # Парсим атрибуты
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            if feature_type is None or feature == feature_type:
                # Проверяем, не является ли это дубликатом
                feature_key = (chrom, start, end, feature)
                if feature_key in feature_set:
                    continue  # Пропускаем дубликат
                
                feature_set.add(feature_key)
                
                # Генерируем уникальный ID
                if 'ID' in attr_dict and attr_dict['ID']:
                    feature_id = attr_dict['ID']
                else:
                    base_id = f"{chrom}:{start}-{end}"
                    if base_id not in id_counter:
                        id_counter[base_id] = 0
                        feature_id = base_id
                    else:
                        id_counter[base_id] += 1
                        feature_id = f"{base_id}_{id_counter[base_id]}"
                
                features.append({
                    'chrom': chrom,
                    'feature': feature,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'attributes': attr_dict,
                    'ID': feature_id
                })
    
    return features

def load_phyloP_data(phylop_file, chromosome=None):
    """Load phyloP data from file"""
    data = pd.read_csv(phylop_file, sep='\t')
    
    if chromosome:
        data = data[data['chrom'] == chromosome]
    
    # Создаем быстрый поиск по позициям
    phylop_dict = {}
    for chrom in data['chrom'].unique():
        chrom_data = data[data['chrom'] == chrom]
        phylop_dict[chrom] = dict(zip(chrom_data['position'], 
                                      zip(chrom_data['phyloP_score'], 
                                          chrom_data['status'].str.strip(), 
                                          chrom_data['significant'].str.strip())))
    
    return phylop_dict

def analyze_continuous_regions(feature, phylop_dict, min_length, region_type='all'):
    """Find continuous significant regions within a feature"""
    chrom = feature['chrom']
    start = feature['start']
    end = feature['end']
    
    if chrom not in phylop_dict:
        return []
    
    regions = []
    current_region = None
    
    for pos in range(start, end + 1):
        if pos in phylop_dict[chrom]:
            score, status, significant = phylop_dict[chrom][pos]
            
            if significant == 'yes' and (region_type == 'all' or status == region_type):
                if current_region is None:
                    current_region = {
                        'start': pos,
                        'end': pos,
                        'status': status,
                        'length': 1
                    }
                else:
                    # Проверяем, что позиция идет непосредственно за предыдущей
                    if pos == current_region['end'] + 1 and status == current_region['status']:
                        current_region['end'] = pos
                        current_region['length'] += 1
                    else:
                        if current_region['length'] >= min_length:
                            regions.append(current_region)
                        current_region = {
                            'start': pos,
                            'end': pos,
                            'status': status,
                            'length': 1
                        }
            else:
                if current_region is not None:
                    if current_region['length'] >= min_length:
                        regions.append(current_region)
                    current_region = None
        else:
            # Если позиция отсутствует, завершаем текущий регион
            if current_region is not None:
                if current_region['length'] >= min_length:
                    regions.append(current_region)
                current_region = None
    
    if current_region is not None and current_region['length'] >= min_length:
        regions.append(current_region)
    
    return regions

def analyze_window_regions(feature, phylop_dict, window_size, min_mean, overlap, region_type='all'):
    """Analyze feature using sliding window"""
    chrom = feature['chrom']
    start = feature['start']
    end = feature['end']
    
    if chrom not in phylop_dict:
        return []
    
    windows = []
    step = window_size - overlap
    
    for window_start in range(start, end - window_size + 2, step):
        window_end = min(window_start + window_size - 1, end)
        scores = []
        statuses = []
        
        for pos in range(window_start, window_end + 1):
            if pos in phylop_dict[chrom]:
                score, status, significant = phylop_dict[chrom][pos]
                if significant == 'yes' and (region_type == 'all' or status == region_type):
                    scores.append(abs(score))
                    statuses.append(status)
        
        if scores:
            mean_score = np.mean(scores)
            if mean_score >= min_mean:
                # Определяем преобладающий статус
                status_counts = defaultdict(int)
                for s in statuses:
                    status_counts[s] += 1
                predominant_status = max(status_counts, key=status_counts.get)
                
                windows.append({
                    'start': window_start,
                    'end': window_end,
                    'mean_score': mean_score,
                    'status': predominant_status
                })
    
    return windows

def analyze_regions_with_gaps(feature, phylop_dict, min_length, max_gaps, region_type='all'):
    """Find regions allowing gaps"""
    chrom = feature['chrom']
    start = feature['start']
    end = feature['end']
    
    if chrom not in phylop_dict:
        return []
    
    regions = []
    current_region = None
    gap_count = 0
    
    for pos in range(start, end + 1):
        if pos in phylop_dict[chrom]:
            score, status, significant = phylop_dict[chrom][pos]
            
            if significant == 'yes' and (region_type == 'all' or status == region_type):
                if current_region is None:
                    current_region = {
                        'start': pos,
                        'end': pos,
                        'status': status,
                        'significant_count': 1,
                        'total_length': 1
                    }
                    gap_count = 0
                else:
                    if status == current_region['status']:
                        current_gap = pos - current_region['end'] - 1
                        gap_count += current_gap
                        
                        if gap_count <= max_gaps:
                            current_region['end'] = pos
                            current_region['significant_count'] += 1
                            current_region['total_length'] = pos - current_region['start'] + 1
                        else:
                            if current_region['total_length'] >= min_length:
                                regions.append(current_region)
                            current_region = {
                                'start': pos,
                                'end': pos,
                                'status': status,
                                'significant_count': 1,
                                'total_length': 1
                            }
                            gap_count = 0
                    else:
                        if current_region['total_length'] >= min_length:
                            regions.append(current_region)
                        current_region = {
                            'start': pos,
                            'end': pos,
                            'status': status,
                            'significant_count': 1,
                            'total_length': 1
                        }
                        gap_count = 0
        else:
            if current_region is not None:
                if gap_count < max_gaps:
                    gap_count += 1
                    current_region['end'] = pos
                    current_region['total_length'] = pos - current_region['start'] + 1
                else:
                    if current_region['total_length'] >= min_length:
                        regions.append(current_region)
                    current_region = None
                    gap_count = 0
    
    if current_region is not None and current_region['total_length'] >= min_length:
        regions.append(current_region)
    
    return regions

def calculate_feature_mean(feature, phylop_dict, region_type='all'):
    """Calculate mean score for entire feature"""
    chrom = feature['chrom']
    start = feature['start']
    end = feature['end']
    
    if chrom not in phylop_dict:
        return 0, 0
    
    scores = []
    for pos in range(start, end + 1):
        if pos in phylop_dict[chrom]:
            score, status, significant = phylop_dict[chrom][pos]
            if significant == 'yes' and (region_type == 'all' or status == region_type):
                scores.append(abs(score))
    
    if scores:
        return np.mean(scores), len(scores)
    else:
        return 0, 0

def plot_conservation_pattern(feature, phylop_dict, output_file):
    """Plot conservation pattern for a feature"""
    chrom = feature['chrom']
    start = feature['start']
    end = feature['end']
    
    if chrom not in phylop_dict:
        return
    
    positions = []
    scores = []
    colors = []
    
    for pos in range(start, end + 1):
        if pos in phylop_dict[chrom]:
            score, status, significant = phylop_dict[chrom][pos]
            if significant == 'yes':
                positions.append(pos)
                scores.append(score)
                colors.append('red' if status == 'accelerated' else 'blue')
    
    if not positions:
        return
    
    plt.figure(figsize=(12, 6))
    plt.scatter(positions, scores, c=colors, alpha=0.7)
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.title(f'Conservation pattern for {feature["ID"]}')
    plt.xlabel('Genomic position')
    plt.ylabel('PhyloP score')
    
    # Добавляем легенду
    blue_patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Constrained')
    red_patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Accelerated')
    plt.legend(handles=[blue_patch, red_patch])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze conservation patterns in genomic features')
    parser.add_argument('--phylop', required=True, help='PhyloP TSV file')
    parser.add_argument('--gff', required=True, help='GFF3 file with features')
    parser.add_argument('--feature-type', default='CDS', help='Feature type to analyze (default: CDS)')
    parser.add_argument('--chromosome', help='Analyze only this chromosome')
    parser.add_argument('--region-type', choices=['all', 'constrained', 'accelerated'], default='all',
                        help='Type of regions to analyze (default: all)')
    
    # Параметры для анализа
    parser.add_argument('--min-length', type=int, default=20, 
                        help='Minimum length for continuous regions')
    parser.add_argument('--window-size', type=int, default=100, 
                        help='Window size for sliding window analysis')
    parser.add_argument('--min-mean', type=float, default=1.5,
                        help='Minimum mean score for window analysis')
    parser.add_argument('--overlap', type=int, default=50, 
                        help='Window overlap')
    parser.add_argument('--max-gaps', type=int, default=3, 
                        help='Maximum allowed gaps')
    
    # Дополнительные параметры
    parser.add_argument('--specific-features', nargs='+', help='Specific feature positions (format: chrom:start)')
    parser.add_argument('--plot-patterns', action='store_true', 
                        help='Plot conservation patterns for features')
    parser.add_argument('--output-prefix', default='feature_analysis', 
                        help='Prefix for output files')
    
    args = parser.parse_args()
    
    try:
        # Загружаем данные phyloP
        print("Loading phyloP data...")
        phylop_dict = load_phyloP_data(args.phylop, args.chromosome)
        
        # Загружаем features из GFF
        print(f"Loading features from GFF ({args.feature_type})...")
        features = parse_gff3(args.gff, args.feature_type)
        
        # Фильтруем по хромосоме, если указано
        if args.chromosome:
            features = [f for f in features if f['chrom'] == args.chromosome]
        
        # Добавляем конкретные features, если указаны
        if args.specific_features:
            existing_coords = {(f['chrom'], f['start']) for f in features}
            
            for spec in args.specific_features:
                chrom, pos = spec.split(':')
                start = int(pos)
                
                # Проверяем, не существует ли уже такой feature
                if (chrom, start) not in existing_coords:
                    features.append({
                        'chrom': chrom,
                        'feature': 'specific',
                        'start': start,
                        'end': start + 100,  # Предполагаем длину 100
                        'strand': '+',
                        'attributes': {},
                        'ID': f"{chrom}:{start}_specific"
                    })
                    existing_coords.add((chrom, start))
        
        print(f"Analyzing {len(features)} features...")
        
        # Счетчики для статистики
        continuous_stats = {
            'total_features_with_continuous': 0,
            'constrained_regions': 0,
            'accelerated_regions': 0,
            'constrained_total_length': 0,
            'accelerated_total_length': 0,
            'features_with_constrained': set(),  # Добавлены новые счетчики
            'features_with_accelerated': set()   # Добавлены новые счетчики
        }
        
        window_stats = {
            'count': 0,
            'constrained_count': 0,
            'accelerated_count': 0,
            'mean_scores': [],
            'features_with_constrained': set(),  # Добавлены новые счетчики
            'features_with_accelerated': set()   # Добавлены новые счетчики
        }
        
        gaps_stats = {
            'count': 0,
            'constrained_count': 0,
            'accelerated_count': 0,
            'total_length': 0,
            'features_with_constrained': set(),  # Добавлены новые счетчики
            'features_with_accelerated': set()   # Добавлены новые счетчики
        }
        
        feature_means = []
        
        # Сохраняем подробные результаты
        results_file = f"{args.output_prefix}_detailed_results.tsv"
        with open(results_file, 'w') as f:
            f.write("Feature_ID\tChromosome\tStart\tEnd\t"
                    "Continuous_regions_constrained\tContinuous_regions_accelerated\t"
                    "Continuous_constrained_coords\tContinuous_accelerated_coords\t"
                    "Window_regions\tWindow_constrained\tWindow_accelerated\t"
                    "Gaps_regions\tGaps_constrained\tGaps_accelerated\t"
                    "Mean_score\n")
            
            # Анализируем каждый feature
            for i, feature in enumerate(tqdm(features, desc="Processing features")):
                feature_id = feature['ID']  # Для отслеживания уникальных фичей
                
                # 1. Непрерывные регионы
                continuous_regions = analyze_continuous_regions(feature, phylop_dict, args.min_length, args.region_type)
                
                # Разделяем continuous regions по типу
                constrained_regions = [r for r in continuous_regions if r['status'] == 'constrained']
                accelerated_regions = [r for r in continuous_regions if r['status'] == 'accelerated']
                
                if continuous_regions:
                    continuous_stats['total_features_with_continuous'] += 1
                
                if constrained_regions:
                    continuous_stats['constrained_regions'] += len(constrained_regions)
                    continuous_stats['features_with_constrained'].add(feature_id)
                    for region in constrained_regions:
                        continuous_stats['constrained_total_length'] += region['length']

                if accelerated_regions:
                    continuous_stats['accelerated_regions'] += len(accelerated_regions)
                    continuous_stats['features_with_accelerated'].add(feature_id)
                    for region in accelerated_regions:
                        continuous_stats['accelerated_total_length'] += region['length']
                
                # Формируем координаты регионов
                constrained_coords = [f"{r['start']}:{r['end']}" for r in constrained_regions]
                accelerated_coords = [f"{r['start']}:{r['end']}" for r in accelerated_regions]
                
                # 2. Window analysis
                window_regions = analyze_window_regions(feature, phylop_dict, 
                                                      args.window_size, args.min_mean, args.overlap, args.region_type)
                
                window_constrained = [w for w in window_regions if w['status'] == 'constrained']
                window_accelerated = [w for w in window_regions if w['status'] == 'accelerated']
                
                if window_regions:
                    window_stats['count'] += 1
                    
                    if window_constrained:
                        window_stats['features_with_constrained'].add(feature_id)
                    
                    if window_accelerated:
                        window_stats['features_with_accelerated'].add(feature_id)
                    
                    for window in window_regions:
                        window_stats['mean_scores'].append(window['mean_score'])
                        if window['status'] == 'constrained':
                            window_stats['constrained_count'] += 1
                        else:
                            window_stats['accelerated_count'] += 1
                
                # 3. Регионы с пропусками
                gaps_regions = analyze_regions_with_gaps(feature, phylop_dict, 
                                                        args.min_length, args.max_gaps, args.region_type)
                
                gaps_constrained = [r for r in gaps_regions if r['status'] == 'constrained']
                gaps_accelerated = [r for r in gaps_regions if r['status'] == 'accelerated']
                
                if gaps_regions:
                    gaps_stats['count'] += 1
                    
                    if gaps_constrained:
                        gaps_stats['features_with_constrained'].add(feature_id)
                    
                    if gaps_accelerated:
                        gaps_stats['features_with_accelerated'].add(feature_id)
                    
                    for region in gaps_regions:
                        gaps_stats['total_length'] += region['total_length']
                        if region['status'] == 'constrained':
                            gaps_stats['constrained_count'] += 1
                        else:
                            gaps_stats['accelerated_count'] += 1
                
                # 4. Среднее по feature
                mean_score, sig_count = calculate_feature_mean(feature, phylop_dict, args.region_type)
                if sig_count > 0:
                    feature_means.append(mean_score)
                
                # Записываем результаты
                f.write(f"{feature['ID']}\t{feature['chrom']}\t{feature['start']}\t{feature['end']}\t"
                        f"{len(constrained_regions)}\t{len(accelerated_regions)}\t"
                        f"{','.join(constrained_coords) if constrained_coords else '-'}\t"
                        f"{','.join(accelerated_coords) if accelerated_coords else '-'}\t"
                        f"{len(window_regions)}\t{len(window_constrained)}\t{len(window_accelerated)}\t"
                        f"{len(gaps_regions)}\t{len(gaps_constrained)}\t{len(gaps_accelerated)}\t"
                        f"{mean_score:.3f}\n")
                
                # 5. Построение графика, если указано
                if args.plot_patterns and i < 10:  # Ограничиваем количество графиков
                    output_file = f"{args.output_prefix}_pattern_{i}.png"
                    plot_conservation_pattern(feature, phylop_dict, output_file)
        
        # Выводим статистику
        print("\n=== Analysis Results ===")
        print(f"Total features analyzed: {len(features)}")
        print(f"Region type filter: {args.region_type}")
        
        print("\n1. Continuous regions:")
        print(f"   Features with continuous regions: {continuous_stats['total_features_with_continuous']} "
              f"({continuous_stats['total_features_with_continuous']/len(features)*100:.1f}%)")
        
        total_regions = continuous_stats['constrained_regions'] + continuous_stats['accelerated_regions']
        if total_regions > 0:
            print(f"   Total continuous regions: {total_regions}")
            print(f"   Constrained regions: {continuous_stats['constrained_regions']}")
            print(f"   Accelerated regions: {continuous_stats['accelerated_regions']}")
            print(f"   Features with constrained regions: {len(continuous_stats['features_with_constrained'])} "
                 f"({len(continuous_stats['features_with_constrained'])/len(features)*100:.1f}%)")
            print(f"   Features with accelerated regions: {len(continuous_stats['features_with_accelerated'])} "
                 f"({len(continuous_stats['features_with_accelerated'])/len(features)*100:.1f}%)")
            
            if continuous_stats['constrained_regions'] > 0:
                avg_constrained_len = continuous_stats['constrained_total_length'] / continuous_stats['constrained_regions']
                print(f"   Average constrained region length: {avg_constrained_len:.1f}")
            
            if continuous_stats['accelerated_regions'] > 0:
                avg_accelerated_len = continuous_stats['accelerated_total_length'] / continuous_stats['accelerated_regions']
                print(f"   Average accelerated region length: {avg_accelerated_len:.1f}")
        
        print("\n2. Window analysis:")
        print(f"   Features with significant windows: {window_stats['count']} "
              f"({window_stats['count']/len(features)*100:.1f}%)")
        if window_stats['count'] > 0:
            print(f"   Constrained windows: {window_stats['constrained_count']}")
            print(f"   Accelerated windows: {window_stats['accelerated_count']}")
            print(f"   Features with constrained windows: {len(window_stats['features_with_constrained'])} "
                 f"({len(window_stats['features_with_constrained'])/len(features)*100:.1f}%)")
            print(f"   Features with accelerated windows: {len(window_stats['features_with_accelerated'])} "
                 f"({len(window_stats['features_with_accelerated'])/len(features)*100:.1f}%)")
            print(f"   Average window score: {np.mean(window_stats['mean_scores']):.3f}")
        
        print("\n3. Regions with gaps:")
        print(f"   Features with gapped regions: {gaps_stats['count']} "
              f"({gaps_stats['count']/len(features)*100:.1f}%)")
        if gaps_stats['count'] > 0:
            print(f"   Constrained regions: {gaps_stats['constrained_count']}")
            print(f"   Accelerated regions: {gaps_stats['accelerated_count']}")
            print(f"   Features with constrained gapped regions: {len(gaps_stats['features_with_constrained'])} "
                 f"({len(gaps_stats['features_with_constrained'])/len(features)*100:.1f}%)")
            print(f"   Features with accelerated gapped regions: {len(gaps_stats['features_with_accelerated'])} "
                 f"({len(gaps_stats['features_with_accelerated'])/len(features)*100:.1f}%)")
            print(f"   Average region length: {gaps_stats['total_length']/gaps_stats['count']:.1f}")
        
        print("\n4. Feature-wide analysis:")
        if feature_means:
            print(f"   Average feature score: {np.mean(feature_means):.3f}")
            print(f"   Median feature score: {np.median(feature_means):.3f}")
            print(f"   Min feature score: {min(feature_means):.3f}")
            print(f"   Max feature score: {max(feature_means):.3f}")
        
        print(f"\nDetailed results saved to: {results_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()