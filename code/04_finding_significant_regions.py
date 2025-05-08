import pandas as pd
import numpy as np
import sys
import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import argparse

def count_lines(filename):
    """Fast line count in file"""
    with open(filename, 'rb') as f:
        lines = 0
        buf_size = 1024 * 1024
        read_f = f.raw.read
        
        buf = read_f(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = read_f(buf_size)
    
    return lines

def get_chunk_boundaries(file_path, n_chunks):
    """Divides file into approximately equal parts for processing"""
    file_size = os.path.getsize(file_path)
    chunk_size = file_size // n_chunks
    boundaries = []
    
    with open(file_path, 'r') as f:
        # Считываем заголовок
        header = f.readline()
        boundaries.append(len(header))
        
        # Находим границы чанков
        for i in range(1, n_chunks):
            target_pos = i * chunk_size
            f.seek(target_pos)
            f.readline()  # Пропускаем неполную строку
            boundaries.append(f.tell())
    
    # Добавляем конец файла
    boundaries.append(file_size)
    return boundaries

def find_regions_in_chunk(args):
    """Finds significant regions in a file chunk"""
    file_path, start_pos, end_pos, mode, params, status_type = args
    
    # Для хранения результатов
    regions = []
    lines_data = []
    
    with open(file_path, 'r') as f:
        # Перемещаемся к началу чанка
        f.seek(start_pos)
        
        # Если это первый чанк, то пропускаем заголовок
        if start_pos == 0:
            header = f.readline().strip().split('\t')
        else:
            # Для не первого чанка определяем столбцы
            with open(file_path, 'r') as header_file:
                header = header_file.readline().strip().split('\t')
        
        columns = {name: i for i, name in enumerate(header)}
        
        # Индексы нужных столбцов
        chrom_idx = columns.get('chrom', 0)
        pos_idx = columns.get('position', 1)
        sig_idx = columns.get('significant', 6)  # Исправил индекс для вашей таблицы
        status_idx = columns.get('status', 3)
        score_idx = columns.get('phyloP_score', 2)
        
        while f.tell() < end_pos:
            line = f.readline()
            if not line:
                break
                
            parts = line.strip().split('\t')
            if len(parts) <= max(chrom_idx, pos_idx, sig_idx, status_idx, score_idx):
                continue
                
            try:
                chrom = parts[chrom_idx]
                pos = int(parts[pos_idx])
                is_significant = parts[sig_idx] == 'yes'
                status = parts[status_idx]
                score = float(parts[score_idx])
                
                lines_data.append({
                    'chrom': chrom,
                    'position': pos,
                    'significant': is_significant,
                    'status': status,
                    'score': score
                })
            except (ValueError, IndexError):
                continue
    
    # Обработка в зависимости от режима
    if mode == 'continuous':
        regions = find_continuous_regions(lines_data, params['min_length'], status_type)
    elif mode == 'window':
        regions = find_window_regions(lines_data, params['window_size'], params['min_mean'], 
                                      params['overlap'], status_type)
    elif mode == 'gaps':
        regions = find_regions_with_gaps(lines_data, params['min_length'], 
                                         params['max_gaps'], status_type)
    
    return regions

def find_continuous_regions(data, min_length, status_type):
    """Find continuous significant regions"""
    regions = []
    current_region = None
    
    for i, line in enumerate(data):
        if line['significant'] and (status_type == 'all' or line['status'] == status_type):
            if current_region is None:
                current_region = {
                    'chrom': line['chrom'],
                    'start': line['position'],
                    'end': line['position'],
                    'status': line['status']
                }
            else:
                if (line['chrom'] == current_region['chrom'] and 
                    line['position'] == current_region['end'] + 1 and 
                    line['status'] == current_region['status']):
                    current_region['end'] = line['position']
                else:
                    if (current_region['end'] - current_region['start'] + 1) >= min_length:
                        regions.append(current_region)
                    current_region = {
                        'chrom': line['chrom'],
                        'start': line['position'],
                        'end': line['position'],
                        'status': line['status']
                    }
        else:
            if current_region is not None:
                if (current_region['end'] - current_region['start'] + 1) >= min_length:
                    regions.append(current_region)
                current_region = None
    
    if current_region is not None and (current_region['end'] - current_region['start'] + 1) >= min_length:
        regions.append(current_region)
    
    return regions

def find_window_regions(data, window_size, min_mean, overlap, status_type):
    """Find regions using sliding window"""
    regions = []
    step = window_size - overlap
    
    # Group by chromosome
    from itertools import groupby
    data_sorted = sorted(data, key=lambda x: (x['chrom'], x['position']))
    
    for chrom, chrom_data in groupby(data_sorted, key=lambda x: x['chrom']):
        chrom_data = list(chrom_data)
        positions = [d['position'] for d in chrom_data]
        scores = [d['score'] if d['significant'] and (status_type == 'all' or d['status'] == status_type) else 0 
                  for d in chrom_data]
        
        for i in range(0, len(positions) - window_size + 1, step):
            window_scores = scores[i:i + window_size]
            mean_score = np.mean([abs(s) for s in window_scores if s != 0])
            
            if mean_score >= min_mean:
                # Determine predominant status
                significant_entries = [chrom_data[j] for j in range(i, i + window_size) 
                                       if scores[j] != 0]
                if significant_entries:
                    status_counts = {}
                    for entry in significant_entries:
                        status_counts[entry['status']] = status_counts.get(entry['status'], 0) + 1
                    predominant_status = max(status_counts, key=status_counts.get)
                    
                    regions.append({
                        'chrom': chrom,
                        'start': positions[i],
                        'end': positions[i + window_size - 1],
                        'status': predominant_status,
                        'mean_score': mean_score
                    })
    
    return regions

def find_regions_with_gaps(data, min_length, max_gaps, status_type):
    """Find regions allowing gaps"""
    regions = []
    current_region = None
    gap_count = 0
    
    for i, line in enumerate(data):
        if line['significant'] and (status_type == 'all' or line['status'] == status_type):
            if current_region is None:
                current_region = {
                    'chrom': line['chrom'],
                    'start': line['position'],
                    'end': line['position'],
                    'status': line['status'],
                    'significant_count': 1
                }
                gap_count = 0
            else:
                if (line['chrom'] == current_region['chrom'] and 
                    line['position'] <= current_region['end'] + max_gaps + 1 and 
                    line['status'] == current_region['status']):
                    
                    # Подсчитываем пропуски
                    current_gap = line['position'] - current_region['end'] - 1
                    gap_count += current_gap
                    
                    if gap_count <= max_gaps:
                        current_region['end'] = line['position']
                        current_region['significant_count'] += 1
                    else:
                        # Завершаем текущий регион и начинаем новый
                        if (current_region['end'] - current_region['start'] + 1) >= min_length:
                            regions.append(current_region)
                        current_region = {
                            'chrom': line['chrom'],
                            'start': line['position'],
                            'end': line['position'],
                            'status': line['status'],
                            'significant_count': 1
                        }
                        gap_count = 0
                else:
                    if (current_region['end'] - current_region['start'] + 1) >= min_length:
                        regions.append(current_region)
                    current_region = {
                        'chrom': line['chrom'],
                        'start': line['position'],
                        'end': line['position'],
                        'status': line['status'],
                        'significant_count': 1
                    }
                    gap_count = 0
        else:
            if current_region is not None:
                # Проверяем, можем ли продолжить с пропуском
                if (line['chrom'] == current_region['chrom'] and 
                    line['position'] == current_region['end'] + 1):
                    if gap_count < max_gaps:
                        current_region['end'] = line['position']
                        gap_count += 1
                    else:
                        if (current_region['end'] - current_region['start'] + 1) >= min_length:
                            regions.append(current_region)
                        current_region = None
                        gap_count = 0
                else:
                    if (current_region['end'] - current_region['start'] + 1) >= min_length:
                        regions.append(current_region)
                    current_region = None
                    gap_count = 0
    
    if current_region is not None and (current_region['end'] - current_region['start'] + 1) >= min_length:
        regions.append(current_region)
    
    return regions

def merge_boundary_regions(chunk_results):
    """Merges regions at chunk boundaries"""
    all_regions = []
    
    for regions in chunk_results:
        all_regions.extend(regions)
    
    # Sort and merge overlapping regions
    if not all_regions:
        return []
    
    all_regions.sort(key=lambda x: (x['chrom'], x['start']))
    merged = [all_regions[0]]
    
    for current in all_regions[1:]:
        previous = merged[-1]
        if (current['chrom'] == previous['chrom'] and 
            current['start'] <= previous['end'] + 1 and 
            current['status'] == previous['status']):
            previous['end'] = max(previous['end'], current['end'])
            if 'significant_count' in current:
                previous['significant_count'] = previous.get('significant_count', 0) + current['significant_count']
            if 'mean_score' in current:
                previous['mean_score'] = (previous.get('mean_score', 0) + current['mean_score']) / 2
        else:
            merged.append(current)
    
    return merged

def process_file_parallel(file_path, mode, params, status_type='all', n_processes=None):
    """Parallel file processing to find significant regions"""
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)
    
    print(f"Using {n_processes} processes")
    
    # Разделяем файл на чанки
    boundaries = get_chunk_boundaries(file_path, n_processes)
    chunk_args = [(file_path, boundaries[i], boundaries[i+1], mode, params, status_type) 
                  for i in range(len(boundaries)-1)]
    
    # Обрабатываем чанки параллельно
    results = []
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [executor.submit(find_regions_in_chunk, arg) for arg in chunk_args]
        
        for future in tqdm(futures, desc="Processing chunks"):
            result = future.result()
            results.append(result)
    
    # Объединяем результаты
    print("Merging chunk results...")
    merged_regions = merge_boundary_regions(results)
    
    return merged_regions

def main():
    parser = argparse.ArgumentParser(description='Finding significant regions in phyloP file')
    parser.add_argument('input_file', help='Path to TSV file with results')
    parser.add_argument('--mode', choices=['continuous', 'window', 'gaps'], default='continuous',
                        help='Mode of region finding')
    parser.add_argument('--status', choices=['constrained', 'accelerated', 'all'], default='all',
                        help='Status type to consider')
    
    # Параметры для режима continuous
    parser.add_argument('--min-length', type=int, default=20, 
                        help='Minimum region length for continuous mode')
    
    # Параметры для режима window
    parser.add_argument('--window-size', type=int, default=100, 
                        help='Window size for window mode')
    parser.add_argument('--min-mean', type=float, default=1.5, 
                        help='Minimum mean score for window mode')
    parser.add_argument('--overlap', type=int, default=50, 
                        help='Window overlap for window mode')
    
    # Параметры для режима gaps
    parser.add_argument('--min-length-gaps', type=int, default=20, 
                        help='Minimum region length for gaps mode')
    parser.add_argument('--max-gaps', type=int, default=3, 
                        help='Maximum allowed gaps')
    
    parser.add_argument('--processes', type=int, default=None, 
                        help='Number of processes (default: CPU count - 1)')
    parser.add_argument('--output', help='Path to output file (default: auto)')
    
    args = parser.parse_args()
    
    # Определяем имя выходного файла
    output_file = args.output
    if not output_file:
        output_file = args.input_file.replace('.tsv', '').replace('.txt', '') + f'_regions_{args.mode}.tsv'
    
    try:
        # Проверяем наличие файла
        if not os.path.exists(args.input_file):
            print(f"Error: File {args.input_file} not found")
            sys.exit(1)
        
        # Подготавливаем параметры в зависимости от режима
        if args.mode == 'continuous':
            params = {'min_length': args.min_length}
        elif args.mode == 'window':
            params = {
                'window_size': args.window_size,
                'min_mean': args.min_mean,
                'overlap': args.overlap
            }
        elif args.mode == 'gaps':
            params = {
                'min_length': args.min_length_gaps,
                'max_gaps': args.max_gaps
            }
        
        # Ищем значимые участки
        print(f"Finding significant regions in {args.mode} mode...")
        regions = process_file_parallel(args.input_file, args.mode, params, args.status, args.processes)
        
        # Сохраняем результаты
        print(f"Saving results to {output_file}...")
        with open(output_file, 'w') as f:
            if args.mode == 'window':
                f.write("Chromosome\tStart\tEnd\tStatus\tMean_Score\tLength\n")
                for region in regions:
                    length = region['end'] - region['start'] + 1
                    f.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t{region['status']}\t{region['mean_score']:.3f}\t{length}\n")
            elif args.mode == 'gaps':
                f.write("Chromosome\tStart\tEnd\tStatus\tSignificant_Count\tLength\n")
                for region in regions:
                    length = region['end'] - region['start'] + 1
                    f.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t{region['status']}\t{region['significant_count']}\t{length}\n")
            else:
                f.write("Chromosome\tStart\tEnd\tStatus\tLength\n")
                for region in regions:
                    length = region['end'] - region['start'] + 1
                    f.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t{region['status']}\t{length}\n")
        
        # Выводим статистику
        constrained = [r for r in regions if r['status'] == 'constrained']
        accelerated = [r for r in regions if r['status'] == 'accelerated']
        
        print("\nStatistics:")
        print(f"Total significant regions: {len(regions)}")
        if regions:
            print(f"Constrained regions: {len(constrained)} ({len(constrained)/len(regions)*100:.1f}%)")
            print(f"Accelerated regions: {len(accelerated)} ({len(accelerated)/len(regions)*100:.1f}%)")
            
            # Распределение длин
            lengths = [r['end'] - r['start'] + 1 for r in regions]
            print(f"\nLength distribution:")
            print(f"Minimum length: {min(lengths)}")
            print(f"Maximum length: {max(lengths)}")
            print(f"Average length: {sum(lengths)/len(lengths):.1f}")
            print(f"Median length: {sorted(lengths)[len(lengths)//2]}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
