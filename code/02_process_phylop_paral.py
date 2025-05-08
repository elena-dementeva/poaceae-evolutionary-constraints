import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse
import re
from scipy.stats import norm
import matplotlib.pyplot as plt
import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from functools import partial
import time
import psutil
import logging

# Настройка базового логирования
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_header(line):
    """Parse a fixedStep header line and extract chromosome, position and step."""
    chrom_match = re.search(r'chrom=(\S+)', line)
    pos_match = re.search(r'start=(\d+)', line)
    step_match = re.search(r'step=(\d+)', line)
    
    if chrom_match and pos_match:
        current_chrom = chrom_match.group(1)
        current_pos = int(pos_match.group(1))
        current_step = int(step_match.group(1)) if step_match else 1
        return current_chrom, current_pos, current_step
    return None, None, None

def analyze_file_structure(input_file):
    """
    Анализирует структуру WIG файла, определяя все блоки fixedStep.
    Возвращает список блоков с их позициями и заголовками.
    """
    blocks = []
    total_lines = 0
    data_lines = 0
    header_lines = 0
    
    with open(input_file, 'r') as f:
        file_pos = 0
        current_block = None
        
        for line_num, line in enumerate(f):
            total_lines += 1
            
            if line.startswith('fixedStep'):
                header_lines += 1
                chrom, start, step = parse_header(line)
                
                # Сохраняем предыдущий блок, если он существует
                if current_block:
                    current_block['end_pos'] = file_pos
                    current_block['line_count'] = line_num - current_block['start_line']
                    blocks.append(current_block)
                
                # Создаем новый блок
                current_block = {
                    'header': line.strip(),
                    'chrom': chrom,
                    'start': start,
                    'step': step,
                    'file_pos': file_pos,
                    'start_line': line_num,
                    'values': []
                }
            elif current_block is not None:
                # Это строка с данными
                data_lines += 1
                
            file_pos += len(line)
        
        # Добавляем последний блок
        if current_block:
            current_block['end_pos'] = file_pos
            current_block['line_count'] = line_num - current_block['start_line'] if 'line_num' in locals() else 0
            blocks.append(current_block)
    
    logger.info(f"File analysis: {len(blocks)} blocks, {total_lines} total lines, {data_lines} data lines")
    return blocks, total_lines, data_lines

def read_block_data(input_file, block_info):
    """
    Читает данные из одного блока WIG файла.
    
    Args:
        input_file: Путь к файлу
        block_info: Словарь с информацией о блоке (из analyze_file_structure)
        
    Returns:
        Кортеж (chrom, positions, scores) или None в случае ошибки
    """
    chrom = block_info['chrom']
    start_pos = block_info['start']
    step = block_info['step']
    file_pos = block_info['file_pos']
    
    try:
        with open(input_file, 'r') as f:
            f.seek(file_pos)
            # Пропускаем заголовок
            header = f.readline()
            
            positions = []
            scores = []
            current_pos = start_pos
            
            # Читаем данные до следующего заголовка или конца файла
            while True:
                line = f.readline()
                if not line or line.startswith('fixedStep'):
                    break
                    
                line = line.strip()
                if not line:
                    continue
                    
                try:
                    score = float(line)
                    positions.append(current_pos)
                    scores.append(score)
                    current_pos += step
                except ValueError:
                    continue
            
            return chrom, positions, scores
    except Exception as e:
        logger.error(f"Error reading block at position {file_pos}: {e}")
        return None

def process_blocks(input_file, blocks_to_process, process_id=0):
    """
    Обрабатывает набор блоков из WIG файла.
    
    Args:
        input_file: Путь к файлу
        blocks_to_process: Список блоков для обработки (из analyze_file_structure)
        process_id: ID процесса для логирования
        
    Returns:
        Список словарей с результатами обработки
    """
    results = []
    
    for block in blocks_to_process:
        block_data = read_block_data(input_file, block)
        if not block_data:
            continue
            
        chrom, positions, scores = block_data
        
        # Преобразуем в массивы NumPy для ускорения
        positions_array = np.array(positions)
        scores_array = np.array(scores)
        
        # Создаем результаты
        block_results = []
        for j in range(len(scores_array)):
            status = "constrained" if scores_array[j] > 0 else "accelerated" if scores_array[j] < 0 else "neutral"
            
            result = {
                "chrom": chrom,
                "position": positions_array[j],
                "phyloP_score": scores_array[j],
                "status": status
            }
            block_results.append(result)
        
        results.extend(block_results)
    
    return results

def distribute_blocks(blocks, n_processes):
    """
    Распределяет блоки между процессами, чтобы обеспечить равномерную нагрузку.
    Возвращает список списков блоков для каждого процесса.
    """
    # Создаем пустые списки для каждого процесса
    process_assignments = [[] for _ in range(n_processes)]
    
    # Распределяем блоки между процессами поочередно
    for i, block in enumerate(blocks):
        process_idx = i % n_processes
        process_assignments[process_idx].append(block)
    
    return process_assignments

def apply_statistical_analysis(df, alpha=0.05, fdr_combined=False):
    """Apply statistical analysis on the complete dataset."""
    logger.info(f"Applying statistical analysis with alpha={alpha}, fdr_combined={fdr_combined}")
    
    # Get all scores as numpy array
    scores = df['phyloP_score'].values
    
    # Split data into positive and negative scores
    positive_mask = scores > 0
    negative_mask = scores < 0
    
    positive_count = np.sum(positive_mask)
    negative_count = np.sum(negative_mask)
    neutral_count = np.sum(scores == 0)
    
    logger.info(f"Score distribution: {positive_count} positive, {negative_count} negative, {neutral_count} neutral")
    
    # Initialize arrays for results
    p_values = np.full_like(scores, np.nan, dtype=float)
    q_values = np.full_like(scores, np.nan, dtype=float)
    significant = np.zeros_like(scores, dtype=bool)
    
    # Convert phyloP scores to p-values
    # For both positive and negative scores, p-value = 10^(-|phyloP|)
    p_values = np.power(10, -np.abs(scores))
    
    if fdr_combined:
        # Process all scores together regardless of sign
        logger.info("Using combined FDR correction for all scores")
        _, q_values_all, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
        q_values = q_values_all
        significant = q_values <= alpha
        
        logger.info(f"Combined FDR: {np.sum(significant)} significant positions out of {len(significant)}")
    else:
        # Process positive scores (constrained regions)
        if np.any(positive_mask):
            p_values_positive = p_values[positive_mask]
            _, q_values_positive, _, _ = multipletests(p_values_positive, alpha=alpha, method='fdr_bh')
            q_values[positive_mask] = q_values_positive
            significant[positive_mask] = q_values_positive <= alpha
        
        # Process negative scores (accelerated regions)
        if np.any(negative_mask):
            p_values_negative = p_values[negative_mask]
            _, q_values_negative, _, _ = multipletests(p_values_negative, alpha=alpha, method='fdr_bh')
            q_values[negative_mask] = q_values_negative
            significant[negative_mask] = q_values_negative <= alpha
    
    # Add results to DataFrame
    df['p_value'] = p_values
    df['q_value'] = q_values
    df['significant'] = significant
    df['significant'] = df['significant'].map({True: 'yes', False: 'no'})
    
    return df

def generate_plots(df, output_file):
    """Generate diagnostic plots from results DataFrame."""
    plot_dir = os.path.dirname(output_file)
    if not plot_dir:
        plot_dir = "."
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    
    plot_prefix = os.path.splitext(os.path.basename(output_file))[0]
    
    # Convert scores column to numeric if it's not already
    scores = pd.to_numeric(df['phyloP_score'])
    p_values = pd.to_numeric(df['p_value'].dropna())
    
    logger.info(f"Generating plots with {len(scores)} values")
    
    # Plot 1: Histogram of phyloP scores (using pandas for speed)
    plt.figure(figsize=(10, 6))
    scores.plot.hist(bins=50, alpha=0.7)
    plt.axvline(x=0, color='r', linestyle='--')
    plt.title("Distribution of phyloP Scores")
    plt.xlabel("phyloP Score")
    plt.ylabel("Frequency")
    plt.savefig(f"{plot_dir}/{plot_prefix}_score_histogram.png")
    plt.close()
    
    # Plot 2: Distribution of p-values
    plt.figure(figsize=(10, 6))
    p_values.plot.hist(bins=20, alpha=0.7)
    plt.title("Distribution of p-values")
    plt.xlabel("p-value")
    plt.ylabel("Frequency")
    plt.savefig(f"{plot_dir}/{plot_prefix}_pvalue_histogram.png")
    plt.close()
    
    # Plot 3: Q-Q plot of p-values
    plt.figure(figsize=(8, 8))
    sorted_pvals = np.sort(p_values)
    n = len(sorted_pvals)
    expected = np.linspace(0, 1, n+2)[1:-1]
    plt.plot(expected, sorted_pvals, 'o', markersize=2)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.title("Q-Q Plot of p-values")
    plt.xlabel("Expected p-value")
    plt.ylabel("Observed p-value")
    plt.savefig(f"{plot_dir}/{plot_prefix}_pvalue_qq.png")
    plt.close()
    
    # Plot 4: Distribution of phyloP scores by significance
    plt.figure(figsize=(12, 6))
    df_sig = df[df['significant'] == 'yes']['phyloP_score']
    df_not_sig = df[df['significant'] == 'no']['phyloP_score']
    plt.hist([df_sig, df_not_sig], bins=50, alpha=0.7, label=['Significant', 'Not significant'])
    plt.axvline(x=0, color='r', linestyle='--')
    plt.title("Distribution of phyloP Scores by Significance")
    plt.xlabel("phyloP Score")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(f"{plot_dir}/{plot_prefix}_score_by_significance.png")
    plt.close()
    
    logger.info(f"Diagnostic plots saved to {plot_dir}")

def print_summary(df):
    """Print summary statistics from results DataFrame."""
    total_positions = len(df)
    significant = df[df['significant'] == 'yes']
    constrained = significant[significant['status'] == 'constrained']
    accelerated = significant[significant['status'] == 'accelerated']
    
    logger.info("\nAnalysis summary:")
    logger.info(f"Total positions analyzed: {total_positions}")
    logger.info(f"Significant positions: {len(significant)} ({len(significant)/total_positions*100:.2f}%)")
    logger.info(f"  - Constrained: {len(constrained)} ({len(constrained)/total_positions*100:.2f}%)")
    logger.info(f"  - Accelerated: {len(accelerated)} ({len(accelerated)/total_positions*100:.2f}%)")
    
    # Дополнительная статистика
    status_counts = df['status'].value_counts()
    logger.info("\nStatus breakdown (before significance filtering):")
    for status, count in status_counts.items():
        logger.info(f"  - {status}: {count} ({count/total_positions*100:.2f}%)")
    
    # Выводим те же данные и на консоль
    print("\nAnalysis summary:")
    print(f"Total positions analyzed: {total_positions}")
    print(f"Significant positions: {len(significant)} ({len(significant)/total_positions*100:.2f}%)")
    print(f"  - Constrained: {len(constrained)} ({len(constrained)/total_positions*100:.2f}%)")
    print(f"  - Accelerated: {len(accelerated)} ({len(accelerated)/total_positions*100:.2f}%)")

def parallel_process_phylop_file(input_file, output_file, alpha=0.05, plot=False, n_processes=None, fdr_combined=False):
    """
    Обрабатывает файл с phyloP-скорами, используя параллельную обработку.
    
    Args:
        input_file (str): Путь к входному файлу с phyloP скорами
        output_file (str): Путь к выходному файлу с результатами
        alpha (float): Порог значимости для q-значений (по умолчанию: 0.05)
        plot (bool): Создавать ли диагностические графики (по умолчанию: False)
        n_processes (int): Количество процессов (по умолчанию: число ядер CPU)
        fdr_combined (bool): Применять ли FDR-коррекцию ко всем скорам независимо от знака (по умолчанию: False)
    """
    start_time = time.time()
    
    # Определяем количество процессов
    if n_processes is None:
        n_processes = min(mp.cpu_count(), 8)  # Ограничиваем 8 процессами по умолчанию
    
    logger.info(f"Processing phyloP file: {input_file} | Output: {output_file} | Using {n_processes} processes | FDR mode: {'combined' if fdr_combined else 'separate'}")
    
    # Фаза 1: Анализ структуры файла
    blocks, total_lines, data_lines = analyze_file_structure(input_file)
    
    # Фаза 2: Распределение блоков между процессами
    process_assignments = distribute_blocks(blocks, n_processes)
    
    # Фаза 3: Параллельная обработка блоков
    all_results = []
    
    if n_processes > 1:
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            futures = [executor.submit(process_blocks, input_file, assignments, i) 
                      for i, assignments in enumerate(process_assignments)]
            
            for i, future in enumerate(futures):
                try:
                    chunk_results = future.result()
                    all_results.extend(chunk_results)
                except Exception as e:
                    logger.error(f"Error processing chunk {i+1}: {e}")
    else:
        # Последовательная обработка для одного процесса
        for i, assignments in enumerate(process_assignments):
            chunk_results = process_blocks(input_file, assignments, i)
            all_results.extend(chunk_results)
    
    logger.info(f"Data extraction completed. Total positions: {len(all_results)} (expected: {data_lines}, coverage: {len(all_results) / data_lines * 100:.2f}%)")
    
    # Преобразуем результаты в DataFrame
    df = pd.DataFrame(all_results)
    
    if df.empty:
        logger.error("No valid data found in the file.")
        return df
    
    # Фаза 4: Статистический анализ
    df = apply_statistical_analysis(df, alpha=alpha, fdr_combined=fdr_combined)
    
    # Создаем диагностические графики, если требуется
    if plot:
        try:
            generate_plots(df, output_file)
        except Exception as e:
            logger.error(f"Warning: Could not generate plots: {e}")
    
    # Записываем результаты в выходной файл
    try:
        df.to_csv(output_file, sep='\t', index=False, na_rep="NA", 
                 float_format='%.6e')
        logger.info(f"Results written to {output_file}")
    except Exception as e:
        logger.error(f"Error writing output file: {e}")
    
    # Печатаем статистику
    print_summary(df)
    
    elapsed_time = time.time() - start_time
    logger.info(f"Total processing time: {elapsed_time:.2f} seconds | Peak memory usage: {psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024:.2f} MB")
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Process phyloP scores to identify significant regions.')
    parser.add_argument('input_file', help='Input phyloP scores file')
    parser.add_argument('output_file', help='Output results file')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance threshold for q-values (default: 0.05)')
    parser.add_argument('--plot', action='store_true', help='Generate diagnostic plots (default: False)')
    parser.add_argument('--processes', type=int, help='Number of processes to use (default: CPU count)')
    parser.add_argument('--fdr-combined', action='store_true', 
                      help='Calculate FDR correction on all scores regardless of sign (default: False)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Настройка уровня логирования
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    parallel_process_phylop_file(
        args.input_file, 
        args.output_file, 
        alpha=args.alpha, 
        plot=args.plot,
        n_processes=args.processes,
        fdr_combined=args.fdr_combined
    )

if __name__ == "__main__":
    main()