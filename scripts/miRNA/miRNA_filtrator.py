import re
import argparse

def parse_rnafold_output(filename):
    with open(filename) as f:
        data = f.read().split('>')[1:]
    candidates = []
    for record in data:
        lines = record.strip().split('\n')
        if len(lines) < 2:
            continue
        header = lines[0]
        seq = lines[1]
        # Ищем dG в строке типа ".... ( -3.10)"
        dG_match = re.search(r"\(([-\s]\d+\.\d+)\)", lines[-1])  # Ищем числа в скобках
        if not dG_match:
            print(f"Ошибка: не найдено dG в строке: {lines[-1]}")
            continue
        dG = float(dG_match.group(1).strip())  # Конвертируем в float
        candidates.append((header, seq, dG))
    return candidates

def main():
    parser = argparse.ArgumentParser(description="Фильтрация miRNA кандидатов по dG и длине.")
    parser.add_argument("-i", "--input", required=True, help="Входной файл (результат RNAfold)")
    parser.add_argument("-o", "--output", required=True, help="Выходной файл (отфильтрованные кандидаты)")
    args = parser.parse_args()

    candidates = parse_rnafold_output(args.input)

    with open(args.output, "w") as out:
        for header, seq, dG in candidates:
            if dG < -20 and 60 <= len(seq):
                out.write(f">{header}\n{seq}\n")

if __name__ == "__main__":
    main()