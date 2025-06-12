def calculate_coverage(maf_file, output_file):
    coverage_dict = {}
    ref_name = None

    with open(maf_file, "r") as f:
        block = []
        for line in f:
            if line.startswith("s"):
                parts = line.split()
                chrom = parts[1]
                start = int(parts[2])
                seq = parts[6]
                if ref_name is None:
                    ref_name = chrom
                block.append((chrom, start, seq))

            elif line.strip() == "" and block:
                for chrom, start, seq in block:
                    for pos in range(len(seq)):
                        if seq[pos] != "-":
                            coord = start + pos
                            coverage_dict[(chrom, coord)] = (
                                coverage_dict.get((chrom, coord), 0) + 1
                            )
                block = []

    with open(output_file, "w") as out:
        for (chrom, position), coverage in sorted(coverage_dict.items()):
            out.write(f"{chrom}\t{position}\t{position+1}\t{coverage}\n")


calculate_coverage("chrM.maf", "coverage_positions.txt")
