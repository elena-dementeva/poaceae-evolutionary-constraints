#!/bin/bash

while IFS=$'\t' read -r accession col2 species col4; do
    # replace space in species names
    species_clean=$(echo "$species" | tr ' ' '_')

    echo "Downloading $accession for species $species_clean..."
    output_dir="/home/labs/alevy/petrzhu/Wheat/Assemblies/Download/${accession}_${species_clean}"

    ~/Prog/ncbi_datasets/datasets download genome accession "$accession" --include gff3,rna,cds,protein,genome,seq-report --filename "$output_dir"
done < accessions_redo.txt
