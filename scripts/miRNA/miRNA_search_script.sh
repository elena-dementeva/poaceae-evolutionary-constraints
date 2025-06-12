#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 -i <input_genome.fa>"
    echo "Options:"
    echo "  -i    Input genome file in FASTA format"
    exit 1
}

while getopts ":i:" opt; do
    case $opt in
        i)
            INPUT_GENOME="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

if [ -z "${INPUT_GENOME:-}" ]; then
    echo "Error: Input genome file not specified!" >&2
    usage
fi

check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 not found!" >&2
        exit 1
    fi
}

check_cmd() {
    if [ $? -ne 0 ]; then
        echo "Error: Command failed: $1" >&2
        exit 1
    fi
}

MAX_THREADS=$(sysctl -n hw.ncpu)

# Расчет свободной памяти (в GB) для macOS
FREE_MEM_GB=$(vm_stat | awk '/Pages free/ {free=$3; gsub(/\./, "", free)} /Pages speculative/ {spec=$3; gsub(/\./, "", spec)} END {total=(free + spec) * 4096 / 1024 / 1024 / 1024; printf "%.0f", total}')

SAFE_THREADS=$((FREE_MEM_GB / 2))  # 2GB на поток
THREADS=$(( SAFE_THREADS < MAX_THREADS ? SAFE_THREADS : MAX_THREADS-1 ))
THREADS=$(( THREADS < 1 ? 1 : THREADS ))  # Минимум 1 поток

echo "Using $THREADS threads (available: $MAX_THREADS cores, $FREE_MEM_GB GB free RAM)"

INPUT_GENOME="Chr1A_candidates_1.fa"
HAIRPIN_FILE="hairpin.fa"
MATURE_FILE="mature.fa"
PLANT_LIST="Acacia auriculiformis|Acacia mangium|Aegilops tauschii|Amborella trichopoda|Aquilegia caerulea|Arabidopsis lyrata|Arabidopsis thaliana|Arachis hypogaea|Asparagus officinalis|Avicennia marina|Brachypodium distachyon|Brassica napus|Brassica oleracea|Brassica rapa|Bruguiera cylindrica|Bruguiera gymnorhiza|Camelina sativa|Carica papaya|Chlamydomonas reinhardtii|Citrus clementina|Citrus reticulata|Citrus sinensis|Citrus trifoliata|Cucumis melo|Cucumis sativus|Cunninghamia lanceolata|Cynara cardunculus|Digitalis purpurea|Elaeis guineensis|Eugenia uniflora|Festuca arundinacea|Fragaria vesca|Glycine max|Glycine soja|Gossypium arboreum|Gossypium herbaceum|Gossypium hirsutum|Gossypium raimondii|Helianthus annuus|Helianthus argophyllus|Helianthus ciliaris|Helianthus exilis|Helianthus paradoxus|Helianthus petiolaris|Helianthus tuberosus|Hevea brasiliensis|Hordeum vulgare|Linum usitatissimum|Lotus japonicus|Malus domestica|Manihot esculenta|Medicago truncatula|Nicotiana tabacum|Oryza sativa|Paeonia lactiflora|Panax ginseng|Phaseolus vulgaris|Physcomitrella patens|Picea abies|Pinus densata|Pinus taeda|Populus euphratica|Populus trichocarpa|Prunus persica|Rehmannia glutinosa|Ricinus communis|Saccharum officinarum|Saccharum sp|Salicornia europaea|Salvia miltiorrhiza|Salvia sclarea|Selaginella moellendorffii|Solanum lycopersicum|Solanum tuberosum|Sorghum bicolor|Theobroma cacao|Triticum aestivum|Triticum turgidum|Vigna unguiculata|Vitis vinifera|Vriesea carinata|Zea mays"

check_file "$INPUT_GENOME"
check_file "$MATURE_FILE"
check_file "$HAIRPIN_FILE"
check_file "miRNA_filtrator.py"
check_file "100_identity_filtrator.py"

eval "$(conda shell.bash hook)"
conda activate mirna_analysis

conda install -y ViennaRNA blast parallel
check_cmd "conda install"

echo "Step 1: Processing FASTA entries with RNAfold"
check_file "$INPUT_GENOME"

time awk '/^>/ { if (seq) print seq; printf("%s\n", $0); seq=""; next } { seq = seq $0 } END { if (seq) print seq }' "$INPUT_GENOME" | \
parallel --jobs $THREADS --pipe -N 2 "RNAfold" > hairpin_structures.txt

check_cmd "RNAfold"

echo "Step 2: Move auxiliary files to a separate folder"
mkdir -p RNAfold_data
find . -maxdepth 1 -name "*_ss.ps" -print0 | xargs -0 -I {} mv {} "RNAfold_data/"
check_cmd "Moving auxiliary files"

echo "Step 3: Filtering sequences by free energy (<-20kcalmol^-1): " 
check_file "hairpin_structures.txt"
python miRNA_filtrator.py -i hairpin_structures.txt -o filtered_hairpin.fa
check_cmd "miRNA_filtrator.py"

echo "Step 4: Selecting plant hairpin sequences"
check_file "$HAIRPIN_FILE"

time awk -v pattern="$PLANT_LIST" '
    BEGIN { RS = ">"; FS = "\n"; ORS = "" } 
    $1 ~ pattern { print ">" $0 }
' "$HAIRPIN_FILE" > plant_hairpin.fa

check_cmd "awk plant sequences"
check_file "plant_hairpin.fa"

echo "Step 5: Creating a database of plant hairpin sequences"

makeblastdb -in plant_hairpin.fa -dbtype nucl -out miRbase_plants_hairpin -parse_seqids
check_cmd "makeblastdb hairpin"

echo "Step 6: Searching for miRNA hits in the hairpin database"
time blastn -query filtered_hairpin.fa -db miRbase_plants_hairpin -outfmt 6 \
     -evalue 1e-5 -word_size 4 -gapopen 5 -gapextend 2 \
     -num_threads $THREADS > hairpin_plant_miRNA_hits.txt
check_cmd "blastn search"
check_file "hairpin_plant_miRNA_hits.txt"

echo “Total found hairpin sequences:”

wc -l hairpin_plant_miRNA_hits.txt

echo "Step 7: Filtering sequences with % Identity > 90:" 
python 100_identity_filtrator.py -i hairpin_plant_miRNA_hits.txt -o filtered_hairpin.fa


check_cmd "100_identity_filtrator.py"
check_file "filtered_hairpin.fa_unique_list.txt"

echo “Total unique hairpin structures:”
wc -l filtered_hairpin.fa_unique_list.txt

echo "Step 8: Selecting plant mature miRNA sequences"
check_file "$MATURE_FILE"

time awk -v pattern="$PLANT_LIST" '
    BEGIN { RS = ">"; FS = "\n"; ORS = "" } 
    $1 ~ pattern { print ">" $0 }
' "$MATURE_FILE" > full_plant_mature.fa

check_cmd "awk mature sequences"
check_file "full_plant_mature.fa"

echo "Step 9: Creating a database of mature sequences:" 
makeblastdb -in full_plant_mature.fa -dbtype nucl -out miRbase_plants_mature_full -parse_seqids
check_cmd "makeblastdb mature"

echo "Step 10: Extracting mature miRNA sequences from the database"
check_file "filtered_hairpin.fa_unique_list.txt"
blastdbcmd -db miRbase_plants_mature_full -entry_batch filtered_hairpin.fa_unique_list.txt \
           -out mature_unique_miRNA_sequence.fa
check_cmd "blastdbcmd"
check_file "mature_unique_miRNA_sequence.fa"

echo "Total unique mature miRNA sequences: "
wc -l mature_unique_miRNA_sequences_.fa

echo "Pipeline completed successfully in $SECONDS seconds!"