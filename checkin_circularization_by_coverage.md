# Getting the 20kb window and end (10kb upstream and 10kb downstream)
cd /media/lgbio-nas1/hectorromao/Anacardium/checking_circularization
mkdir -p {junction/chloroplast,junction/mitochondria}

for i in genomes/Chloroplast/*.fa; do
    seqkit subseq -r 1:10000 "$i" > "junction/chloroplast/$(basename "$i").temp1"
    seqkit subseq -r -10000:-1 "$i" > "junction/chloroplast/$(basename "$i").temp2"
    seqkit concat "junction/chloroplast/$(basename "$i").temp2" "junction/chloroplast/$(basename "$i").temp1" >  "junction/chloroplast/$(basename "$i")_20kb_junction.fasta"
done


for i in genomes/Mitochondria/*.fa; do
    seqkit subseq -r 1:10000 "$i" > "junction/mitochondria/$(basename "$i").temp1"
    seqkit subseq -r -10000:-1 "$i" > "junction/mitochondria/$(basename "$i").temp2"
    seqkit concat "junction/mitochondria/$(basename "$i").temp2" "junction/mitochondria/$(basename "$i").temp1" >  "junction/mitochondria/$(basename "$i")_20kb_junction.fasta"
done

rm junction/*/*temp*


# Aligning the reads against the genome

## Organizing the environtment

mkdir -p mapping/{chloroplast,mitochondria}
mkdir -p stats/{chloroplast,mitochondria}
mkdir -p tablet/{chloroplast,mitochondria}


REFDIR_mt=/media/lgbio-nas1/hectorromao/Anacardium/checking_circularization/junction/mitochondria
OUTDIR_mt=/media/lgbio-nas1/hectorromao/Anacardium/checking_circularization/mapping/mitochondria
FASTQ=/media/lgbio-nas1/renatadias/genoma-cajuzinho/2.Filtering_raw_reads/Ahu_trimmed_q20_l500.fastq.gz


REFDIR_pt=/media/lgbio-nas1/hectorromao/Anacardium/checking_circularization/junction/chloroplast
OUTDIR_pt=/media/lgbio-nas1/hectorromao/Anacardium/checking_circularization/mapping/chloroplast
FASTQ=/media/lgbio-nas1/renatadias/genoma-cajuzinho/2.Filtering_raw_reads/Ahu_trimmed_q20_l500.fastq.gz



## aligning

for ref in "$REFDIR_pt"/*.fasta; do

    name=$(basename "$ref" .fasta)

    echo "========================================"
    echo "Mapping: $name"
    echo "========================================"

    minimap2 \
        -ax map-ont \
        -t 48 \
        "$ref" \
        "$FASTQ" \
    | samtools sort \
        -@ 48 \
        -o "$OUTDIR_pt/${name}.bam" -

    samtools index \
        "$OUTDIR_pt/${name}.bam"

done


for ref in "$REFDIR_mt"/*.fasta; do

    name=$(basename "$ref" .fasta)

    echo "========================================"
    echo "Mapping: $name"
    echo "========================================"

    minimap2 \
        -ax map-ont \
        -t 48 \
        "$ref" \
        "$FASTQ" \
    | samtools sort \
        -@ 48 \
        -o "$OUTDIR_mt/${name}.bam" -

    samtools index \
        "$OUTDIR_mt/${name}.bam"

done

## How many reads were aligned against each genome

for bam in mapping/chloroplast/*.bam; do

    name=$(basename "$bam" .bam)

    samtools flagstat  -@ 24 "$bam" \
        > "stats/chloroplast/${name}.flagstat.txt"

done

for bam in mapping/mitochondria/*.bam; do

    name=$(basename "$bam" .bam)

    samtools flagstat -@ 24 "$bam" \
        > "stats/mitochondria/${name}.flagstat.txt"

done


## Getting the depth

for bam in mapping/chloroplast/*.bam; do

    name=$(basename "$bam" .bam)

    samtools depth  -@ 24 "$bam" \
        > "stats/chloroplast/${name}.depth.txt"

done

for bam in mapping/mitochondria/*.bam; do

    name=$(basename "$bam" .bam)

    samtools depth -@ 24 "$bam" \
        > "stats/mitochondria/${name}.depth.txt"

done


## Checking how many reads goes through the junction

### To chloroplast

declare -A CHLOROPLAST

for fasta in genomes/Chloroplast/*.fa; do
    name=$(basename "$fasta" .fa)
    contig=$(grep '^>' "$fasta" | sed 's/^>//' | awk '{print $1}')

    CHLOROPLAST["$name"]="$contig"
done


for name in "${!CHLOROPLAST[@]}"; do
    contig="${CHLOROPLAST[$name]}"

    echo "Processing: $name"
    echo "Contig: $contig"

    samtools view -F 0x900 \
        "mapping/chloroplast/${name}.fa_20kb_junction.bam" \
        "${contig}:9900-9999" |
        cut -f1 |
        sort -u > "/tmp/${name}_left.txt"

    samtools view -F 0x900 \
        "mapping/chloroplast/${name}.fa_20kb_junction.bam" \
        "${contig}:10001-10100" |
        cut -f1 |
        sort -u > "/tmp/${name}_right.txt"

    comm -12 \
        "/tmp/${name}_left.txt" \
        "/tmp/${name}_right.txt" |
        wc -l > "stats/chloroplast/${name}_reads_across_junction.txt"
done


### For mitochondria

declare -A MITOCHONDRIA

for fasta in genomes/Mitochondria/*.fa; do
    name=$(basename "$fasta" .fa)
    contig=$(grep '^>' "$fasta" | sed 's/^>//' | awk '{print $1}')

    MITOCHONDRIA["$name"]="$contig"
done

for fasta in genomes/Chloroplast/*.fa; do
    name=$(basename "$fasta" .fa)
    contig=$(grep '^>' "$fasta" | sed 's/^>//' | awk '{print $1}')
    MITOCHONDRIA["$name"]="$contig"
done

for name in "${!MITOCHONDRIA[@]}"; do
    contig="${MITOCHONDRIA[$name]}"

    echo "Processing: $name"
    echo "Contig: $contig"

    samtools view -F 0x900 \
        "mapping/mitochondria/${name}.fa_20kb_junction.bam" \
        "${contig}:9900-9999" |
        cut -f1 |
        sort -u > "/tmp/${name}_left.txt"

    samtools view -F 0x900 \
        "mapping/mitochondria/${name}.fa_20kb_junction.bam" \
        "${contig}:10001-10100" |
        cut -f1 |
        sort -u > "/tmp/${name}_right.txt"

    comm -12 \
        "/tmp/${name}_left.txt" \
        "/tmp/${name}_right.txt" |
        wc -l > "stats/mitochondria/${name}_reads_across_junction.txt"
done



## IGV graphics

#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# CONFIGURAÇÕES — ajuste só esses caminhos, se precisar
# ============================================================
BASE="/media/lgbio-nas1/hectorromao/Anacardium/checking_circularization"
JUNCTION_DIR="$BASE/junction"
MAPPING_DIR="$BASE/mapping"
SNAPSHOT_DIR="$BASE/snapshots"
BATCH_FILE="$BASE/igv_batch_all.bat"

mkdir -p "$SNAPSHOT_DIR"
: > "$BATCH_FILE"   # cria/zera o arquivo de batch

# ============================================================
# FUNÇÃO: escreve um bloco de comandos IGV para UMA montagem
# (edite aqui as preferências/zooms — a mudança vale pra
#  todas as combinações automaticamente)
# ============================================================
add_block () {
    local organelle="$1"
    local fasta="$2"
    local bam="$3"
    local label
    label=$(basename "$fasta" .fasta)

    cat >> "$BATCH_FILE" <<EOF
new
genome $fasta
load $bam
snapshotDirectory $SNAPSHOT_DIR
maxPanelHeight 50000
preference SAM.DOWNSAMPLE_READS false
preference SAM.SHOW_COV_TRACK true
preference SAM.SHOW_ALIGNMENT_TRACK true
preference SAM.SHOW_INSERTION_MARKERS false
preference SAM.SHOW_CENTER_LINE true
preference SAM.SHOW_SOFT_CLIPPED true
preference SAM.HIDE_SMALL_INDEL true
preference SAM.SMALL_INDEL_BP_THRESHOLD 3
expand
colorBy MAPPED_SIZE
goto 9000-11000
sort POSITION 10000
snapshot ${organelle}_${label}_2kb.png
goto 9800-10200
sort POSITION 10000
snapshot ${organelle}_${label}_400bp.png
EOF
}

# ============================================================
# LOOP: para cada organela (subpasta de mapping/), casa cada
# BAM presente com o FASTA de mesmo nome em junction/
# ============================================================
n_ok=0
n_skip=0

for organelle_dir in "$MAPPING_DIR"/*/; do
    organelle=$(basename "$organelle_dir")

    for bam in "$organelle_dir"*.bam; do
        [ -e "$bam" ] || continue   # pasta sem nenhum .bam, ignora

        base=$(basename "$bam" .bam)
        fasta="$JUNCTION_DIR/$organelle/${base}.fasta"

        if [ ! -f "$fasta" ]; then
            echo "AVISO: fasta não encontrado para $bam" >&2
            echo "       esperado em: $fasta" >&2
            n_skip=$((n_skip + 1))
            continue
        fi

        echo "  [$organelle] $base"
        add_block "$organelle" "$fasta" "$bam"
        n_ok=$((n_ok + 1))
    done
done

echo "exit" >> "$BATCH_FILE"

echo ""
echo "Blocos adicionados: $n_ok"
echo "Combinações puladas (sem fasta correspondente): $n_skip"
echo "Batch gerado em: $BATCH_FILE"
echo ""
echo "Para rodar:"
echo "  igv.sh -b $BATCH_FILE"



xvfb-run --auto-servernum IGV_Linux_2.19.2/igv.sh -b