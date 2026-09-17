# Circularization checkup by coverage

This workflow was designed to evaluate alternative genome assemblies and to estimate which assembly provides the strongest evidence of circularization. The analysis compares candidate assemblies for chloroplast and mitochondrial genomes by measuring support for the junction that connects the ends of the circular molecule.

The main idea is simple: if a genome is truly circular, reads should support the junction between the end and the beginning of the sequence. By measuring that support with robust metrics, we can compare assemblies and select the one that best fits the expected circular structure.

## Overview

The workflow does the following:

1. Extracts the terminal 10 kb from both ends of each assembled contig.
2. Concatenates them to generate a synthetic 20 kb circularization region.
3. Maps long reads against each candidate junction sequence.
4. Summarizes alignment quality using `samtools`.
5. Counts how many reads span the junction.
6. Produces IGV snapshots to visually confirm the signal around the junction.

---

## Requirements

The pipeline expects the following tools to be installed and available in the environment:

- `seqkit`
- `minimap2`
- `samtools`
- `IGV`
- `xvfb-run` (when running IGV in a headless Linux session)

---

## Project structure

```text
.
├── README.md
├── checkin_circularization_by_coverage.md
├── genomes/
│   ├── Chloroplast/
│   └── Mitochondria/
├── junction/
│   ├── chloroplast/
│   └── mitochondria/
├── mapping/
│   ├── chloroplast/
│   └── mitochondria/
├── snapshots/
├── stats/
│   ├── chloroplast/
│   └── mitochondria/
└── .git/
```

---

## 1) Build the 20 kb circularization window

The first step creates a synthetic 20 kb region by taking the last 10 kb of the sequence, the first 10 kb of the same sequence, and joining them in the order expected for a circular molecule.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="$(pwd)"
mkdir -p "$BASE/junction/chloroplast" "$BASE/junction/mitochondria"

# Build the 20 kb circularization window for chloroplast assemblies.
for fasta in "$BASE"/genomes/Chloroplast/*.fa; do
    sample_name=$(basename "$fasta" .fa)

    # Extract the first 10 kb and the last 10 kb.
    seqkit subseq -r 1:10000 "$fasta" > "$BASE/junction/chloroplast/${sample_name}.temp1"
    seqkit subseq -r -10000:-1 "$fasta" > "$BASE/junction/chloroplast/${sample_name}.temp2"

    # Concatenate the two ends to recreate a circular-junction representation.
    seqkit concat \
        "$BASE/junction/chloroplast/${sample_name}.temp2" \
        "$BASE/junction/chloroplast/${sample_name}.temp1" \
        > "$BASE/junction/chloroplast/${sample_name}_20kb_junction.fasta"
done

# Build the 20 kb circularization window for mitochondrial assemblies.
for fasta in "$BASE"/genomes/Mitochondria/*.fa; do
    sample_name=$(basename "$fasta" .fa)

    seqkit subseq -r 1:10000 "$fasta" > "$BASE/junction/mitochondria/${sample_name}.temp1"
    seqkit subseq -r -10000:-1 "$fasta" > "$BASE/junction/mitochondria/${sample_name}.temp2"

    seqkit concat \
        "$BASE/junction/mitochondria/${sample_name}.temp2" \
        "$BASE/junction/mitochondria/${sample_name}.temp1" \
        > "$BASE/junction/mitochondria/${sample_name}_20kb_junction.fasta"
done

# Remove temporary files after the final junction FASTA has been generated.
rm -f "$BASE"/junction/*/*.temp1 "$BASE"/junction/*/*.temp2
```

---

## 2) Align reads against each junction sequence

The next step maps the read set against each circularization candidate and stores the results in `mapping/`.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="$(pwd)"
FASTQ="/path/to/reads.fastq.gz"

mkdir -p "$BASE/mapping/chloroplast" "$BASE/mapping/mitochondria"
mkdir -p "$BASE/stats/chloroplast" "$BASE/stats/mitochondria"

# Map reads against chloroplast junction candidates.
for ref in "$BASE"/junction/chloroplast/*.fasta; do
    sample_name=$(basename "$ref" .fasta)

    echo "========================================"
    echo "Mapping chloroplast: $sample_name"
    echo "========================================"

    minimap2 \
        -ax map-ont \
        -t 48 \
        "$ref" \
        "$FASTQ" \
    | samtools sort \
        -@ 48 \
        -o "$BASE/mapping/chloroplast/${sample_name}.bam" -

    samtools index "$BASE/mapping/chloroplast/${sample_name}.bam"
done

# Map reads against mitochondrial junction candidates.
for ref in "$BASE"/junction/mitochondria/*.fasta; do
    sample_name=$(basename "$ref" .fasta)

    echo "========================================"
    echo "Mapping mitochondria: $sample_name"
    echo "========================================"

    minimap2 \
        -ax map-ont \
        -t 48 \
        "$ref" \
        "$FASTQ" \
    | samtools sort \
        -@ 48 \
        -o "$BASE/mapping/mitochondria/${sample_name}.bam" -

    samtools index "$BASE/mapping/mitochondria/${sample_name}.bam"
done
```

---

## 3) Collect BAM summary statistics

This step stores alignment metrics that describe how many reads map confidently and how much coverage the reads provide to each candidate assembly.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="$(pwd)"

# Compute flagstat summaries for chloroplast mappings.
for bam in "$BASE"/mapping/chloroplast/*.bam; do
    sample_name=$(basename "$bam" .bam)
    samtools flagstat -@ 24 "$bam" > "$BASE/stats/chloroplast/${sample_name}.flagstat.txt"
done

# Compute flagstat summaries for mitochondrial mappings.
for bam in "$BASE"/mapping/mitochondria/*.bam; do
    sample_name=$(basename "$bam" .bam)
    samtools flagstat -@ 24 "$bam" > "$BASE/stats/mitochondria/${sample_name}.flagstat.txt"
done

# Compute depth profiles for chloroplast mappings.
for bam in "$BASE"/mapping/chloroplast/*.bam; do
    sample_name=$(basename "$bam" .bam)
    samtools depth -@ 24 "$bam" > "$BASE/stats/chloroplast/${sample_name}.depth.txt"
done

# Compute depth profiles for mitochondrial mappings.
for bam in "$BASE"/mapping/mitochondria/*.bam; do
    sample_name=$(basename "$bam" .bam)
    samtools depth -@ 24 "$bam" > "$BASE/stats/mitochondria/${sample_name}.depth.txt"
done
```

---

## 4) Count reads spanning the circularization junction

The most informative metric for circularization is the number of reads that map to both sides of the synthetic junction. A stronger circularized assembly typically shows a larger number of reads spanning the boundary.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="$(pwd)"

# Define a helper function to count reads spanning the junction.
count_reads_across_junction() {
    local organelle="$1"
    local sample_name="$2"
    local contig="$3"

    # Reads mapping just before the junction.
    samtools view -F 0x900 \
        "$BASE/mapping/${organelle}/${sample_name}.bam" \
        "${contig}:9900-9999" \
    | cut -f1 \
    | sort -u > "/tmp/${sample_name}_left.txt"

    # Reads mapping just after the junction.
    samtools view -F 0x900 \
        "$BASE/mapping/${organelle}/${sample_name}.bam" \
        "${contig}:10001-10100" \
    | cut -f1 \
    | sort -u > "/tmp/${sample_name}_right.txt"

    # Keep only reads detected in both regions.
    comm -12 \
        "/tmp/${sample_name}_left.txt" \
        "/tmp/${sample_name}_right.txt" \
    | wc -l > "$BASE/stats/${organelle}/${sample_name}_reads_across_junction.txt"
}

# Evaluate chloroplast assemblies.
declare -A CHLOROPLAST
for fasta in "$BASE"/genomes/Chloroplast/*.fa; do
    sample_name=$(basename "$fasta" .fa)
    contig=$(grep '^>' "$fasta" | sed 's/^>//' | awk '{print $1}')
    CHLOROPLAST["$sample_name"]="$contig"
done

for sample_name in "${!CHLOROPLAST[@]}"; do
    contig="${CHLOROPLAST[$sample_name]}"
    echo "Processing chloroplast: $sample_name"
    count_reads_across_junction "chloroplast" "$sample_name" "$contig"
done

# Evaluate mitochondrial assemblies.
declare -A MITOCHONDRIA
for fasta in "$BASE"/genomes/Mitochondria/*.fa; do
    sample_name=$(basename "$fasta" .fa)
    contig=$(grep '^>' "$fasta" | sed 's/^>//' | awk '{print $1}')
    MITOCHONDRIA["$sample_name"]="$contig"
done

for sample_name in "${!MITOCHONDRIA[@]}"; do
    contig="${MITOCHONDRIA[$sample_name]}"
    echo "Processing mitochondria: $sample_name"
    count_reads_across_junction "mitochondria" "$sample_name" "$contig"
done
```

---

## 5) Generate IGV snapshots

This section produces visual evidence around the putative circularization site. The final goal is to inspect read support around the junction in a genome browser and compare the different assemblies side by side.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="$(pwd)"
JUNCTION_DIR="$BASE/junction"
MAPPING_DIR="$BASE/mapping"
SNAPSHOT_DIR="$BASE/snapshots"
BATCH_FILE="$BASE/igv_batch_all.bat"

mkdir -p "$SNAPSHOT_DIR"
: > "$BATCH_FILE"

add_block() {
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

n_ok=0
n_skip=0

for organelle_dir in "$MAPPING_DIR"/*/; do
    organelle=$(basename "$organelle_dir")

    for bam in "$organelle_dir"*.bam; do
        [ -e "$bam" ] || continue

        sample_name=$(basename "$bam" .bam)
        fasta="$JUNCTION_DIR/$organelle/${sample_name}.fasta"

        if [ ! -f "$fasta" ]; then
            echo "Warning: FASTA not found for $bam" >&2
            echo "Expected path: $fasta" >&2
            n_skip=$((n_skip + 1))
            continue
        fi

        echo "[$organelle] $sample_name"
        add_block "$organelle" "$fasta" "$bam"
        n_ok=$((n_ok + 1))
    done
done

echo "exit" >> "$BATCH_FILE"

echo "Blocks added: $n_ok"
echo "Missing FASTA files: $n_skip"
echo "Batch file generated at: $BATCH_FILE"
echo "Run with:"
echo "  igv.sh -b $BATCH_FILE"

# Example for headless execution:
# xvfb-run --auto-servernum IGV_Linux_2.19.2/igv.sh -b "$BATCH_FILE"
```

---

## Interpretation of the results

The best assembly is not necessarily the one with the longest sequence or the highest total read coverage. In this project, the assembly is considered more promising when it shows:

- strong read support at the circularization boundary;
- high coverage across the junction region;
- many reads spanning both sides of the junction;
- consistent signal in the IGV snapshots;
- no obvious read dropout or abnormal coverage interruption near the junction.

This repository is meant to support assembly selection with quantitative evidence, not only visual inspection.

---

## Final note

This project is a comparative framework for assembly evaluation. It is especially useful when the goal is to choose the most convincing circularized genome reconstruction from multiple candidates produced with different assembly methods or parameter sets.
