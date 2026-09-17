# Circularization Checkup

This repository is dedicated to evaluating alternative genome assemblies and selecting the one that provides the strongest evidence of circularization. The workflow compares candidate assemblies for plastid and mitochondrial genomes by measuring whether sequencing reads support the contig ends and the reconstructed circular junction.

The main goal is to move beyond visual inspection and produce quantitative metrics that help answer a practical question:

> Which assembly is the most convincing circularized representation of the genome?

## Project objective

The analyses in this repository generate metrics for the comparison of assembly candidates, including:

- read mapping depth across the assembled sequence;
- overall alignment quality via flagstat summaries;
- number of reads spanning the circularization junction;
- visual evidence from IGV-style snapshot images;
- comparison of junction support between different assembly strategies.

This makes the repository useful for benchmarking assemblies produced with different tools or parameter sets (for example, Flye, HiFiasm, or reference-guided reassembly), and for choosing the best assembly based on reproducible evidence.

## Main workflow

The pipeline follows a standard comparison strategy:

1. Build synthetic junction sequences by extracting the terminal 10 kb from each end of each assembly and concatenating them into a 20 kb circularization window.
2. Map long reads to each junction sequence.
3. Summarize alignment statistics with samtools.
4. Measure read depth and support at the putative circularization boundary.
5. Count reads that map to both sides of the junction.
6. Generate visual snapshots to confirm the signal in a genome browser.

## Repository structure

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

## What the outputs represent

### `junction/`
Contains the artificial 20 kb circularization sequences built around the ends of each assembly.

### `mapping/`
Stores BAM files produced by mapping reads against each candidate junction sequence.

### `stats/`
Contains quantitative summaries such as:

- read counts;
- depth distributions;
- flagstat quality metrics;
- counts of reads crossing the circularization junction.

### `snapshots/`
Stores visual views around the junction region for each assembly. These images help confirm whether the supporting reads cluster at the expected circularization site.

## Expected use

This project is intended for exploratory and comparative assembly validation. It is especially useful when the objective is to identify which reconstruction presents the most consistent evidence of a true circular molecule.

In practice, an assembly is considered more promising when it shows:

- strong, consistent read coverage around the junction;
- a substantial number of reads spanning both sides of the junction;
- clean mapping behavior in the region with no obvious breakpoints or low-coverage gaps;
- clear visual support in the junction image.

## Typical interpretation

A candidate with higher junction-spanning read support and more stable coverage around the circularization boundary is usually a stronger candidate for the final circularized assembly. However, the final decision should consider both the numerical metrics and the visual evidence, because biological or technical artifacts can sometimes inflate one metric without representing a true circularization event.

## Notes

- This repository is designed for reproducible assembly evaluation.
- Most results are generated from long-read data and low-level genome comparison scripts.
- The workflow is intended to support model selection and decision-making rather than to replace biological validation.

## License

This project is intended for internal analysis and research use. If you plan to share it more broadly, please document the source of the data and the corresponding analysis context.
