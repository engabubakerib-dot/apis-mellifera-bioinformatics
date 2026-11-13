# Apis mellifera Stress Transcriptomics: Vitellogenin (Vg) and Major Royal Jelly Proteins (MRJPs)

## Project Overview

This repository contains source code developed for the bioinformatics analysis of *Apis mellifera* stress-related transcriptomic data, with particular focus on vitellogenin (Vg), major royal jelly proteins (MRJPs), and associated molecular pathways.

The project includes scripts for differential expression analysis, gene mapping, transcriptomic data processing, literature analysis, and biological interpretation.

## Repository Contents

- `scripts/R/` – R scripts for transcriptomic analysis and supporting workflows.
- `scripts/differential_expression/` – Differential expression re-analysis scripts.
- `scripts/download/` – Scripts for retrieving public datasets.
- `scripts/gene_mapping/` – Gene mapping and ortholog verification.
- `scripts/interpretation/` – Scripts for biological interpretation and analysis of Vg/MRJP-associated patterns.
- `scripts/literature/` – Literature processing and analysis utilities.
- `scripts/meta_analysis/` – Scripts for meta-analysis and figure generation.
- `scripts/rawreads/` – Raw-read processing and quantification workflows.
- `scripts/report/` – Automated report-generation code.
- `scripts/validation/` – Validation and analytical assessment scripts.
- `tools/` – Supporting third-party tool documentation and licenses.

## Data Availability

This repository focuses on source code. Raw sequencing data, intermediate files, generated results, and manuscript materials are not included in this public repository.

Where applicable, the analysis scripts refer to publicly available datasets and external resources. Required input files and software dependencies must be obtained and configured separately.

## Software and Dependencies

The project uses Python, R, and supporting command-line tools. The `environment.yml` file provides the Conda environment specification for the Python components.

Additional R and Bioconductor dependencies may be required for specific workflows.

## Scope and Reproducibility

The repository is provided to support methodological transparency and inspection of the computational workflows associated with this research.

It is not, by itself, a complete standalone reproducibility package because the underlying datasets, intermediate files, and generated outputs are not distributed here.

## License

Please refer to the license information accompanying the relevant third-party tools. The licensing terms for the project source code should be specified by the authors before reuse.
