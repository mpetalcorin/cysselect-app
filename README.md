# CysSelect

CysSelect is a browser-based Streamlit app for chemoproteomic hit calling, enantiomer comparison, prioritization, simulation, and next-library design for covalent fragment discovery.
<img width="1537" height="1017" alt="Screenshot 2026-04-18 at 04 46 45" src="https://github.com/user-attachments/assets/256615d4-7759-4a27-b44b-de5face28e59" />
## Overview

CysSelect helps researchers move from raw chemoproteomic screening tables to interpretable, decision-oriented outputs. The app is designed for research, method development, and educational use, with support for both real screening datasets and simulation-based benchmarking.

The current workflow is organized into three main modes:

- **Analyze Real Data**, for upload, quality control, hit calling, enantiomer comparison, and pair-level prioritization
- **Simulation and Benchmarking**, for generating synthetic datasets and stress-testing thresholds and workflow behavior
- **Library Design**, for scoring virtual compounds and selecting the next batch for follow-up screening

## Main Features

- Upload and quality-check chemoproteomic screening data
- Calculate core screening metrics, including CR, percent competition, AUC, and ŒîAUC
- Compare R and S enantiomer behavior at the cysteine and pair levels
- Summarize pair-level promiscuity and stereoselectivity
- Identify sweet-spot fragment pairs
- Generate simulated benchmark datasets with configurable noise and missingness
- Score virtual libraries for next-round compound selection
- Export processed results as CSV files

## App Structure

### Analyze Real Data
This mode is intended for routine analysis of uploaded screening tables.

Tabs include:
- Home
- Upload and QC
- Hit Calling
- Enantiomer Comparison
- Pair-Level Summary
- Hit Prioritization
- Export

### Simulation and Benchmarking
This mode is intended for method development, teaching, and threshold stress-testing.

### Library Design
This mode is intended for ranking candidate compounds and selecting a next-round screening batch.

## Expected Input

Required columns:

- `pair_id`
- `enantiomer`
- `cysteine_id`
- `concentration_uM`
- `replicate`
- `intensity_dmso`
- `intensity_condition`

Optional columns:

- `p_value`
- `protein_id`
- `gene_name`
- `scaffold_cluster`

## Core Concepts

### Competition Ratio
CysSelect calculates the competition ratio as:

`CR = intensity_dmso / intensity_condition`

### Percent Competition
Percent competition is calculated as:

`%comp = (1 - 1 / CR) * 100`

### Enantiomer Comparison
For two-concentration screening tables, the app estimates enantiomer-specific response summaries and computes ŒîAUC to compare R and S behavior.

### Pair-Level Summary
For each enantiopair, the app summarizes:
- Promiscuity
- Stereoselectivity fraction
- Sweet-spot status

## Running Locally

### 1. Create or activate an environment
```bash
conda create -n cysselect python=3.11 -y
conda activate cysselect
```

### 2. Install dependencies
```bash
pip install -r requirements.txt
```

## Requirements

Minimum dependencies:

- streamlit
- pandas
- numpy
- plotly

## Repository Files

Suggested repository contents:

- `app.py`
- `requirements.txt`
- `README.md`
- `modelcard.md`
- `datasheet.md`
- `.gitignore`

## Intended Use

CysSelect is intended for:
- chemoproteomic screening interpretation
- enantiomer-aware hit prioritization
- benchmark simulation
- early-stage library design support
- teaching and demonstration of chemoproteomic workflows

## Limitations

- The current implementation is intended as a research and prototyping tool
- Outputs should be interpreted alongside experimental validation
- Threshold-based hit calling depends on user settings and data quality
- Virtual library scoring is a prioritization aid, not a substitute for medicinal chemistry judgment

## Disclaimer

CysSelect is intended for research, method development, and educational use. App outputs are computational summaries and prioritization aids and should be interpreted alongside experimental validation and scientific judgment.

## License
MIT

## Citation
**Petalcorin, M.I.R.** (2026). CysSelect: a browser-based platform for chemoproteomic hit calling, enantiomer comparison, simulation, and next-library design for covalent fragment discovery. 
https://github.com/mpetalcorin/cysselect-app
