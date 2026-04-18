# Model Card: CysSelect

## Model Name

CysSelect

## Model Type

Browser-based chemoproteomic analysis and prioritization application with heuristic and simulation-supported scoring components.

## Version

Version 1.0, initial public GitHub release

## Summary

CysSelect is a Streamlit-based research app designed to analyze covalent fragment chemoproteomic screening data. It supports hit calling, enantiomer comparison, pair-level prioritization, simulation-based benchmarking, and virtual library ranking for next-round compound selection.

The app is not a single monolithic predictive model. Instead, it combines:
- deterministic metric calculations
- threshold-based classification logic
- pair-level aggregation
- simulation-based synthetic data generation
- heuristic prioritization and acquisition scoring

## Intended Users

- chemical biologists
- proteomics researchers
- covalent fragment discovery teams
- medicinal chemists
- computational biologists
- students and trainees in chemoproteomics

## Intended Use

CysSelect is intended for:
- exploratory analysis of chemoproteomic screening tables
- enantiomer-aware prioritization of fragment pairs
- simulation-based workflow testing
- preliminary virtual library triage
- educational demonstration of covalent fragment discovery workflows

## Out of Scope Use

CysSelect is not intended for:
- clinical decision-making
- patient diagnosis
- regulatory decision support
- autonomous medicinal chemistry optimization
- use as a substitute for experimental validation

## Inputs

### Required Inputs
- `pair_id`
- `enantiomer`
- `cysteine_id`
- `concentration_uM`
- `replicate`
- `intensity_dmso`
- `intensity_condition`

### Optional Inputs
- `p_value`
- `protein_id`
- `gene_name`
- `scaffold_cluster`

### Simulation Inputs
Simulation mode allows user-defined parameters such as:
- number of enantiopairs
- number of cysteine sites
- assay noise
- missing-data rate
- random seed

## Outputs

CysSelect produces:
- competition ratio
- percent competition
- enantiomer-level response summaries
- ŒîAUC
- liganded hit labels
- enantioselective class labels
- pair-level promiscuity
- pair-level stereoselectivity fraction
- sweet-spot status
- ranked priority tables
- simulated benchmark datasets
- virtual library scores
- selected next-round batches

## Core Decision Logic

### Metric Layer
The application computes CR and derived quantities from uploaded intensity tables.

### Comparison Layer
The application compares R and S enantiomer behavior at the site level and summarizes differences using ŒîAUC.

### Aggregation Layer
The application summarizes pair-level promiscuity and stereoselectivity.

### Prioritization Layer
The application ranks pairs or candidate compounds using sweet-spot logic and acquisition scoring.

## Factors Influencing Outputs

Outputs depend on:
- data quality
- missingness
- threshold settings
- replicate consistency
- simulation settings
- heuristic scoring design
- optional metadata quality

## Strengths

- easy browser-based workflow
- integrates upload, QC, hit calling, prioritization, simulation, and design
- emphasizes enantiomer-resolved interpretation
- useful for exploratory and educational workflows
- deployable without specialized infrastructure

## Limitations

- threshold-based and heuristic components may be sensitive to user settings
- simulation outputs are synthetic and depend on model assumptions
- current virtual library scoring is a prioritization aid rather than a fully validated production predictor
- performance depends on the structure and quality of the input table
- no claim is made that outputs generalize to all chemoproteomic platforms without adaptation

## Risks

Potential risks include:
- overinterpretation of weak or noisy hits
- misuse of outputs as definitive evidence without validation
- threshold choices that are too permissive or too strict
- misformatted uploads leading to misleading summaries

## Mitigations

Users should:
- inspect uploaded data carefully
- review QC output before trusting downstream results
- compare results across multiple thresholds when appropriate
- validate prioritized hits experimentally
- treat simulation and virtual design outputs as hypothesis-generating, not final truth

## Evaluation

Current evaluation is based on:
- workflow functionality
- internal consistency of computed metrics
- usability for exploratory chemoproteomic analysis
- deployment and interactive use in a browser

Future evaluation should include:
- benchmarking against real chemoproteomic datasets
- calibration of virtual library ranking against follow-up screening outcomes
- robustness testing across assay formats and sparse-label settings

## Ethical and Safety Considerations

CysSelect should not be used in high-stakes clinical or regulatory settings. It is a research tool intended to support scientific reasoning, not replace it.

## Maintenance

Recommended maintenance tasks:
- dependency updates
- threshold and scoring review
- simulation model refinement
- expanded annotation support
- stronger model validation with public datasets

## Citation

Petalcorin, M.I.R. (2026). CysSelect: a browser-based platform for chemoproteomic hit calling, enantiomer comparison, simulation, and next-library design for covalent fragment discovery. https://github.com/mpetalcorin/cysselect-app