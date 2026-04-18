# Datasheet: CysSelect

## Dataset / System Name

CysSelect Input and Output Data Framework

## Purpose

This datasheet documents the kinds of data handled by the CysSelect application. The goal is to clarify what information the app expects, what it generates, and how these data should be interpreted in a research setting.

## Motivation

CysSelect was built to provide a practical interface for chemoproteomic screening analysis, especially for workflows involving covalent fragment enantiopairs, cysteine engagement, pair-level prioritization, and virtual library design.

The system supports both:
- real uploaded screening tables
- internally generated simulation datasets

## Data Handled by the App

### 1. Real Uploaded Screening Tables
These are user-provided tabular datasets representing chemoproteomic screening measurements.

### 2. Simulated Benchmark Datasets
These are synthetic tables generated inside the app for benchmarking, teaching, and workflow testing.

### 3. Derived Analytical Tables
These include processed hit tables, pair-level summaries, and ranked outputs derived by the app.

### 4. Virtual Library Candidate Tables
These are user-uploaded or app-generated candidate compound tables used for ranking and next-batch selection.

## Required Input Fields

For uploaded screening tables, the following columns are required:

- `pair_id`
- `enantiomer`
- `cysteine_id`
- `concentration_uM`
- `replicate`
- `intensity_dmso`
- `intensity_condition`

These fields define the minimum information needed to calculate core screening metrics.

## Optional Input Fields

Optional metadata columns include:

- `p_value`
- `protein_id`
- `gene_name`
- `scaffold_cluster`

These fields improve interpretation and filtering but are not required for core calculations.

## Data Generation and Transformation

### Real Data
Uploaded raw tables are transformed into:
- competition ratios
- percent competition
- pair‚Äìcysteine summaries
- enantiomer comparison outputs
- pair-level summary tables
- priority-ranked tables

### Simulated Data
Simulation mode generates synthetic screening tables using configurable assumptions such as:
- number of enantiopairs
- number of cysteines
- assay noise
- missing-data rate
- latent interaction behavior

### Virtual Library Data
Virtual library mode generates or accepts candidate tables and produces:
- predicted promiscuity
- predicted stereoselectivity fraction
- soft sweet-spot scores
- acquisition scores
- selected next-round batches

## Data Quality Considerations

Input data quality strongly influences output quality.

Important factors include:
- correct column names
- presence of both enantiomers
- valid concentration values
- nonzero and nonmissing intensities
- reasonable replicate structure
- trustworthy optional metadata

The app includes QC checks to identify obvious structural problems before analysis.

## Recommended Use Conditions

The data framework is suitable for:
- exploratory chemoproteomic analysis
- prototyping and method development
- educational demonstration
- workflow benchmarking
- early-stage hit triage

## Not Recommended For

The framework is not intended for:
- clinical records
- patient identifiers
- protected health information
- confidential human subject data without proper governance
- final therapeutic decision-making without validation

## Privacy and Sensitivity

Users are responsible for the data they upload. Public deployments of the app should avoid using sensitive or confidential datasets unless appropriate access controls are in place.

Demo and simulation data should be:
- synthetic
- de-identified
- nonconfidential

## Potential Sources of Bias or Error

Potential issues include:
- noisy assay measurements
- incomplete replicates
- misassigned enantiomer labels
- threshold sensitivity
- simplified simulation assumptions
- overreliance on pair-level heuristics
- sparse positives in prioritization workflows

## Data Outputs

The app can export:
- uploaded raw tables
- processed hit tables
- pair-level summaries
- simulated benchmark tables
- scored virtual libraries
- selected candidate batches

Users should inspect exported tables before downstream use.

## Update and Versioning Guidance

When the data schema changes, update:
- the README
- the model card
- this datasheet
- any deployment instructions
- demo files and templates

## Recommended Future Improvements

- support for richer annotation schemas
- stronger schema validation
- versioned demo datasets
- clearer provenance tracking for simulated outputs
- optional integration with external cysteine knowledge resources

## Summary

CysSelect handles structured chemoproteomic screening inputs and converts them into interpretable analytical outputs for hit calling, enantiomer comparison, simulation, prioritization, and virtual library design. The system is best suited for research and educational workflows and should be used with appropriate attention to input quality, threshold choice, and experimental validation.
