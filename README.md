# internal-tide-eddy-cascade
Code for analyzing internal tide–eddy interactions and forward energy cascades in oceanic flows.
Code and analysis for the paper "Internal tide-eddy interactions enhance the forward energy cascade of oceanic flows".

This repository contains the model configuration files, parallel structure-function calculation codes, data preprocessing scripts, analysis functions, and figure-generation scripts used in the study. The main workflow is designed to quantify how internal tide-eddy interactions modify flow structures and enhance the forward kinetic energy cascade in oceanic flows.

## Repository structure

- `model config/`
  - Model configuration files used for the numerical experiments.

- `SF calc python/`
  - Parallel Python codes for calculating structure functions on high-performance computing platforms.

- `SF analysis/`
  - Codes for data preprocessing, post-processing, diagnostic analysis, and figure generation.
  - `preprocessing/`
    - Interpolation from model S-coordinates to z-coordinates.
    - Extraction of the target eddy-region data.
    - Low-pass filtering.
    - Single-layer, 3-layer, and 15-layer preprocessing scripts.
  - `function/`
    - Common functions used in the analysis.
    - Functions for spectral flux, cospectrum, bootstrapping, and other diagnostics.
  - `calc analysis/`
    - Post-processing of structure-function results calculated by `SF calc python/`.
    - Errorbar estimation.
    - Cospectrum calculation.
  - `print figures/`
    - Plotting scripts used to generate the figures in the paper.

## Directory description

### `model config/`

This directory contains the configuration files for the numerical model experiments used in the study.

### `SF calc python/`

This directory contains the parallel Python codes developed for calculating velocity and vorticity structure functions on supercomputing platforms.

The structure-function calculation requires pairwise computations between grid points. Its computational complexity is approximately `O(N^2)`, where `N` is the total number of grid points. For large model datasets, direct serial computation is prohibitively expensive. Therefore, these scripts were developed for parallel execution on high-performance computing platforms.

### `SF analysis/`

This directory contains the main data analysis and figure-generation workflow.

#### `SF analysis/preprocessing/`

This folder contains scripts for preprocessing model outputs before structure-function and spectral analyses, including:

- interpolation from model S-coordinates to z-coordinates;
- extraction of the target eddy region;
- low-pass filtering;
- preparation of single-layer, 3-layer, and 15-layer datasets.

The single-layer, 3-layer, and 15-layer scripts are separated because they are used for different scientific purposes:

- single-layer data: horizontal structure-function and energy-cascade analyses;
- 3-layer data: vertical shear-related diagnostics;
- 15-layer data: vertical-structure analyses.

#### `SF analysis/function/`

This folder contains commonly used functions for data analysis, including routines for:

- spectral flux calculation;
- cospectrum calculation;
- bootstrapping-based uncertainty estimation;
- other helper functions used throughout the analysis.

#### `SF analysis/calc analysis/`

This folder contains scripts used after the structure functions are computed by `SF calc python/`, including:

- post-processing of structure-function results;
- error-bar calculation;
- cospectrum calculation.

#### `SF analysis/print figures/`

This folder contains the plotting scripts used to generate the figures in the paper.

## Workflow overview

The general workflow of this repository is:

1. Prepare model outputs and configuration files.
2. Preprocess the model data, including coordinate interpolation, target-region extraction, and filtering.
3. Calculate structure functions using the parallel codes in `SF calc python/`.
4. Post-process the structure-function results and estimate uncertainties.
5. Calculate spectral fluxes, cospectra, and related diagnostics.
6. Generate the figures used in the paper.

## Notes

The codes in this repository were developed for the analyses presented in the accompanying manuscript. Some scripts may require users to modify file paths, input file names, environment settings, or job-submission settings before running them on a different machine or cluster.
