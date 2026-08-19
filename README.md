# Population-Based-IVT-Model-Example-MATLAB-Implementation

This repository provides a representative MATLAB implementation of the population-based mechanistic in vitro transcription (IVT) process model described in:
**Mendis, P. and Kis, Z. "A Population-Based Mechanistic Process Model of In Vitro Transcription with Explicit Tracking of Intermediate Transcription Complexes"**

The example corresponds to Case Study II of the manuscript, based on experimental data for a 12-mer RNA construct reported by Young et al. (1997).

## Description

The code simulates the concentration-time profiles of runoff RNA and total abortive products using the optimized model parameters reported in the manuscript. Three experimental conditions with different NTP and Mg concentrations are simulated and compared with the corresponding experimental data.

This implementation is provided as an example for transparency and model simulation. It does not reproduce the parameter-estimation procedure used to obtain the optimized parameter values.

## Requirements

- MATLAB
- `gridded_lookup_table_5_50_50_50.mat`

The lookup table is used for calculation of the ionic equilibria during model simulation.

## Running the example

Place the MATLAB script and lookup-table file in the same directory and run the main script.

The script:
1. loads the ionic-equilibrium lookup table;
2. defines the experimental conditions and optimized parameters;
3. solves the population-model equations using `ode15s`; and
4. generates a figure comparing the model predictions with the experimental data.

## Reference
Young, J. S.; Ramirez, W. F.; Davis, R. H. Modeling and Optimization of a Batch Process for In Vitro RNA Production. Biotechnol. Bioeng. 1997, 56 (2), 210–220.

For the complete mathematical formulation, assumptions, parameter definitions, and parameter-estimation methodology, please refer to the manuscript and Supporting Information.
