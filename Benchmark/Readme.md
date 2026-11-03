The source data of RNA 3D structure predicitons for 

## HuRNA: A Hybrid framework for universal RNA 3D Structure Prediction

This repository contains source data (`3D_Structure_predicitons.zip`) for the HuRNA benchmarks, including prediction inputs, reference structures, HuRNA predictions, and, for some datasets, predictions from other methods.

## Contents in `3D_Structure_predicitons.zip`:

- [`CASP15/`](CASP15/) and [`CASP16/`](CASP16/): CASP RNA benchmark datasets.
- [`RNA-Puzzles/`](RNA-Puzzles/): RNA-Puzzles benchmark data.
- [`DRfold2_test_set/`](DRfold2_test_set/): DRfold2 test set, with native 2D and sequence-only inputs.
- [`RNA_Dimer/`](RNA_Dimer/): RNA dimer cases.
- [`alternative_structure/`](alternative_structure/): Alternative-structure cases.
- [`circular_RNA/`](circular_RNA/), [`designed_RNA/`](designed_RNA/), and [`synthetic_RNAs/`](synthetic_RNAs/): Circular, designed, and synthetic RNA cases.
- [`end-to-end/`](end-to-end/): End-to-end predictions of selected targets from RNA Puzzles.
- [`prediction_with_restriants/`](prediction_with_restriants/): Predictions using restraints.
- [`chain_sel_configuration/`](chain_sel_configuration/): JSON configurations selecting reference and predicted PDB chains and residue ranges for metric calculation.

Dataset folders generally organize prediction inputs and reference structures under `input/` and predicted structures under `predictions/`. Predictions from other methods, where available, are grouped by method.

********************************************************
Citation:
Wang et al. HuRNA: A Hybrid framework for universal RNA 3D Structure Prediction

Please contact us at zhangd_iqb@zju.edu.cn for more information.
