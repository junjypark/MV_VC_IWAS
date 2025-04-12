# MV-Modality-IWAS Analysis Pipeline

## Data Acquisition

The data acquisition process is described in detail in [DataCleaning.md](DataCleaning.md), with all related files located in the DataCleaning folder. It consists of three main steps:

1. **Data Collection**: Obtain the UKB, 1000G, and IGAP datasets.
2. **Data Reprocessing**: Intersect the three datasets to ensure that only common SNPs are included.
3. **Data Quality Control**: Retain only SNPs with a MAF greater than 0.01 and perform clumping using an \( R^2 \) threshold of 0.5.

## Methods

The analysis methods are implemented in the R folder.

## Simulations

Simulation studies are available in the Simulation folder. These include simulations for both one-sample and two-sample scenarios, covering binary and continuous outcomes.

# Real Data Analysis

Real data analyses are performed in the DataAnalysis folder, with separate analyses for IGAP and UKB. Additionally, the MV-IWAS analysis used for comparison with the proposed method can be found in the MVIWAS folder.