# MV-Modality-IWAS Analysis Pipeline

## Data Acquisition

The data acquisition process is described in detail in [DataCleaning.md](https://github.com/junjypark/MV_VC_IWAS/blob/main/DataCleaning/DataCleaning.md), with all related files located in the DataCleaning folder. Data cleaning consists of three main steps:

1. **Data Collection**: Obtain the UKB, 1000G, and IGAP datasets.
2. **Data Reprocessing**: Intersect the three datasets to ensure that only common SNPs are included.
3. **Data Quality Control**: Retain only SNPs with a MAF greater than 0.01 and perform clumping using an \( R^2 \) threshold of 0.5.

## Methods

The proposed methods are implemented in the R folder. The main function is `mv_vc_iwas`, which can be used to perform the MV-Modality-IWAS analysis. The function takes the following parameters:
- `corr`: p times 1 correction vector (from AD GWAS)
- `LD`: p times p LD matrix of the variants
- `A1`: p times q1 matrix of coefficients predicting imaging data (in modality of testing interest) from genotypes (from IDP GWAS)
- `A2`: p times q2 matrix of coefficients predicting imaging data from genotypes (from IDP GWAS)
- `method`: Either "davies" or "Liu". Davies method is used as a default.


## Simulations

Simulation studies are available in the Simulation folder. These include simulations for both one-sample and two-sample scenarios, covering binary and continuous outcomes.

## Real Data Analysis

Real data analyses scripts are included in the DataAnalysis folder, with separate analyses for IGAP and UKB. 

- Step 1: The first step of the analysis is performed using the `run_IGAP_part1.R` and `run_UKB_part1.R` scripts. These scripts perform the preprocessing of the data, including the calculation of the correlation vector and LD matrix. The results are to be used in Step 2

- Step 2: The second step of the analysis is performed using the `run_IGAP_part2_S.R`, `run_IGAP_part2_D.R`, `run_IGAP_part2_F.R` and `run_UKB_part2_S.R`, `run_UKB_part2_D.R`, `run_UKB_part2_F.R` scripts. Each script corresponds to a specific dataset (IGAP or UKB) and a particular MRI modality of testing interest (structural, diffusion, or functional). These scripts perform the MV-Modality-IWAS analysis using the results from step 1, together with IDP GWAS data. 

Additionally, the MV-IWAS analysis used for comparison with the proposed method can be found in the MVIWAS folder.
