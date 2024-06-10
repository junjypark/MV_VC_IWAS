# Work track
0. Data location
   * Results located in: `/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS`
   * IDP id, description and category data located in: `/gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/IDP_category/IDP_description.csv`
     structural MRI ID: 1-1451, diffusion MRI ID: 1452-2126
   * Stage 2 final pvalue read in from: `/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/IGAP_results`
   * Stage 1 IDPs being included in Stage 2 analysis read in: `/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/`

1. `IDP_results.R`
   Obtain MRI IDP ID, category, MRI significance results


# Current questions
  **Not sure if the category I created for IDP is too general**
  
  **dMRI_reprocess.csv seems not correct**: max diffusion MRI ID should be 2100 according to the current dMRI_reprocess.csv and sMRI_reprocess.csv file. But for Chr1, gene 1, diffusion MRI IDP include ID 2110

   **Need to check if gene index match**: Also have some concern that if gene number in Stage 1 IDPs being included in Stage 2 analysis read in: `/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/` match the gene names in Stage 2 final pvalue read in from: `/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/IGAP_results`
