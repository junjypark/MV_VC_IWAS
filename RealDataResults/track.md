# Work track
0. Data location
   * Results located in: `/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS`
   * IDP id, description and category data located in: `/gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/IDP_category/IDP_description.csv`
     structural MRI ID: 1-1436, diffusion MRI ID: 1437-2126
   * Stage 2 final pvalue read in from: `/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/IGAP_results`
   * Stage 1 IDPs being included in Stage 2 analysis read in: `/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/`

1. `IDP_results.R`
   Obtain MRI IDP ID, category, MRI significance results, results located in `/gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/IDP_category` by chr

# Significant IDP categories
1. Structural MRI are categorized to regional and tissue volume (647), cortical area (372), cortical thickness (306), cortical grey-white contrast (70), regional and tissue intensity (41).
2. Difussion MRI are categorized to white matter hyperintensity volume (1), regional T2* (14), WM tract FA (75), WM tract MO (75), WM tract diffusivity (300), WM tract ICVF (75), WM tract OD (75), WM tract ISOVF (75)
3. In Chr 8, 5 gene tested positive through regional and tissue volume IDPS, 5 gene tested positive through cortical area, 5 gene tested positive through  cortical thickness, 1 gene tested positive through cortical grey-white contrast, 1 gene tested positive through regional and tissue intensity

# Current questions
  **Not sure if the category I created for IDP is too general**
  
  **dMRI_reprocess.csv seems not correct**: max diffusion MRI ID should be 2100 according to the current dMRI_reprocess.csv and sMRI_reprocess.csv file. But for Chr1, gene 1, diffusion MRI IDP include ID 2110

   **Need to check if gene index match**: Also have some concern that if gene number in Stage 1 IDPs being included in Stage 2 analysis read in: `/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/` match the gene names in Stage 2 final pvalue read in from: `/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/IGAP_results`
