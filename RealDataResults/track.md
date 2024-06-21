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
3. In Chr 8 for structrual MRI, 5 genes tested positive through regional and tissue volume IDPS, 5 genes tested positive through cortical area, 5 genes tested positive through  cortical thickness, 1 gene tested positive through cortical grey-white contrast, 1 gene tested positive through regional and tissue intensity
4. In Chr 8 for diffusion MRI, 9 genes tested positive through WM tract FA, 9 genes tested positive through WM tract diffusivity, 7 genes tested positive through WM tract ICVF, 9 genes tested positive through WM tract OD
5. In Chr 19 for structrual MRI, 40 genes tested positive through regional and tissue volume IDPS, 40 genes tested positive through cortical area, 40 genes tested positive through  cortical thickness, 3 gene tested positive through cortical grey-white contrast, 40 gene tested positive through regional and tissue intensity
6. In Chr 19 for diffusion MRI, 38 genes tested positive through WM tract FA, 38 genes tested positive through WM tract diffusivity, 38 genes tested positive through WM tract ICVF, 38 genes tested positive through WM tract OD, 38 genes tested positive through WM tract ISOVF

# No adjustment for other MRI
1. UKB diffusion: found 3 significant genes in Chr 22, but no gene significant in Chr 8, 20 genes in Chr 19, all in adjusted
2. IGAP diffusion: found 8 significant genes in Chr 1, 9 significant genes in Chr 2, 41 significant genes in Chr 19 (ZNF227, ZNF233, ZNF235 were not significant with adjustment), but no gene significant in Chr 8; IGAP structural: 9 significant genes in Chr 2, 39 significant genes in Chr 19 (DKFZp434J0226, IGFL2 in adjustment but not here, PPP1R13L here is not significant in adjustment), but no gene significant in Chr 8

# Overlapping genes
1. Chr 19 PVR and MIR4531; BCL3 and MIR8085; APOC4-APOC2; APOC2 and APOC4 overlaps; KLC3 overlaps ERCC2
2. Chr 8, none

# To discuss
1. MIR4287 has negative pvalue for sMRI when using Davies, but fixed when switched to Liu
2. IDP region
3. MV-IWAS but used our method, do we want to use MVIWAS summary statistics version
