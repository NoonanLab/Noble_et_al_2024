# **Custom code and analysis for Noble et al (2024)**

This repository contains directories for the sgRNA library design pipeline and subsequent Perturb-seq analysis. 

## Library Design
Pipeline to design high-specificity sgRNA sequences in two species, human and chimpanzee, targeting Human Accelerated Regions of interest. 

### Contents

**LibraryDesign_1_HAR_FlashFry.sh** – Commands to apply FlashFry (McKena et al., 2018) to set of HARs and chimpanzee orthologs to discover all possible sgRNAs. 

Usage: Import list of HAR genomic positions in Hg38 and PanTro6. 

**LibraryDesign_2_FlashFry_Filter.py** - Filter potential sgRNAs generated for specificity metrics provided by FlashFry. 
```
$ python 03_flashfry_filter.py infile outfile
```

Usage: Apply to output of Step 1. 

**LibraryDesign_3_SpeciesFilter.R** - Filter potential sgRNAs passing the Step 2 filter across species, to reduce the possibility for species-specific genomic background effects. 

Usage: Apply to output of Step 2. 
```
$ R CMD BATCH LibraryDesign_3_SpeciesFilter.R
```


**LibraryDesign_4_SelectGuides.R** - Select the top-scoring sgRNA for each region of interest and the next-highest scoring sgRNA whose center is at least 10bp from the first guide, to generate a library of 2 sgRNAs for each region of interest. 

Usage: Apply to output of Step 3. 
```
$ R CMD BATCH LibraryDesign_4_SelectGuides.R
```

## Perturb-seq analysis 
Analysis pipeline to used to interpret transcriptional effects resulting of HAR perturbation.

### Contents 
**GuideAssingment.R** – Identify the transduced sgRNA in each cell. This pipeline accounts for the distribution of sgRNA counts throughout the cell population, as well as the distribution of counts within cells. 

Usage: Import Perturb-seq CRISPR-capture count matrix and Feature Reference indicating guide sequences for each target. 
```
$ R CMD BATCH GuideAssignment.R
```

**Perturb_DifferentialExpression.R** – Perform Wilxocon Rank Sum differential expression for each perturbation in the Perturb-seq experiment. This compares all cells with bearing sgRNAs against a particular target against all cells bearing NTC sgRNAs. 

Usage: Import Perturb-seq filtered gene expression matrix. 
```
$ R CMD BATCH Perturb_DifferentialExpression.R
```

**LinearRegression_byGuide.R** – Perform Multiple Linear Regression to estimate effect sizes on hdWGCNA module expression, with the ‘perturbation’ variable levels being each individual sgRNA. 

Usage: Import Perturb-seq hdWGCNA module expression matrix. 
```
$ R CMD BATCH LinearRegression_byGuide.R
```

**LinearRegression_byPerturbation.R** – Perform Multiple Linear Regression to estimate effect sizes on hdWGCNA module expression, with the ‘perturbation’ variable levels being each target in the screen (i.e., both sgRNAs for each HAR are analyzed together).

Usage: Import Perturb-seq hdWGCNA module expression matrix. 
```
$ R CMD BATCH LinearRegression_byPerturbation.R
```
![image](https://github.com/NoonanLab/Noble_et_al/assets/100241154/f758f98c-1f3e-40c1-8b03-388e1dde221d)
