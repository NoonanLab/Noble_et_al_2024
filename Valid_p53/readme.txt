README

We identified a background rate of TP53 mutant reads in the perturb-seq. 
In this directory:
Running Vartrix on the Perturb-seq data to identify TP53 mutant cells. 
Assign these calls to cells in the Seurat object. 
Visualize mutant and WT cells. 
Perform DE on mutant and WT cells. 
Re-analyze Batch 3 (unaffected) using linear regression for HAR effects on module expression. 
Results of these analyses. 

STEPS

STEP 1. RUN VARTRIX 

/Scripts
-vartrix_run.sh 
Batch runs VarTrix to identify wt and p53 mutant reads from the bam files for all lines across the perturb experiment. 
-Test_Vartrix_PostCall.R
Script tests one example (one lane) of the data for debugging.
-Analyze_VarTrix.R 
Script combines alt and ref matrices per lane and assigns p53 status to barcodes; updates barcodes so that lanes are 1-24 rather than 3 sets of 1-8. Summarizes by-batch mutation rates. 

/Batch_Vartrix_Outputs
-Contains outputs of Vartrix from all lanes. Intermediate files.
-Ref matricies contain a count matrix for cells with the wt p53 allele. 
-Alt matrices contain a count matrix for cells with the mutant p53 allele. 
-These must be paired together, and the CBCs must be modified to match the structure in the Seurat object, using Analyze_Vartrix.R 

/Outs
-Vartrix_CBC_Results.txt
Calls for every cell barcode of the counts for the p53 mutation and calls for wt versus p53.mutant. 

-Vartrix_Batch_Summary.txt
Summary of mutation rate per batch. 

-Vartrix_Filtered_CBC_p53.txt
Same as Vartrix_CBC_Results, but only for cells in the final filtered count matrix (singles). 

-Vartrix_Batch_Summary_FilteredCells.txt
Same as Vartrix_batch_Summary, but only for cells in the final filtered count matrix (singles). 

STEP 2: ANALYZE RESULTS

/Scripts
-Analyze_p53_Seurat.R
Uses cell classifications to visualize p53 mutant cells in UMAP. 
Runs differential expression using DESeq2 for P53 mutant versus wild-type cells. 

-Analyze_Batch3.R
Re-Runs linear regression analysis on unaffacted Batch 3. 
Compares these results to the full Perturb LR results using Pearson Correlation. 


OUTS 
-DimPlot_p53_vs_WT.pdf
Plot highlighting p53 mutant and wild-type cells from Analyze_p53_Seurat.R. 

-VolcanoPlot_p53_v_WT.pdf 
Volcano plot showing results of differential expression analysis from Analyze_p53_Seurat.R 

-Module_LR_Batch3_Reanalysis.txt
Results of Batch 3 linear regression re-analysis. 
