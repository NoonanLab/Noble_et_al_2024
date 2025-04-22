#!/bin/bash 

# Paths to reference genome and minimal VCF with p53 mutation
FASTA="/gpfs/gibbs/pi/noonan/man59/References/CellRanger/Hg38_dCas9_CellRanger/fasta/genome.fa"
VCF="/home/man59/gibbs/Valid_p53/p53_mutation.vcf"

# Number of threads per VarTrix job
THREADS=4

# Loop over each lane directory
for lane_dir in Rep*_L*; do 
	(
	
	echo "Running VarTrix on $lane_dir"
	
	BAM="${lane_dir}/outs/possorted_genome_bam.bam"
	BARCODE="${lane_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
	OUTNAME="vartrix_output_${lane_dir}"
	
	/gpfs/gibbs/pi/noonan/man59/Valid_p53/./vartrix_linux \
	 --bam="$BAM" \
	 --cell-barcodes="$BARCODE" \
	 --fasta="$FASTA" \
	 --vcf="$VCF" \
	 --out-matrix="${OUTNAME}_alt" \
	 --ref-matrix="${OUTNAME}_ref" \
	 --scoring-method=coverage \
	 --threads="$THREADS"
	
	echo "Done with $lane_dir"
	) & 
done 

# Wait for all parallel jobs to finish
wait
echo "All VarTrix jobs completed"
