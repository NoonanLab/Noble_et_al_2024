#!/usr/bin/python

'''
Usage: python 03_flashfry_filter.py infile outfile
Filters FlashFry outputs based on all outputs for CRISPRi. No off-target shall have a higher specificity score than the target.

Examples:
python 03_flashfry_filter.py HARS3_linars_panTro6_300bp.output.scored hars_panTro6.txt
python 03_flashfry_filter.py HARS3_linars_hg38_300bp.output.scored hars_hg38.txt
'''

# Import command line arguments
# Define f_in as input file in read mode
# Define f_out as output file write mode
import sys
f_in = open(sys.argv [1], 'r')
f_out = open(sys.argv [2], 'w')

# Filter for loop
for line in f_in: 
	
	# For each line, separate by tab
	split_line = line.strip().split('\t')
	
	# If the first column is the header (i.e says "contig"), then print the line
	if split_line[0] == 'contig':
		f_out.write(line)
	
	else:	
		# Get fields
		hsu_value = float(split_line[13])
		jost_ot = float(split_line[17])
		jost_spec = float(split_line[18])
		gc_flag = split_line[10]
		polyt_flag = split_line[11]
	
		# Get mismatch column, separate by comma
		split_ot = (split_line[16].strip().split(',')) 
		mismatch_0 = split_ot[0]
		mismatch_1 = split_ot[1]
		mismatch_2 = split_ot[2]
	
		# Print line if:
		# Hsu value is greater than 65, off-target score less than or equal to 0.4
		# Specifity score greater than or equal to 0.5, no PolyT or GC flags (NONE)
		# There is only 1 perfect match, and no 1 or 2 bp mismatches
		if (hsu_value >= 65.0 and jost_ot <= 0.4 and jost_spec >= 0.5 and
		    gc_flag == 'NONE' and polyt_flag == 'NONE' and
		    mismatch_0 == '1' and mismatch_1 == '0' and mismatch_2 == '0'):
	        	f_out.write(line)

# Close both files
f_in.close()
f_out.close()
