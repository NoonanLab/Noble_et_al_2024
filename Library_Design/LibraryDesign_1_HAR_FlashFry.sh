#Use expanded HARs as regions for sgRNA discovery using FlashFry


#Expand HAR by 300bp each side via BEDTools

bedtools slop -i hars_hg38.bed -g hg38.chrom.sizes -b 300 > hars_hg38_300bp.bed
bedtools slop -i hars_panTro6.bed -g panTro6.chrom.sizes  -b 300 > hars_panTro6_300bp.bed

#Convert to FASTA Format

bedtools getfasta -fi hg38.fa -bed hars_hg38_bp300.bed > hars_hg38_bp300.fa
bedtools getfasta -fi panTro6.fa -bed hars_panTro6_300bp.bed > hars_panTro6_300bp.fa

#Discover guides with FlashFry
#hg38
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
 discover \
--database hg38_cas9ngg_database \
--fasta hars_hg38_300bp.fa \
--output hars_hg38_300bp.output
 
#panTro6
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
discover \
--database panTro6_spcas9ngg_database \
--fasta hars_panTro6_300bp.fa \
--output hars_panTro6_300bp.output 
 
#Score sites
 
#hg38
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
score \
--input hars_hg38_300bp.output \
--output hars_hg38_300bp.output.scored \
--scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,JostandSantos \
--database hg38_cas9ngg_database
 
#panTro6
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
score \
--input hars_panTro6_300bp.output \
--output hars_panTro6_300bp.output.scored \
--scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,JostandSantos \
--database panTro6_spcas9ngg_database
