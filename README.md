# ROHs
The repository contains the scripts required to identify runs of homozygosity (ROHs) using whole-genome sequencing data. The scripts implement the method developed by Bosse et al. (2012). 

# How does it work?
A region of homozygosity (ROHS) is a contiguous genomic stretch of homozygous genotypes characterized by less variation in an individual than is expected based on the genomic average. These genomic regions are inherited from a common ancestor by both parents (identical haplotypes), and therefore indicates a certain level of relatedness. ROHs are here identified using information on heterozygosity calculated in non-overlapping consecutive windows along the genome. 

# How do I run it?
Before identifying ROHs, it is necessary to create the input file (step 1) for calculating heterozygosity (step 2) in non-overlapping windows.
Step 1: make_json.py
The script takes three input files:
1. A file containing the path and name of each individual bam file, one per line (-b option)
2. A file containing the path and name of each individual vcf file, one per line (-v option)
3. A Tab delimited file with information on genome-wide coverage for each individual (-c option)
In order to calculate heterozygosity on high confident variants, it is also necessary to set a minimum depth value to consider a SNP (-m option). In case of whole-genome sequencing data with at least 10x homogenous coverage along the genome, we recommend to set this value to at least 4x. The script will then consider only SNPs with a depth of coverage between 4x and 2*$average genome-wide coverage. 

Step 2: calculate_heterozygosity.py
The script calculates the level of heterozygosity in non-overlapping consecutive windows along the genome for each individual. To do so, the script takes three input parameters:
1. A indexed fasta file with information on chromosome and total length (-f option) (i.e. this indexed fasta file can be generated with samtools faidx using the reference genome in fasta format of your species of interest)
2. The json file of your individual for which you want to calculate binned heterozygosity (-j option)
3. The window size in base pairs (-b option)



Nucleotide diversity was calculated for bins of 10 kbp over the entire genome within each individual. “SNPbin” is the SNP count per 10 kbin, corrected for the number of bases within that bin that was not covered enough for the VCFtools filtering, so that the eventual SNP count per bin (SNPbin) is proportional to 10.000 covered bases. SNPcount = total number of SNPs counted in a bin of 10 kbp. DP = coverage in bp/bin (per base at least depth of 7× and maximum of ∼2*average coverage). Binsize = 10.000. Correction factor = DP/binsize. SNPbin = SNPscount/Correction factor.
