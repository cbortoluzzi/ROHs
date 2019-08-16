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
In order to calculate heterozygosity on high confident variants, it is also necessary to set a minimum depth value to consider a SNP (-m option). In case of whole-genome sequencing data with at least 10x homogenous coverage along the genome, we recommend to set this value to at least 4x. The script will then consider only SNPs with a depth of coverage between 4x and 2*average genome-wide coverage. 

Step 2: calculate_heterozygosity.py
The script calculates the level of heterozygosity in non-overlapping consecutive windows along the genome for each individual. To do so, the script takes three input parameters:
1. A indexed fasta file with information on chromosome and total length (-f option) (i.e. this indexed fasta file can be generated with samtools faidx using the reference genome in fasta format of your species of interest)
2. The json file of your individual for which you want to calculate binned heterozygosity (-j option)
3. The window size in base pairs (-b option)
We recommend to use a window size of 10 kb (i.e. -b 10000) for the calculation of the heterozygosity.


Step 3: identify_rohs.py
The script uses the output file of step 2 (-i option) on binned heterozygosity to identify ROHs. To identify ROHs a set of parameters are required:
1. The window size in base pairs used in step 2 (-b option)
2. The minimum number of well-covered sites identified in the window (-t1 option)
3. The threshold to filter out SNPs within a candidate autozygous stretch (-t2 option)
4. The threshold for the average heterozygosity within a candidate autozygous stretch (-t3 option)
5. Minimum number of consecutive windows (-t4 option)
6. Threshold for local maximum (-t5 option)
For more information on the parameters required to identify ROHs, refer to the paper of Bosse et al. (2012) and Bortoluzzi et al. (2019). 


# How do I interpret the output?
The most important files to look at are Sample_input.txt and Sample_ROHs.txt
The first file (Sample_input.txt) is the output file generated by step 2. It contains the following information:
1. Sample ID
2. Chromosome 
3. Start position of the window
4. End position of the window
5. Total number of well-covered sites  (4x-2*average coverage) calculated from the BAM file
6. Total number of well-covered heterozygous sites (4x-2*average coverage) calculated from the VCF file. 

The second file (Sample_ROHs.txt) is the output file generated by step 3 and contains information on all identified ROHs. The informations are:
1. Sample ID
2. Chromosome
3. Start position of ROHs
4. End position of ROHs
5. ROH size (how many windows of e.g. 10 kb are included in the ROHs)
6. ROH length (how many windows of e.g. 10 kb are included in the ROHs, including those that did not meet the right coverage criteria)
7. SNP count 
Ideally, the ROH size and ROH length are the same. However, in some cases the length is much higher than size, meaning that many uncovered windows make up the ROH and this ROH is therefore of low quality. If this is the case, it is recommended to filter out these ROHs (e.g. if size/length < 0.5, then exclude the ROH). 
The SNP count is the number of heterozygous sites called in a window corrected for the number of sites within that window that were not covered enough. In other words, the corrected number of heterozygous SNPs is calculated as:
SNPcount = (widnow size/number of well-covered sites in a window) * number of well-covered heterozygous sites in a window

# How do I cite it?
If you use this pipeline in your research, please cite the following two papers: 
Bosse et al. "Regions of homozygosity in the porcine genome: consequence of demography and the recombination landscape".  Plos Genetics (2012).
Bortoluzzi et al. "The type of bottleneck matters: insights into the deleterious variation landscape of small managed popualtions". Evolutionary Applications (2019).

# Requirements
Python +3.6

