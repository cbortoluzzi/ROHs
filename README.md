# ROHs
This repository contains the scripts required to identify runs of homozygosity (ROHs) using whole-genome sequencing data. The scripts implement the method developed by Bosse et al. (2012). 

# How does it work?
Runs of homozygosity (ROHs) are identified using information on heterozygosity calculated in non-overlapping consecutive windows across the genome. ROHs are here defined as genomic regions showing lower heterozygosity than expected based on the genome-wide heterozygosity. 

# How do I run it?
Before identifying ROHs, it is necessary to create the input file (step 1) and calculate heterozygosity (step 2) in non-overlapping windows.
Step 1: make_json.py
This script takes three input files:
1. A file containing the path and name of each individual bam file, one per line (-b option)
2. A file containing the path and name of each individual vcf file, one per line (-v option)
3. A Tab delimited file with information on genome-wide coverage for each individual (-c option)
The script will output a json file, one for each individual, with all the information required for step 2. 


