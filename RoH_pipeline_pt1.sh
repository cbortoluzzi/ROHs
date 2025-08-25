#!/bin/bash


# Author : Chiara Bortoluzzi


#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=Heterozygosity
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 3 ]
then
	echo -e "\nusage: `basename $0` <VCF> <bam> <cov>\n"
	echo -e "DESCRIPTION: Calculate a corrected measure of heterozygosity in a 10-Kb window\n\n"
	echo -e "INPUT:       <VCF>         The filtered VCF file"
 	echo -e "			  <bam>			The alignment file in BAM format"
	echo -e "             <coverage>    The coverage file obtained from samtools depth\n"

	echo -e "OUTPUT:      <corrected_heterozygosity>   A tab-delimited file with the corrected measure of heterozygosity in a 10-Kb window\n"
	exit
fi


# EXAMPLE: run python script using default options for minimum depth and window size
python3 calculate_genome_wide_heterozygosity.py --vcf $vcf --bam $bam --cov $cov --o genome_wide_heterpzygosity

# EXAMPLE: run python script changing the minimum depth parameter while maintaining the default value for the window size
python3 calculate_genome_wide_heterozygosity.py --vcf $vcf --bam $bam --cov $cov --dp 10 --o genome_wide_heterpzygosity


