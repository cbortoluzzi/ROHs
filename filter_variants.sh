#!/bin/bash


# Author : Chiara Bortoluzzi


#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=Filter_VCF
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 2 ]
then
	echo -e "\nusage: `basename $0` <VCF> <coverage>\n"
	echo -e "DESCRIPTION: Filter VCF file obtained from running DeepVariant\n\n"
	echo -e "INPUT:       <VCF>         The VCF file obtained from DeepVariant"
	echo -e "             <coverage>    The coverage file obtained from samtools depth\n"

	echo -e "OUTPUT:      <filtered_vcf>   A filtered VCF file\n"
	exit
fi


vcf=$1 
cov=$2


# EXAMPLE: if running the python script using the default parameters for filtering the VCF file
python3 filter_variants.py --vcf $vcf --cov $cov --o filtered_vcf

# EXAMPLE: if running the python script using different parameters from the default ones
python3 filter_variants.py --vcf $vcf --cov $cov --q 20 --dp 8 --gq 30 --o filtered_vcf

