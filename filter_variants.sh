#!/bin/bash


# Author : Chiara Bortoluzzi


#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=Filter_VCF
#SBATCH --output=output_%J
#SBATCH --error=error_%J



# Input files
vcf=$1 # VCF file obtained from DeepVariant
cov=$2 # Coverage file obtained from samtools depth


# EXAMPLE: if running the python script using the default parameters for filtering the VCF file
python3 filter_variants.py --vcf $vcf --cov $cov --o filtered_vcf

# EXAMPLE: if running the python script using different parameters from the default ones
python3 filter_variants.py --vcf $vcf --cov $cov --q 20 --dp 8 --gq 30 --o filtered_vcf
