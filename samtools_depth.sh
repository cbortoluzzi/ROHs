#!/bin/bash


# Author : Chiara Bortoluzzi


#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=deepvariant
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 1 ]
then
	echo -e "\nusage: `basename $0` <bam>\n"
	echo -e "DESCRIPTION: Calculate genome-wide coverage from a bam file\n\n"
	echo -e "INPUT:     <bam>       Aligned sequences in BAM format\n"

	echo -e "OUTPUT:    <cov>       A text file with the calculated genome-wide coverage\n"

	echo -e "REQUIRES:  The script requires samtools (v1.22.1)"
	exit
fi



bam=$1


samtools depth $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $bam".cov"


