#!/bin/bash


# Author : Chiara Bortoluzzi


#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=deepvariant
#SBATCH --output=output_%J
#SBATCH --error=error_%J


export OMP_NUM_THREADS=20



if [ $# -ne 3 ]
then
	echo -e "\nusage: `basename $0` <reference_genome> <bam> <model>\n"
	echo -e "DESCRIPTION: Use DeepVariant to call variants from sequencing reads\n\n"
	echo -e "INPUT:       <reference_genome>    The reference genome of the species of interest"
	echo -e "             <bam>                 Aligned sequences in BAM format"
  	echo -e "             <model>               Model that represents the type of data you are analysing (use WGS for Illumina whole-genome sequencing data, PACBIO for PacBio data, ONT_R104 for Oxford Nanopore R10.4.1 chemistry Simplex and Duplex reads)\n"

	echo -e "OUTPUT:      <vcf>                 A VCF file with called variants"
	echo -e "             <html>                A VCF stats report\n"

	echo -e "REQUIRES: The script requires deepvariant (v1.6.0)"
	exit
fi


ref=$1
bam=$2
model=$3

export TMPDIR="/tmp"

# Make sure that the reference genome and the bam file are indexed

mkdir -p deepvariant

vcf=deepvariant/"`basename $bam .bam`.deepvariant"

if [[ ! -f "$vcf" ]]; then
	echo -e "Run DeepVariant on $ref and $bam ... this is going to take some time\n\n"
	singularity run -B /cluster/work/pausch/cbortoluzzi/cetaceans deepvariant_latest.sif /opt/deepvariant/bin/run_deepvariant \
	--model_type $model \
	--ref=$ref \
	--reads=$bam \
	--output_vcf=$vcf.vcf.gz \
	--output_gvcf=$vcf.g.vcf.gz \
	--num_shards=16 \
	--intermediate_results_dir /tmp
fi
