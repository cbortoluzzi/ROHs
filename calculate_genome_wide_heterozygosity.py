#!/usr/bin/env python



# Author : Chiara Bortoluzzi



import vcf
import argparse
import subprocess
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Calculate genome-wide heterozygosity using a sliding window approach')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--bam', help = 'BAM file')
parser.add_argument('--cov', help = 'Genome-wide coverage')
parser.add_argument('--dp', help = 'Minimum read depth to retain a variant [defulat = 6]', type = int, default = 6)
parser.add_argument('--w', help = 'Window size [default = 10000 bp]', type = int, default = 10000)
parser.add_argument('--o', help = 'Output directory')




def sequences_bam(bam_f):
	mygenome = {}
	# Get BAM index stats: we will retain only the chromosome name and its total sequence length (in bp)
	command = 'samtools idxstats %s | cut -f 1,2' %(bam_f)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			chromosome, length = line.strip().split()
			length = int(length)
			if chromsome.startswith('CM'):
					mygenome[chromosome] = length
			else:
				continue
	return mygenome



def average_genome_coverage(coverage):
	with open(coverage) as f:
		for line in f:
			line = line.split('=')
			avg_genome_coverage = line[1].replace(' ','')
			# Set maximum depth as 2 times the average genome coverage
			max_depth = 2 * float(avg_genome_coverage)
	return max_depth



def calculate_binned_heterozygosity(mygenome, window, min_depth, max_depth, bam_f, vcf_f, filename, path):
	for chromosome in mygenome:
		seq_length = mygenome[chromosome]
		for i in range(0, seq_length, window):
			start = i
			end = i + window
			bam_depth(chromosome, start, end, min_depth, max_depth, bam_f, vcf_f, window, filename, path)



def bam_depth(chromosome, start, end, min_depth, max_depth, bam_f, vcf_f, window, filename, path):
	cov_sites = 0
	# Obtain the read depth of each site from the BAM file
	command = 'samtools depth -r %s:%d-%d %s' %(chromosome, start, end, bam_f)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			chromosome, position, depth = line.strip().split()
			# Retain only sites that meet our minimum and maximum depth 
			if int(depth) >= min_depth and int(depth) <= max_depth:
				cov_sites += 1
	heterozygosity(chromosome, start, end, cov_sites, vcf_f, window, filename, path)



def heterozygosity(chromosome, start, end, cov_sites, vcf_f, window, filename, path):
	vcf_reader = vcf.Reader(filename=vcf_f)
	nhet = 0
	for record in vcf_reader.fetch(chromosome, start, end):
		# Count number of heterozygous sites in a window
		for call in record.samples:
			nhet += record.num_het
	if cov_sites == window + 1:
		cov_sites = window
	try:
		# Calculate corrected heterozygosity
		SNPcount = round((window / cov_sites) * nhet, 3)
	except ZeroDivisionError:
		SNPcount = 0.0
	output_f = Path(path, filename + '.heterozygosity.txt')
	with open(output_f, 'a') as out:
		out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, start, end, cov_sites, nhet, SNPcount))




if __name__ == "__main__":
	args = parser.parse_args()
	filename = Path(args.vcf).stem.replace('.vcf', '')
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	sequences = sequences_bam(args.bam)
	max_depth = average_genome_coverage(args.cov)
	binned_heterozygosity = calculate_binned_heterozygosity(sequences, args.w, args.dp, max_depth, args.bam, args.vcf, filename, args.o)


