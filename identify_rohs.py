#!/usr/bin/env python 


# Chiara Bortoluzzi
# Date: 15 August 2019 (last updates)
# Identify runs of homozygosity (ROHs) using information on heterozygosity calculated in non-overlapping consecutive windows along the genome


# Imports
from __future__ import division
from sys import argv
from collections import OrderedDict
import argparse



parser = argparse.ArgumentParser(description = "Identify runs of homozygosity	(ROHs) along the genome")
parser.add_argument("-i", help = "Input file with information on heterozygosity calculated in non-overlapping windows")
parser.add_argument("-b", help = "Window size in base pairs", type=int)
parser.add_argument("-t1", help = "Minimum number of well-covered sites", type=int)
parser.add_argument("-t2", help = "Threshold to filter the snp count within a stretch", type=float)
parser.add_argument("-t3", help = "Threshold to filter the average diversity within a run", type=float)
parser.add_argument("-t4", help = "Number of consecutive bins", type=int)
parser.add_argument("-t5", help = "Threshold for local maximum", type=int)
parser.add_argument("-roh", help = "Output file where to save the runs of homozygosity")
parser.add_argument("-nonroh", help = "Output file containing the nucleotide diversity outside ROHs")




def is_number(x):
	try:
		xf = float(x)
		return True
	except:
		return False



def genome_nucl_div(data,binsize,thresh1):
	mygw = {}
	n_sites = 0
	n_snp = 0
	for line in data:
		sample,chrom,start,end,nsites,het = line
		nsites = int(nsites)
		het = int(het)
		prop_het = round((binsize/nsites)*het,2)
		if nsites > thresh1:
			n_sites += 1
			n_snp += prop_het
	mygw['total_number_well_covered_sites'] = n_sites
	mygw['sum_of_well_covered_snps'] = n_snp
	mygw['average_genome_wide_nucl_div'] = n_snp/n_sites
	return mygw




def filter_roh(data,mygw,binsize,thresh1,thresh4):
	mybins = []
	mysnp = []
	myroh = OrderedDict()
	mynonroh = OrderedDict()
	pre_chrom = 1
	for line in data:
		cur_chrom = int(line[1])
		if cur_chrom == pre_chrom:
			sample,chrom,start,end,nsites,het = line
			nsites = int(nsites)
			het = int(het)
			prop_het = round((binsize/nsites)*het,2)
			mybins.append([line,prop_het]) 
			if len(mybins) == thresh4:
				mybins,mysnp,myroh,mynonroh = extract_consecutive_bins(line,mygw,mybins,mysnp,myroh,mynonroh,binsize,thresh1,thresh4)
		else:
			mybins,mysnp,myroh,mynonroh = extract_consecutive_bins(line,mygw,mybins,mysnp,myroh,mynonroh,binsize,thresh1,thresh4)
			sample,chrom,start,end,nsites,het = line
			nsites = int(nsites)
			het = int(het)
			prop_het = round((binsize/nsites)*het,2)
			mybins.append([line,prop_het])
			pre_chrom = cur_chrom
	return myroh,mynonroh



def extract_consecutive_bins(line,mygw,mybins,mysnp,myroh,mynonroh,binsize,thresh1,thresh4):	
	for item in mybins:
		# Check that the number of sites is greater to 2000
		if int(item[0][4]) > thresh1:
			mysnp.append(item[1])
	if len(mysnp) != 0:
		avg_bin_div = round(sum(mysnp)/len(mysnp),2)    
	for item in mybins:
		if int(item[0][4]) > thresh1:
			# Check that the average SNPs diversity in the 10 consecutive bins is < or equal to the genome-wide average
			if avg_bin_div < mygw['average_genome_wide_nucl_div']:
				values = "\t".join(map(str,item[0][0:5]))
				myroh[values] = item[1]
			else:
				values_non_roh = "\t".join(map(str,item[0][0:5]))
				mynonroh[values_non_roh] = item[1]
		else:
			values = "\t".join(map(str,item[0][0:5]))
			rm_site = "removed_site"
			myroh[values] = rm_site

	mybins = []
	mysnp = []
	return mybins,mysnp,myroh,mynonroh



def nucleotide_diversity_outROH(mynonroh,output_nonroh):
	for key,value in mynonroh.items():
		sample,chrom,start,end,nsites = key.split("\t")
		output_nonroh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample,chrom,start,end,nsites,value))



def concatenate_roh(myroh,binsize,mygw,thresh2,thresh3,thresh5,output_roh):
	prev_bin = 0
	consecutive_bins = []
	sublist = []
	snp_count_roh = []
	final_roh = OrderedDict()
	for key, value in myroh.items():
		sample,chrom,start,end,nsites = key.split()
		cur_bin = int(end)
		if cur_bin == prev_bin + binsize or not sublist:
			sublist.append((sample,chrom,start,end,value))
			prev_bin = cur_bin
		else:
			consecutive_bins.append(sublist)
			sublist = []
			sublist.append((sample,chrom,start,end, value))
			prev_bin = cur_bin
	for row in consecutive_bins:
		snp_counts = []
		num_removed_sites = []
		for element in row:
			if element[4] == "removed_site":
				pass
			else:
				if float(element[4]) < thresh2*mygw['average_genome_wide_nucl_div'] or float(element[4]) < thresh5*mygw['average_genome_wide_nucl_div']:
					snp_counts.append(float(element[4]))	
		if len(snp_counts) != 0:
			avg_snp_counts_roh = round(sum(snp_counts)/len(snp_counts),2)
			if avg_snp_counts_roh < thresh3*mygw['average_genome_wide_nucl_div']:
				list_roh = [i for i in row if is_number(i[4]) < thresh5*mygw['average_genome_wide_nucl_div']]
				peak_sites = [i for i in row if is_number(i[4]) > thresh5*mygw['average_genome_wide_nucl_div']]
				removed_sites = [i for i in row if i[4] == "removed_site"]
				sample_id = list_roh[0][0]
				chrom = list_roh[0][1]
				start_bin = list_roh[0][2]
				end_bin = list_roh[-1][3]
				size_roh = len(list_roh)
				length_roh = len(list_roh)+len(removed_sites)
				output_roh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_id,chrom,start_bin,end_bin,size_roh,length_roh,avg_snp_counts_roh))





if __name__ == "__main__":
	args = parser.parse_args()
	output_roh = open(args.roh,"w")
	output_nonroh = open(args.nonroh,"w")
	with open(args.i) as source:
		data = source.readlines()
		data = [line.strip().split() for line in data]
	avg_gw_nucl_div = genome_nucl_div(data, args.b, args.t1)
	myroh, mynonroh = filter_roh(data, avg_gw_nucl_div, args.b, args.t1, args.t4)
	avg_nucl_outROH = nucleotide_diversity_outROH(mynonroh, output_nonroh)
	merge_filtered_roh = concatenate_roh(myroh, args.b, avg_gw_nucl_div, args.t2, args.t3, args.t5, output_roh)
	output_roh.close()
	output_nonroh.close()
	

