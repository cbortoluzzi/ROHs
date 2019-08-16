#!/usr/bash/env python


# Author: Chiara Bortoluzzi
# Date: 15 August 2019 (last updates)
# Calculate heterozygosity in non-overlapping consecutive windows (e.g. 10kb)


# Imports
from __future__ import division
import argparse
import vcf
import os
import sys
import subprocess
import json
from collections import OrderedDict


parser = argparse.ArgumentParser(description="Calculate heterozygosity in non-overlapping consecutive windows along the genome")
parser.add_argument("-f", help="Fasta file with chromosomes and length")
parser.add_argument("-j",help="Path to the json file")
parser.add_argument("-b", help="Binsize in base pairs",type=int)



def is_number(x):
        try:
            	xf = float(x)
                return True
        except:
               	return False



def chromosome_length(fasta):
        mychr = OrderedDict()
        for line in fasta:
                info = line.split("\t")
                if is_number(info[0]):
                        chr_num,length = info[0],info[1]
                        mychr[chr_num] = length
        return mychr



def get_samtools_output(read_json,chr_num,start,end):
        bam_path = read_json['bam_path']
        vcf_path = read_json['vcf_path']
        min_cov_indiv = read_json['min_cov']
        max_cov_indiv = read_json['max_cov']
        indiv = read_json['sample']
        input_roh = read_json['output']
        with open(input_roh, "a") as roh:
                nsites = 0
                cov_sites = 0
                nhet = 0
                vcf_name = open(vcf_path,'rb')
                vcf_reader = vcf.Reader(vcf_name)
                for record in vcf_reader.fetch(chr_num,start,end):
                        cov = record.INFO['DP']
                        if cov >= min_cov_indiv and cov <= max_cov_indiv:
                                nsites += 1
                        nhet += int(record.num_het)
                mysam_sites = {}
                command = "samtools depth -r %d:%d-%d %s | awk '($3 >= %d && $3 <= %d){{print}}'" % (chr_num,start,end,bam_path,min_cov_indiv,max_cov_indiv)
                try:
                    	process = subprocess.check_output(command, shell = True)
                except:
                       	raise ValueError("samtools failed, sucka!")
                output_samtools = process.decode().split("\n")
                for site in output_samtools:
                        site = site.strip().split()
                        if site == ['']:
                                chrom,pos,depth = site[0:4]
                                mysam_sites[indiv,chr_num,start,end] = 0,0
                        else:
                             	cov_sites += 1
                                mysam_sites[indiv,chr_num,start,end] = cov_sites,nhet
                for k,v in mysam_sites.items():
                        roh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(k[0],str(k[1]),str(k[2]),str(k[3]),str(v[0]),str(v[1])))


def parameters_samtools(mychr,binsize,read_json):
        for chr_num in mychr:
                for i in range(0,int(mychr[str(chr_num)]),binsize):
                        get_samtools_output(read_json,int(chr_num),i,i+binsize)



if __name__ == "__main__":
        args = parser.parse_args()
        fasta = open(args.f)
        json_file = open(args.j)
        read_json = json.load(json_file)
        chromosomes_length = chromosome_length(fasta)
        param_samtools = parameters_samtools(chromosomes_length,args.b,read_json)
        fasta.close()


