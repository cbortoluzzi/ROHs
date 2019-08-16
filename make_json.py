#!/usr/bin/env python


# Author: Chiara Bortoluzzi
# Date: 15 August 2019 (last updates)
# Generate for each individual a json file with information on BAM, VCF, minimum, and maximum coverage per SNP
# The json file will be used as input file in the ROH pipeline


# Imports
import argparse
import os
import json
import re


parser = argparse.ArgumentParser(description='Generate for each individual a json file with information on BAM, VCF, minimum, and maximum coverage')
parser.add_argument("-b", help="List of bam files, one per line")
parser.add_argument("-v", help="List of vcf files, one per line")
parser.add_argument("-c", help="Tab delimited file with coverage information (i.e. sample, genome-wide coverage)")
parser.add_argument("-m", help="Minimum coverage", type=int)


def make_dictionary(bam, vcf_input, coverage, min_cov):
        myinput = {}
        mybam = {}
        myvcf = {}
        ## Bam files list
        for line in bam:
                line = line.strip()
                if os.path.isfile(line):
                        info = line.split("/")
                        bam_name = info[-1]
                        bam_sample_name = re.split(r'[_.\s]\s*', bam_name)[0]
                        mybam[bam_sample_name] = line
                else:
                     	raise ValueError("BAM file does not exist")
        ## VCF files list
        for line in vcf_input:
                line = line.strip()
                if os.path.isfile(line):
                        info = line.split("/")
                        vcf_name = info[-1]
                        vcf_sample_name = re.split(r'[_.\s]\s*', vcf_name)[0]
                        myvcf[vcf_sample_name] = line
                else:
                     	raise ValueError("VCF file does not exist")
        ## Table of coverage file
        for line in coverage:
                individual, avg_coverage = line.strip().split()
                avg_coverage = avg_coverage.replace('X','')
                avg_coverage = float(avg_coverage)
                min_cov = min_cov
                max_cov = int(avg_coverage*2)
                myinput[individual] = [avg_coverage,min_cov,max_cov,mybam[individual],myvcf[individual]]
        return myinput


def make_json(myinput,indiv):
        avg_cov = myinput[indiv][0]
        min_cov = myinput[indiv][1]
        max_cov = myinput[indiv][2]
        bam_path = myinput[indiv][3]
        vcf_path = myinput[indiv][4]
        output = open(indiv+".json","w")
        data = {'sample': indiv, 'bam_path': bam_path,'vcf_path':vcf_path, 'avg_cov': avg_cov, 'min_cov': min_cov, 'max_cov': max_cov, 'output': "".join([indiv,"_input.txt"])}
        json.dump(data,output,indent=2)
        output.write("\n")
        output.close()


def write_json(myinput):
        for indiv in myinput.keys():
                make_json(myinput,indiv)



if __name__ == "__main__":
        args = parser.parse_args()
        bam, vcf_input, coverage = open(args.b), open(args.v), open(args.c)
        make_dict = make_dictionary(bam, vcf_input, coverage, args.m)
        write_json_file = write_json(make_dict)
        bam.close(), vcf_input.close(), coverage.close()


