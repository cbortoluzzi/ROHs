#!/bin/bash


# Author : Chiara Bortoluzzi


#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=RoH
#SBATCH --output=output_%J
#SBATCH --error=error_%J


het=$1


python runs_of_homozygosity.py --het $het --w 10000 --t1 10 --t2 6000 --t3 0.25 --o RoHs


