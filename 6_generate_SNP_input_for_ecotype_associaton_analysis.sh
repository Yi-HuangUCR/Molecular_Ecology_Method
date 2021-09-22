#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-02:00:00     # 1 hour
#SBATCH --output=vcf.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="vcf"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short


# Print current date
date

# Load software
module load vcftools

#Filter missingdata
vcftools --vcf samtoolLGC.vcf --012  

# Print name of node
hostname
