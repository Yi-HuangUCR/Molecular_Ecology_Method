#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-02:00:00     # 1 hour
#SBATCH --output=vcftool.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="vcftool"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short


# Print current date
date


# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

module load vcftools
vcftools -vcf samtoolsLGC.vcf -012
paste out.012.indv out.012 > PCAinput.txt

# Print name of node
hostname
