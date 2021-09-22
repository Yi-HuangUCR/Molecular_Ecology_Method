#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=18G
#SBATCH --time=2-00:00:00     # 1 hour
#SBATCH --output=freebayes.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="freebayesSNPcalling"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short


# Print current date
date

# Load software
module load freebayes

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

# SNP calling
freebayes -f /bigdata/entm249/yhuan073/RAD-seq2018/all_joined-SR_clipped_passed-re-filter_filtered-clusters_THIS_IS_THE_RAD_REFERENCE_FILE.fasta --ploidy 4 All_sorted_matefixed.bam > freebayes.vcf

# Print name of node
hostname
