#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00     # 1 hour
#SBATCH --output=samtoolSNPLGC.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="samtoolSNPcallLGC"
#SBATCH -p highmem # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short


# Print current date
date

# Load software
module load samtools/0.1.19

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

# SNP calling
samtools mpileup -q 20 -uDf /bigdata/entm249/yhuan073/RAD-seq2018/all_joined-SR_clipped_passed-re-filter_filtered-clusters_THIS_IS_THE_RAD_REFERENCE_FILE.fasta -b bam.txt | bcftools view -vcg - > samtoolLGC.vcf

# Print name of node
hostname
