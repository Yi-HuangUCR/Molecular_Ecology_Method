#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-02:00:00     # 1 hour
#SBATCH --output=bowtieloop.stdout
#SBATCH --mail-user=yhuan073@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="samtool_sort"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short
module load samtools
i=$SLURM_ARRAY_TASK_ID
Samtools index sample_${i}.bam
