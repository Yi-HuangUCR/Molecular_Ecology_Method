#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-02:00:00     # 1 hour
#SBATCH --output=bowtieloop.stdout
#SBATCH --mail-user=yhuan073@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="bowtieloop"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short



# Print current date
date

# Load software
#module load samtools
module load bowtie2/2.3.4.1

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR


#bowtie paired reads in loop
i=$SLURM_ARRAY_TASK_ID
bowtie2 -x /bigdata/entm249/yhuan073/RAD-seq2018/LGCref/LGCref -1 ${i}_R1_clipped_trimmed.fastq -2 ${i}_R2_clipped_trimmed.fastq \
-p 64 --ff --very-sensitive > \
/bigdata/entm249/yhuan073/RAD-seq2018/QualityTrimmed/samfile/sample_${i}.sam


# Print name of node
hostname
