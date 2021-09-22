#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=02:00:00     
#SBATCH --output=build_lib.stdout
#SBATCH --mail-user=yhuan073@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="build_lib"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu


# Print current date
date

# Load modules
module load bowtie2

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

# Main works

bowtie2-build all_joined-SR_clipped_passed-re-filter_filtered-clusters_THIS_IS_THE_RAD_REFERENCE_FILE.fasta LGCref

# Print name of node
hostname

