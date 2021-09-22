#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=62
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-02:00:00     # 1 hour
#SBATCH --output=samtobam.stdout
#SBATCH --mail-user=yhuan073@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="samtobam"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short


# Print current date
date

# Load software
module load samtools
#module load bowtie2

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

#change sam to bam
samtools view -S -b *.sam -o *.bam -@ 63

# Concatenate BAMs
#samtools cat -h header.sam -o out.bam in1.bam in2.bam

#bowtie paired reads
#bowtie2 -x ~/bigdata/bowtieDB/working/f_selysi_v02 -1 KYOC2W1.pear.unassembled.forward.fastq \
#-2 KYOC2W1.pear.unassembled.reverse.fastq --ff --very-sensitive | \
#samtools import ~/bigdata/bowtieDB/working/f_selysi_v02.fasta.fai - - | samtools sort - KYOC2W1.pair 

#bowtie single reads
#bowtie2 -x ~/bigdata/bowtieDB/working/f_selysi_v02 -U KYOC2W1.single.fq.gz --very-sensitive | \
#samtools import ~/bigdata/bowtieDB/working/f_selysi_v02.fasta.fai - - | samtools sort - KYOC2W1.single 

#merge the two bams
#samtools merge KYOC2W1.all.bam KYOC2W1.pair.bam KYOC2W1.single.bam -@ 2
#rm KYOC2W1.pair.bam KYOC2W1.single.bam 

#index the final bam
#samtools index KYOC2W1.all.bam

# Print name of node
hostname
