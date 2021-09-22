#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-02:00:00     # 1 hour
#SBATCH --output=stepone.stdout
#SBATCH --mail-user=yhuan073@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="stepone"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short

grep "#" freebayes.vcf | tail --l 1 | cut -c 2- | tr '-' '_' > header_of_your_vcf_file.txt

grep "TYPE=snp;" freebayes.vcf > biallelic_snps.txt

grep "TYPE=snp,snp;" freebayes.vcf > triallelic_snps.txt

grep "TYPE=snp,snp,snp;" freebayes.vcf > fourallelic_snps.txt

cat header_of_your_vcf_file.txt biallelic_snps.txt triallelic_snps.txt fourallelic_snps.txt > output_from_step_one.txt

