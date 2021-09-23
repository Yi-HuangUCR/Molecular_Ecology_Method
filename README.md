# Molecular_Ecology_Method
Indepdent project for ENTM249
Proceess ddRAD-seq data of Arctostaphylos glandulosa, a plant species composed of multiple subspecies including two endangered ones. Scripts for reads alignment, SNPs calling & filtering, differentiation visualization via Principle Component analysis and ecotype-genotype association analysis.

Tips learned from this project
● Working with polyploid individuals will lead to limited choice of package for analysis (I used Freebayes here) and may require customized scripts for downstream analysis.
● A good reference genome with well-assembled scaffolds is important for downstream analysis. (The reference genome of a relative species led to many miscalling of SNPs)
● Bam files are easily truncated when the storage space in your working directory is limited.
● Submitting array jobs could process the data analysis of multiple samples efficiently and is helpful in population genetics study.
