#!/bin/bash
#SBATCH --job-name=index-fa
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/chloehe/logs/index_fa.out
#SBATCH --error=/scratch/groups/dpwall/personal/chloehe/logs/index_fa.err
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chloehe@stanford.edu
#SBATCH --mem=50G
######SBATCH --cpus-per-task=1
dir=/scratch/groups/dpwall/personal/chloehe/unmapped_reads/ref_genome
echo "Indexing human viruses reference genome..."
bwa index "$dir/virus/human_viruses.fasta"
echo "Indexing bacteria reference genome..."
bwa index "$dir/bacteria/all_seqs.fa"
echo "Indexing bacteriophages reference genome..."
bwa index "$dir/bacteriophage/bacteriophages.fasta"
echo "Indexing alternative haplotypes reference genome..."
bwa index "$dir/alt_haplotype/alt_haplotypes.fa"
echo "All indexing completed"
