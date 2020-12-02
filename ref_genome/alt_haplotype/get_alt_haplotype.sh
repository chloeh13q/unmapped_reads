#!/bin/bash
#SBATCH --job-name=get_alt_haplotype
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/chloehe/logs/get_alt_haplotype.out
#SBATCH --error=/scratch/groups/dpwall/personal/chloehe/logs/get_alt_haplotype.err
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chloehe@stanford.edu
#SBATCH --mem=50G
######SBATCH --cpus-per-task=1
job_directory="${MY_SCRATCH}/unmapped_reads/ref_genome/alt_haplotype"
python3 "${job_directory}/list_accn.py" MH533022 MH534863
python3 "${job_directory}/get_fasta_by_accn.py" "accn.txt" "fasta_by_accn/"
echo "job finished, all alternative haplotypes downloaded"
