#!/bin/bash
#SBATCH --job-name=batch_table
#SBATCH --array=1-54%10
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/chloehe/logs/batch_table_%a.out
#SBATCH --error=/scratch/groups/dpwall/personal/chloehe/logs/batch_table_%a.err
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chloehe@stanford.edu
#SBATCH --mem=40G
######SBATCH --cpus-per-task=1

batch_ids=(batch_00009 batch_00010 batch_00011 batch_00012 batch_00013 batch_00014 batch_00015 batch_00016 batch_00017 batch_00018 batch_00024 batch_00025 batch_00026 batch_00027 batch_00028 batch_00514 batch_00516 batch_00924 batch_00925 batch_00926 batch_00927 batch_00928 batch_00929 batch_00930 batch_00931 batch_00932 batch_00933 batch_00935 batch_01002 batch_01003 batch_01004 batch_01005 batch_01006 batch_01007 batch_01008 batch_01009 batch_01010 batch_01011 batch_01012 batch_01013 batch_01014 batch_01015 batch_01016 batch_01017 batch_01018 batch_01019 batch_01020 batch_01021 batch_01022 batch_01024 batch_01025 batch_01372 batch_01378 batch_01379)
batch="${batch_ids[$(( ${SLURM_ARRAY_TASK_ID-1}-1 ))]}"
echo ${batch}
module load py-numpy/1.18.1_py36
python3 gen_batch_table.py /scratch/users/chloehe/unmapped_reads/bam/${batch}
