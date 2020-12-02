#!/bin/bash
#SBATCH --job-name=ihart_1
#SBATCH --array=1-192%10
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/chloehe/logs/ihart_1_%a.out
#SBATCH --error=/scratch/groups/dpwall/personal/chloehe/logs/ihart_1_%a.err
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chloehe@stanford.edu
#SBATCH --mem=40G
######SBATCH --cpus-per-task=1

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Start processing task ID: ${SLURM_ARRAY_TASK_ID}"

# Get batch ID and job number of the current job
num_samples=(96 192 288 383 477 567 655 745 836 925 1020 1115 1208 1302 1397 1399 1487 1581 1662 1749 1844 1939 2025 2120 2211 2249 2295 2386 2473 2558 2644 2707 2767 2858 2953 3046 3140 3235 3325 3418 3513 3603 3698 3793 3888 3978 4071 4165 4253 4272 4343 4409)
batch_ids=(batch_00009 batch_00010 batch_00011 batch_00012 batch_00013 batch_00014 batch_00015 batch_00016 batch_00017 batch_00018 batch_00024 batch_00025 batch_00026 batch_00027 batch_00028 batch_00516 batch_00924 batch_00925 batch_00926 batch_00927 batch_00928 batch_00929 batch_00930 batch_00931 batch_00932 batch_00933 batch_00935 batch_01002 batch_01003 batch_01004 batch_01005 batch_01006 batch_01007 batch_01008 batch_01009 batch_01010 batch_01011 batch_01012 batch_01014 batch_01015 batch_01016 batch_01017 batch_01018 batch_01019 batch_01020 batch_01021 batch_01022 batch_01024 batch_01025 batch_01372 batch_01378 batch_01379)
for i in "${!num_samples[@]}"; do
   if [[ "${num_samples[$i]}" -ge "${SLURM_ARRAY_TASK_ID-1}" ]]; then
       if [[ $i -eq 0 ]]; then
           job=$(( ${SLURM_ARRAY_TASK_ID-1} ))
       else
           job=$(( ${SLURM_ARRAY_TASK_ID-1}-num_samples[$(($i-1))] ))
       fi
       batch="${batch_ids[$i]}"; break;
   fi
done

echo "${batch} - ${job}"
. /scratch/groups/dpwall/personal/chloehe/unmapped_reads/batch/batch_info.sh 
arr=(${!batch})

# Retrieve current sample from AWS
HOST=${arr[job]}
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Retrieving sample ID $HOST from AWS..."

dir="/scratch/users/chloehe/unmapped_reads/bam/${batch}"
mkdir -p $dir
fullpath="${dir}/${HOST}.final"
aws s3 cp s3://ihart-hg38/cram/$HOST.final.cram "${fullpath}.cram"
aws s3 cp s3://ihart-hg38/cram/$HOST.final.cram.crai "${fullpath}.cram.crai"

ml samtools
ml biology

outpath="${dir}/${HOST}"
mkdir -p $outpath
#mkdir -p ${outpath}/flagstat
prefix=$(echo ${fullpath} | awk -F "/" '{print $NF}')
ID=$(echo ${outpath} | awk -F "/" '{print $NF}')

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Preprocessing file for sample ID ${ID}..."

# Remove non-primary and duplicate reads
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Removing non-primary and duplicate reads..."
samtools view -b -@ 8 -h -F 1280 "${fullpath}.cram" -o "${outpath}/${prefix}.filtered.bam"
#samtools flagstat "${outpath}/${prefix}.filtered.bam" > "${outpath}/flagstat/${prefix}.filtered.flagstat"

# Split primary and supplementary reads
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Splitting primary and supplementary reads..."
samtools view -b -@ 8 -h -F 2048 "${outpath}/${prefix}.filtered.bam" -o "${outpath}/${prefix}.primary.bam"

# Extract mapped/unmapped reads
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Extracting mapped and unmapped reads..."
samtools view -@ 8 -h -f 1 -F 268 -bo "${outpath}/${prefix}.map_map.bam" "${outpath}/${prefix}.primary.bam"
samtools view -@ 8 -f 4 -F 264 -h -bo "${outpath}/${prefix}.unmap_map.bam" "${outpath}/${prefix}.primary.bam"
samtools view -@ 8 -f 8 -F 260 -h -bo "${outpath}/${prefix}.map_unmap.bam" "${outpath}/${prefix}.primary.bam"
samtools view -@ 8 -f 12 -F 256 -h -bo "${outpath}/${prefix}.unmap_unmap.bam" "${outpath}/${prefix}.primary.bam"

# Extract reads that are improperly paired and split these reads into high MAPQ and low MAPQ buckets
# Filtering out properly paired reads will leave improperly paird reads + UM + UU
# So need to use 0x4 flag here to filter out "read unmapped" (UM + UU) and deal with these reads in the next step
#samtools view -h -F 262 "${outpath}/${prefix}.primary.bam" |
#tee >(awk '{if($1 ~ /^@/) {print $0} else if($5<10) {print $0}}' |
#samtools view -Sbo "${outpath}/${prefix}.improper_lowq.bam" -) >(awk '{if($1 ~ /^@/) {print $0} else if($5>=10) {print $0}}' |
#samtools view -Sbo "${outpath}/${prefix}.improper_highq.bam" -) >/dev/null
#samtools flagstat "${outpath}/${prefix}.improper_lowq.bam" > "${outpath}/${prefix}.improper_lowq_flagstat"
#samtools flagstat "${outpath}/${prefix}.improper_highq.bam" > "${outpath}/${prefix}.improper_highq_flagstat"

# Instead of filtering by MAPQ, filter by alignment score
# Safer and more efficient to extract directly from map_map.bam as both ends need to be mapped in order to be considered improperly paired
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Filtering reads with low alignment scores..."
samtools view -@ 8 -h -F 262 "${outpath}/${prefix}.map_map.bam" |
tee >(egrep '^@|AS:i:1[0-9]{2,}' |
samtools view -@ 8 -Sbo "${outpath}/${prefix}.improper_highq.bam" -) >(grep -v 'AS:i:1[0-9]{2,}' |
samtools view -@ 8 -Sbo "${outpath}/${prefix}.improper_lowq.bam" -) >/dev/null

# Separate low score reads into paired-end and single-end
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Splitting low score reads..."
samtools sort -n -@ 8 -m 4G "${outpath}/${prefix}.improper_lowq.bam" > "${outpath}/${prefix}.improper_lowq.sorted.bam"
samtools view -@ 8 "${outpath}/${prefix}.improper_lowq.sorted.bam" |
cut -f1 |
uniq -c |
sed -e 's/^[ \t]*//' |
tee >(grep '1\s' > "${outpath}/${prefix}.improper_lowq.single") >(grep '2\s' > "${outpath}/${prefix}.improper_lowq.paired") >/dev/null
samtools view -h -@ 8 ${outpath}/${prefix}.improper_lowq.sorted.bam |
tee >(awk 'FNR==NR{a[$2];next}{if($1 ~ /^@/) {print $0} else if($1 in a) {print $0}}' ${outpath}/${prefix}.improper_lowq.paired - |
samtools view -@ 8 -bo "${outpath}/${prefix}.improper_lowq.paired.bam" -) >(awk 'FNR==NR{a[$2];next}{if($1 ~ /^@/) {print $0} else if($1 in a) {print $0} else if($1 in a) {print $0}}' ${outpath}/${prefix}.improper_lowq.single - |
samtools view -@ 8 -bo "${outpath}/${prefix}.improper_lowq.single.bam" -) >/dev/null

# Sort reads by name and merge
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Sorting reads..."
samtools sort -n -@ 8 -m 4G "${outpath}/${prefix}.unmap_map.bam" > "${outpath}/${prefix}.unmap_map.sorted.bam"
samtools sort -n -@ 8 -m 4G "${outpath}/${prefix}.unmap_unmap.bam" > "${outpath}/${prefix}.unmap_unmap.sorted.bam"
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Merging reads..."
samtools merge -f -n -@ 8 "${outpath}/${prefix}.single_to_aln.bam" "${outpath}/${prefix}.improper_lowq.single.bam" "${outpath}/${prefix}.unmap_map.sorted.bam"
samtools merge -f -n -@ 8 "${outpath}/${prefix}.paired_to_aln.bam" "${outpath}/${prefix}.improper_lowq.paired.bam" "${outpath}/${prefix}.unmap_unmap.sorted.bam" 

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Done extracting reads for realignment"

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Start realigning..."
ref=/scratch/groups/dpwall/personal/chloehe/unmapped_reads/ref_genome
echo "Realign to reference genome at ${ref}"

# Send paired-end reads for realignment
samtools fastq -@ 8 "${outpath}/${prefix}.paired_to_aln.bam" |
bwa mem -t 24 -p -T 100 $ref/combined.fa - |
samtools sort -@ 8 -m 4G - |
samtools view -@ 8 -h -b -o "${outpath}/${prefix}.paired.aln_all.bam" -
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Done realigning paired-end reads"

# Send single-end reads for realignment
samtools fastq -@ 8 "${outpath}/${prefix}.single_to_aln.bam" |
bwa mem -t 24 -T 100 $ref/combined.fa - |
samtools sort -@ 8 -m 4G - |
samtools view -@ 8 -h -b -o "${outpath}/${prefix}.single.aln_all.bam" -
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Done realigning single-end reads"

# Generate alignment summary table
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Generating alignment summary table..."
job_directory=/scratch/groups/dpwall/personal/chloehe/unmapped_reads/scripts_fam
module load py-numpy/1.18.1_py36
python3 "${job_directory}/gen_alignment_table.py" "${outpath}"

# Remove intermediate files
rm "${fullpath}.cram" "${fullpath}.cram.crai" "${outpath}/${prefix}.filtered.bam" "${outpath}/${prefix}.primary.bam" "${outpath}/${prefix}.map_map.bam" "${outpath}/${prefix}.unmap_map.bam" "${outpath}/${prefix}.unmap_unmap.bam" "${outpath}/${prefix}.single_to_aln.bam" "${outpath}/${prefix}.paired_to_aln.bam"
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Task ID ${SLURM_ARRAY_TASK_ID} is complete"
