#!/bin/bash
#SBATCH --job-name=fam1_align_all
#SBATCH --array=1-5
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/chloehe/logs/fam1_align_all_%a.out
#SBATCH --error=/scratch/groups/dpwall/personal/chloehe/logs/fam1_align_all_%a.err
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chloehe@stanford.edu
#SBATCH --mem=50G
######SBATCH --cpus-per-task=1

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Start processing task ID: ${SLURM_ARRAY_TASK_ID}"

dir=/scratch/groups/dpwall/personal/chloehe/unmapped_reads/bam/fam1
ml samtools
ml biology
i=0
fullpath=''
for file in ${dir}/*.cram
do
  i=$(( i + 1 ))
  if [ $SLURM_ARRAY_TASK_ID -eq $i ]; then
    fullpath=${file%.cram}
    echo "Processing ${file}"
    break
  fi
done

outpath=$(echo ${fullpath} | cut -d. -f1)
mkdir -p $outpath
mkdir -p ${outpath}/flagstat
prefix=$(echo ${fullpath} | awk -F "/" '{print $NF}')
ID=$(echo ${outpath} | awk -F "/" '{print $NF}')

#samtools flagstat "${fullpath}.cram" > "${outpath}/${prefix}_flagstat"

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Preprocessing file for sample ID ${ID}..."

# Remove non-primary and duplicate reads
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Removing non-primary and duplicate reads..."
samtools view -b -@ 8 -h -F 1280 "${fullpath}.cram" -o "${outpath}/${prefix}.filtered.bam"
samtools flagstat "${outpath}/${prefix}.filtered.bam" > "${outpath}/flagstat/${prefix}.filtered.flagstat"

# Split primary and supplementary reads
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Splitting primary and supplementary reads..."
samtools view -b -@ 8 -h -F 2048 "${outpath}/${prefix}.filtered.bam" -o "${outpath}/${prefix}.primary.bam" -U "${outpath}/${prefix}.supplementary.bam"

# Subsampling
#samtools view -bs 42.1 "$rawinfile.bam" > "$infile.bam"

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
tee >(egrep '^@|AS:i:1[0-9]*' |
samtools view -@ 8 -Sbo "${outpath}/${prefix}.improper_highq.bam" -) >(grep -v 'AS:i:1[0-9]*' |
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

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Sorting reads..."
samtools sort -n -@ 8 -m 4G "${outpath}/${prefix}.unmap_map.bam" > "${outpath}/${prefix}.unmap_map.sorted.bam"
samtools sort -n -@ 8 -m 4G "${outpath}/${prefix}.unmap_unmap.bam" > "${outpath}/${prefix}.unmap_unmap.sorted.bam"
echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Merging reads..."
samtools merge -f -n -@ 8 "${outpath}/${prefix}.single_to_aln.bam" "${outpath}/${prefix}.improper_lowq.single.bam" "${outpath}/${prefix}.unmap_map.sorted.bam"
samtools merge -f -n -@ 8 "${outpath}/${prefix}.paired_to_aln.bam" "${outpath}/${prefix}.improper_lowq.paired.bam" "${outpath}/${prefix}.unmap_unmap.sorted.bam" 
#rm "${outpath}/${prefix}.lowq.paired.bam" "${outpath}/${prefix}.lowq.single.bam" "${outpath}/${prefix}.unmap_map.bam" "${outpath}/${prefix}.map_unmap.bam" "${outpath}/${prefix}.unmap_unmap.bam"

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Done extracting reads for realignment"

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Start realigning..."
ref=/scratch/groups/dpwall/personal/chloehe/unmapped_reads/ref_genome
echo "Realignment to reference genome at ${ref}"

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

echo "[`date "+%Y-%m-%d %H:%M:%S"`]" "Task ID ${SLURM_ARRAY_TASK_ID} is complete"
