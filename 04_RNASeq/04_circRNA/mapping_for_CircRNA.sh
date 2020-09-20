#!/usr/bin/sh

## Remapping the polyA-/total seq samples using fastq extracted from split bam files
## Input fastq are trimmed 

mkdir ../circRNA_mapping
indir=../../02_autoMapping/analysis_results/polyA_minus/Trim_Galore
outdir=../circRNA_mapping

#hg38
hgff=/data/repository/organisms/GRCh38_ensembl/ensembl/release-78/genes.gtf
hgstar=/data/repository/organisms/GRCh38_ensembl/STARIndex

## MAP using STAR
module load slurm

for file in $(ls ${indir} | grep _R1.fastq.gz$)
do base=$(basename $file _R1.fastq.gz)
mkdir ${outdir}/star_${base}
echo "#!/usr/bin/sh" > ${base}_slurm.sh
echo "/package/STAR-2.4.2a/bin/STAR --runThreadN 15 --genomeDir ${hgstar} --readFilesCommand zcat --readFilesIn ${indir}/${base}_R1.fastq.gz ${indir}/${base}_R2.fastq.gz --sjdbGTFfile $hgff --chimSegmentMin 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outdir}/star_${base}/${base} " >> ${base}_slurm.sh
chmod 777 ${base}_slurm.sh
SlurmEasy -t 15 ${base}_slurm.sh
done


