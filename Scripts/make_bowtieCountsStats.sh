
## Code to make bowtie count stats for depth correction in mapping pipeline
# write the numbers in a file (folder : ${outdir}/bowtie_wholeGenome/)

logfolder=${1}
ext=${2}

for file in $(ls ${logfolder}| grep ${ext})
do
base=$(basename ${file} ${ext})
echo ${base} >> ${logfolder}/filenames.txt
awk 'NR > 5 && NR < 10 {print $1}' ${logfolder}/${file} | tr '[:space:]' '\t' >> ${logfolder}/readCount.txt
echo "" >> ${logfolder}/readCount.txt
done

echo -e "sample"'\t'"totalreads"'\t'"unmapped"'\t'"aligned.once"'\t'"aligned.more"'\t' > ${logfolder}/MappedReadCounts.txt
paste ${logfolder}/filenames.txt ${logfolder}/readCount.txt >> ${logfolder}/MappedReadCounts.txt
rm ${logfolder}/filenames.txt ${logfolder}/readCount.txt
