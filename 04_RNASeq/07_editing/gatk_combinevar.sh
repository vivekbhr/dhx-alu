#!/usr/bin/bash

varcall=${1} #varcall_polyA-minus
outdir=${2} #editing

module load GATK/3.5

# combine varients CONTROL
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 \
-V:ctrl_1 ${varcall}/HEK_ctrl_1/call.filtered.vcf \
-V:ctrl_2 ${varcall}/HEK_ctrl_2/call.filtered.vcf \
-V:ctrl_3 ${varcall}/HEK_ctrl_3/call.filtered.vcf \
-V:ctrl_4 ${varcall}/HEK_ctrl_4/call.filtered.vcf \
-o ${outdir}/combined_vars_ctrl.vcf --minimumN 3 \
-genotypeMergeOptions UNIQUIFY \
-R genome/GRCh38.fa

# combine varients TEST
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 \
-V:test19_1 ${varcall}/HEK_s4019_1/call.filtered.vcf \
-V:test19_2 ${varcall}/HEK_s4019_2/call.filtered.vcf \
-V:test19_3 ${varcall}/HEK_s4019_3/call.filtered.vcf \
-V:test19_4 ${varcall}/HEK_s4019_4/call.filtered.vcf \
-V:test20_1 ${varcall}/HEK_s4020_1/call.filtered.vcf \
-V:test20_2 ${varcall}/HEK_s4020_2/call.filtered.vcf \
-V:test20_3 ${varcall}/HEK_s4020_3/call.filtered.vcf \
-V:test20_4 ${varcall}/HEK_s4020_4/call.filtered.vcf \
-o ${outdir}/combined_vars_test.vcf --minimumN 6 \
-genotypeMergeOptions UNIQUIFY \
-R genome/GRCh38.fa

# select only SNPs
for var in ctrl test
do
java -Xmx2g -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
-R genome/GRCh38.fa \
-T SelectVariants \
-V ${outdir}/combined_vars_${var}.vcf \
-o ${outdir}/snps_${var}.vcf \
--selectTypeToInclude SNP & done


# combine varients CONTROL+TEST
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 \
-V ${outdir}/combined_vars_ctrl.vcf -V ${outdir}/combined_vars_test.vcf \
-genotypeMergeOptions UNIQUIFY \
-R genome/GRCh38.fa -o ${outdir}/combined_vars_all.vcf

