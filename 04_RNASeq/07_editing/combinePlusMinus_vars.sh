#!/usr/bin/bash

module load GATK/3.5

## MERGE Vaiants which are filtered
# combine test uniq vars
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 -V:polyA_minus editing_polyA-minus/testUniq_filtered.vcf -V:polyA_plus editing_polyA-plus/testUniq_filtered.vcf -genotypeMergeOptions UNIQUIFY -R genome/GRCh38.fa -o testUniq_vars_polyAminus_polyAplus.vcf &
# combine ctrl uniq vars
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 -V:polyA_minus editing_polyA-minus/controlUniq_filtered.vcf -V:polyA_plus editing_polyA-plus/controlUniq_filtered.vcf -genotypeMergeOptions UNIQUIFY -R genome/GRCh38.fa -o ctrlUniq_vars_polyAminus_polyAplus.vcf &

# GET variants which are known editng events
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 -V:polyA_minus editing_polyA-minus/testUniq_filtered_known_edits.vcf -V:polyA_plus editing_polyA-plus/testUniq_filtered_known_edits.vcf -genotypeMergeOptions UNIQUIFY -R genome/GRCh38.fa -o testUniq_vars_polyAminus_polyAplus.known.vcf &
# combine ctrl uniq vars
java -jar /package/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineVariants -nt 20 -V:polyA_minus editing_polyA-minus/controlUniq_filtered_known_edits.vcf -V:polyA_plus editing_polyA-plus/controlUniq_filtered_known_edits.vcf -genotypeMergeOptions UNIQUIFY -R genome/GRCh38.fa -o ctrlUniq_vars_polyAminus_polyAplus.known.vcf &
