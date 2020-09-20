#!/usr/bin/sh

## Annotate circRNA from STAR mapping outputs (CIRCexplorer 2.1 )

## New feature in version 2.1 is that it can directly find circRNAs coming from splicing events, 
## if I provide the poly-A plus Tophat output.

## output saved in circRNA_annotation directory
mkdir circRNA_annotation

#hg38
hg=/data/repository/organisms/GRCh38_ensembl/genome_fasta/genome.fa
hCircRef=hg38_ens.circRef.txt 
# human
for file in $(ls circRNA_mapping | grep star_)
do base=$(basename ${file} | sed 's/star\_//')
# parse STAR output
CIRCexplorer2 parse -t STAR circRNA_mapping/star_${base}/${base}Chimeric.out.junction -o circRNA_annotation/${base}
# annotate circRNA
CIRCexplorer2 annotate -g $hg -r $hCircRef circRNA_annotation/${base}
# denovo find splice junctions related to circRNAs
CIRCexplorer2 denovo --as --rpkm -g $hg -r $hCircRef -a ../02_autoMapping/analysis_results/polyA_plus/TopHat2/${base} \
--tophat-dir ../02_autoMapping/analysis_results/polyA_minus/TopHat2/${base} circRNA_annotation/${base}
done

