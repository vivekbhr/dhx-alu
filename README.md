# DHX9-Alu 

Scripts and pipelines corresponding to the article : **DHX9 suppresses RNA processing defects originating from the Alu invasion of human genome** (Aktas et. al., Nature 2017).

This directory contains scripts related to : 

1. Repeat enrichment after mapping to a library of repeat sequences.
2. Repeat clustering and annotation from fasta files.
3. RNA-Seq analysis of DHX9 knockdown.


## Repeat enrichment after mapping to a library of repeat sequences

We followed the strategy described in [Dey. et. al., 2010](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-6-r69) to calculate enrichment of repeat classes from DHX9 uvCLAP data, for Human (hg38), Mouse (mm10) and Fly (dm6) genome. In short :

+ Repeat library was created using sequences of cannonical repeats and repeat instances. [Script](./Dey_et.al_pipeline/make_repeat_assembly.pl)

+ uvCLAP reads were mapped to the library using BWA, and counts unique to each repeat class was obtained.

+ Mapping was done again using bowtie2, and counts from BWA were normalized by library sizes obtained 
  from the bowtie2 mapping. [Script](./make_bowtieCountsStats.sh)

+ Maximum Likelihood Estimates (MLE) was calculated for each repeat class, in each sample. [Scripts](./Rscripts_Dhx-Alu/estimate_repenr.R) 

+ All samples were clustered using MLE estimates and plotted. [Scripts](./Rscripts_Dhx-Alu/batchEnrichment.R)


## Repeat clustering and annotation from fasta files.

We followed a graph-based clustering strategy originaly described in [Novák et al., 2010](http://www.biomedcentral.com/1471-2105/11/378), and implemented as a pipeline in [Novák et al., 2013](https://academic.oup.com/bioinformatics/article/29/6/792/184407/RepeatExplorer-a-Galaxy-based-web-server-for). In short : 

+ We sampled 100,000 fastq reads from the uvCLAP experiments and converted to fasta.

+ Ran the **repeatexplorer** pipeline obtained from [bitbucket](https://bitbucket.org/repeatexplorer/repeatexplorer), using *seqclust_cmd.py*.

+ Alu/B1 SINE IDs were extracted from repeatmasker annotated clusters [Script]() and [SeqGrapheR](https://github.com/vivekbhr/SeqGrapheR) was used to visualize the clusters.


## RNA-Seq analysis of DHX9 knockdown

scripts coming soon..