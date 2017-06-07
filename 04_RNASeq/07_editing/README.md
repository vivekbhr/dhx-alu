
# Variant calling and editing detection for DHX9 KD over control from RNA-Seq data

- Vivek Bhardwaj, Jan 2017

**For the workflow, have a look at the image in : Extended Data Figure 10c**

### Follow these steps to reproduce the output

1. Run `Snakefile.varcall` to run variant calling on `polyA_plus` or `polyA_minus` samples.
2. Run `Snakefile.editing` to run RNA-editing-specific filtering on `polyA_plus` or `polyA_minus` samples.
3. Run `gatk_combinevar.sh <varcall_folder>` to combine vars from control/test replicates, and filter out SNPs only.
4. Run `combinePlusMinus_vars.sh` to combine `polyA_plus` and `polyA_minus` samples.
5. Run `Identify_editing.Rmd` to annotate variants, collect some stats and make plots.

**Note for external users :**
1. Please place the file GRCh38.fa inside the directory `genome`.
2. Place the fastq files downloaded from GEO into any fastq directory and provide
this under variable `FASTQ_DIR` before running `Snakefile.varcall`
3. Replace the path to commands `/package/<command>` with your own path.
