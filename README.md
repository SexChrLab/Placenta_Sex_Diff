# Placenta_Sex_Diff
Sex differences in full-term, >= 36 weeks,  uncomplicated human placentas

## Sex differences in placenta expression 
#### Step 1: Process RNA sequence data
- FastQC, MultiQC, and bamstats for visualizing quality. bbduk to trim reads and remove adaptors. 
- Reads aligned to a GRCh38 sex chromosomome complement reference genome. See https://github.com/SexChrLab/XY_RNAseq for more details on aligning to a sex chromosomome complement reference genome
- Use the snakemake file `rna.snakefile` that includes the following steps:
 - trimming, alignmnet, obtaining counts, and generating bam stats. 

#### Step 2: Run Limma/Voom 
- Use the r file `R_limmaVooom_placentaSexDiff.r` that includes the following steps:
  - You will need a counts table generated from step 1 FeatureCounts and a pheno type file
  - run limma/voom
