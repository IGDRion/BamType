## Scripts
This repository contains the scripts used to make quality control of the annotation of new transcripts. BamType.R is used to extract caracteristics of reads/alignment of a LR-RNAseq sample. The other scripts are used to generate figures with the output of BamType. The README.md of this directory gives more information about all the scripts presented.

## Plots
The plots directory contains all the figures generated during the analyses.

# BamType

Tool for quality control of newly discovered transcript annotation in long-read (Oxford Nanopore Technology (ONT)) RNAseq analysis

## Description

BamType is a R script that allow the extraction of diverse features of the alignements of a BAM file obtained from long-read RNAseq analysis (features can be about the read (i.e read length, Quality Score), the transcript (i.e transcript length), or the alignment itself (i.e transcript coverage, read coverage, alignement score)). 
BamType distinguishes alignements that target newly discovered transcripts / known transcripts ; and mRNAs / lncRNAs. 
The purpose of this script is to facilitate the comparison of alignments targeting new transcripts and known transcripts, ensuring that the reads aligned to newly discovered transcripts, as well as the alignments and transcripts themselves, meet high-quality standards to confirm the quality of annotation of new transcripts.

## Prerequisites

- R version > 4.2

## Installation

Download BamType.R
Download all R packages needed for BamType : 
* *GenomicAlignments* (1.34.0)
* *dplyr* (1.1.0)
* *ggplot2* (3.4.2)
* *hexbin* (1.28.3)
* *gridExtra* (2.3)
* *rtracklayer* (1.58)

## Usage

1. Preparation of BAM files

Before using BamType, you need to create a new alignment file (BAM) for your sample:
- Create a multifasta that contains known AND new transcripts sequences : extract all exons from the extended annotation (example on how to do this : [**gtf_exon.sh**](https://github.com/IGDRion/BamType/blob/master/scripts/gtf_exon.sh)) and create the multifasta (you can use *gffreads*, example on how to do this : [**gtf2fasta.sh**](https://github.com/IGDRion/BamType/blob/master/scripts/gtf2fasta.sh)).
- Use an alignment program (*minimap2* (2.24) was used for the tests) to align fastq sequences of your sample to the multifasta just created.


2. Launch & Input
`
Rscript Bamtype.R \\\

--type \\\

--bam \\\

--gtf \\\

--seqsum \\\

--output \\\  
`

- *--type* type of sequencing: cdna or rna
- *--bam* BAM file of the sample to analyse
- *--gtf* ANNEXA output: extended annotation with known and novel transcript, with biotype of novel transcripts
- *--seqsum* sequencing summary file corresponding to the sample
- *--output* path for output files

3. Output

Two output directories are created : **csv/** and **plots/**  
- csv:
    - **all_data.csv*** : contains all alignements (including all secondary alignements) with their features
    - **primary_data.csv*** : contains only primary alignments
    - **subsample.csv*** : random subsample of 10% of the **primary_data.csv**
    - **sample_features.csv*** : contains median/mean features for the primary alignments of the entire sample (distinguishing the known/novel transcripts and mRNAs/lncRNAs)
    - **per_unique_transcript.csv*** : contains names of every transcript associated with their median coverage
- plots:
    - accuracy repartition
    - density of coverage depending of transcript length
    - coverage fraction distribution
    - repartition of the number of secondary alignements
    - transcript length distribution
    - number of read aligned on known vs novel transcripts
    - proportions of biotypes in the known vs novel transcripts

Scripts using the .csv files (**primary_data.csv** and **subsample.csv**) to generate figures are shown as example in the "script/" directory of this Gitlab repository.
To merge all the subsamples of the analysis and make comparisons across samples, **merge_sample_bamtype.R**(https://github.com/IGDRion/BamType/blob/master/scripts/merge_sample_bamtype.R) was used.

## Example

Rscript Bamtype.R --type="cdna" --bam="path/to/501Mel_1-3.bam" --gtf="path/to/extended_annotation.gtf" --seqsum="path/to/501Mel_1-3_seqsum.txt" --output="bamtype_results/"

## Troubleshooting

BamType can cause *out-of-memory* errors: be sure to allocate enough ressources when launching the script.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/IGDRion/BamType/blob/master/LICENSE)

## Acknowledgments

BamType was developed from [BamSlam](https://github.com/josiegleeson/BamSlam/tree/master), a script written by **Josie Gleeson** to analyse ONT sequencing alignment.