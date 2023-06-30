# Scripts

This directory countains scripts used for quality control of newly discovered transcripts annotation. 

## BamType.R
BamType.R is used to extract features of alignements depending of the nature of the transcript (known/new and mRNA/lncRNA).

## gtf_exon.sh & gtf2fasta.sh
Scripts used before BamType to create the multifasta to which the reads are aligned.

## merge_sample_bamtype.R
Merges all sample.csv files of the target directory into one .csv file (useful for cross sample comparisons).

## fig_across_samples.ipynb
Notebook to create various plots comparing alignment features over several samples. It take the output of **merge_sample_bamtype.R** in input.

## read_count_across_samples.ipynb
Creates a plot representing the number of read (mapped and unmapped) of all the samples of the analysis. It takes [**count_output.csv**](https://github.com/IGDRion/BamType/blob/master/input/count_output.csv) in input, which is obtain by keeping unmapped read in the alignment (*minimap2: --sam-hit-only=FALSE*) step and counting both unmapped reads and mapped reads with samtools (see **read_count_launch.sh** for an example on how to do that).

## sankey_diagram.ipynb
Creates a sankey diagram representing the repartition of reads ont the different nature of transcripts (known/new and mRNA/lncRNA/others) within the replicates of a sample.  
It uses output of **merge_csv.py** as input. **merge_csv.py** need the path to the *sample_features* directory (containing features means and medians of a sample) and the output of **count_output.csv** as input.