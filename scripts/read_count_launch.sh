#! /bin/bash

# sbatch spec
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G

. /local/env/envsamtools-1.15.sh

# Get the sample name
sample_name=$(basename $1 .bam)


# ---------------------------------------------------------------------------- #
#                       Counting mapped / unmapped reads                       #
# ---------------------------------------------------------------------------- #

# Run samtools view -c on the BAM file, redirect output to a temp file, for primary aligned mapped reads
# Cannot just use samtools view -c because of supplementary alignments
samtools view -F 260 $1 | awk '{ print $1 }' | sort -n | uniq | wc -l > "/scratch/vlebars/tmp-$sample_name-mapped.txt"

# Run samtools view -c on the BAM file, redirect output to a temp file, for unmapped reads
samtools view -c -f 4 $1 > "/scratch/vlebars/tmp-$sample_name-unmapped.txt"

# Append the count and sample name to the csv file
echo "$sample_name-mapped,"$(cat /scratch/vlebars/tmp-$sample_name-mapped.txt) >> /groups/dog/stage/victor/count_output.csv
echo "$sample_name-unmapped,"$(cat /scratch/vlebars/tmp-$sample_name-unmapped.txt) >> /groups/dog/stage/victor/count_output.csv

# Clean up the temp files
rm /scratch/vlebars/tmp-$sample_name*

