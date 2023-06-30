#!/bin/sh

# sbatch spec
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

INFILE="/groups/dog/nanopore/lncrna_resist_cgo/secondary/1_annexa/results/final/extended_annotations.filter.gtf"
OUTFILE="/scratch/vlebars/extended_annotations_exon.filter.gtf"

awk '$3=="exon"' $INFILE > $OUTFILE