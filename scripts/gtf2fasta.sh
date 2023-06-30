#!/bin/sh

# sbatch spec
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

GENOME="/groups/dog/data/hg38_GRCh38/sequence/softmasked/Homo_sapiens.GRCh38.dna_sm.primary_assembly_withChr.fa"
INFILE_GTF="/scratch/vlebars/extended_annotations_exon.filter.gtf"
OUTFILE_FA="/scratch/vlebars/new_transcriptome.filter.fa"

gffread -w $OUTFILE_FA -g $GENOME $INFILE_GTF 