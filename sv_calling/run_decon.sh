#!/bin/bash
set -euo pipefail

EXON_BED=$1 # The Exon BED file
BAM_LIST=$2 # List of BAMs to process
FASTA=$3 # The FASTA used BAMs were aligned to
WORKSHEET=$4 # Worksheet ID - used for naming output
VERSION=$5 # The pipeline version

#Make directories for storing data
mkdir sv_analysis
mkdir sv_analysis/plots

# Take BAMs and generate coverage data for each exon/bin
Rscript /home/joseph/Documents/decon/GermlineEnrichment/sv_calling/sv_calling/ReadInBams.R \
	--bams $BAM_LIST \
	--bed $EXON_BED \
	--fasta $FASTA \
	--out sv_analysis/"$WORKSHEET"

# Calculate QC stats and produce report showing failed samples and regions
Rscript /home/joseph/Documents/decon/GermlineEnrichment/sv_calling/sv_calling/IdentifyFailures.R \
	--Rdata sv_analysis/"$WORKSHEET".RData \
	--mincorr 0.98 \
	--mincov 160 \
	--out sv_analysis/"$WORKSHEET"

# Call CNVs using coverage data - make plots in plots directory.
Rscript /home/joseph/Documents/decon/GermlineEnrichment/sv_calling/sv_calling/makeCNVcalls.R \
	--transProb 0.05 \
	--Rdata sv_analysis/"$WORKSHEET".RData \
	--out sv_analysis/"$WORKSHEET" -plot All --plotFolder sv_analysis/plots/

# Remove files no longer needed
rm sv_analysis/"$WORKSHEET".RData

