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

COUNT=1

for i in $(cat $BAM_LIST); do

	grep -v $i $BAM_LIST > sv_analysis/BAM_LIST_"$COUNT".txt

	echo $COUNT

	# Take BAMs and generate coverage data for each exon/bin
	Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$VERSION"/sv_calling/ReadInBams.R \
		--bams sv_analysis/BAM_LIST_"$COUNT".txt \
		--bed $EXON_BED \
		--fasta $FASTA \
		--out sv_analysis/"$WORKSHEET"_"$COUNT"

	# Calculate QC stats and produce report showing failed samples and regions
	Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$VERSION"/sv_calling/IdentifyFailures.R \
		--Rdata sv_analysis/"$WORKSHEET"_"$COUNT".RData \
		--mincorr 0.98 \
		--mincov 160 \
		--out sv_analysis/"$WORKSHEET"_"$COUNT"

	# Call CNVs using coverage data - make plots in plots directory.
	Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$VERSION"/sv_calling/makeCNVcalls.R \
		--transProb 0.05 \
		--Rdata sv_analysis/"$WORKSHEET"_"$COUNT".RData \
		--out sv_analysis/"$WORKSHEET"_"$COUNT" -plot All --plotFolder sv_analysis/plots/

	# Remove files no longer needed
	rm sv_analysis/"$WORKSHEET"_"$COUNT".RData
	rm sv_analysis/BAM_LIST_"$COUNT".txt

	COUNT=$(($COUNT + 1 ))



done