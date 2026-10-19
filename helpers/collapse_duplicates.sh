#!/bin/bash

sample=$1


python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/collapse_duplicates.py \
	-v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".no_clustered_variants.tsv  \
	-b /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/annotation_files/move_haplotypes/surjected_outputs/"$sample"/"$sample"-SurjOnto-"$sample"_1_and_"$sample"_2_MERGED_sorted_extracted.bed  \
	-f /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/reference/DSA_"$sample"_v1.0.0.fasta \
	-c /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/cram_processing/retag/hifi/"$sample".fibroblast.retagged.cram \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".collapsed_duplicates.tsv
