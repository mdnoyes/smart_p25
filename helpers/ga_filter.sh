#!/bin/bash

sample=$1

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/ga_richness_filter.py \
	-v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".passed_strand_bias.tsv \
	-f /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/reference/DSA_"$sample"_v1.0.0.fasta \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".no_homopolymer.tsv
