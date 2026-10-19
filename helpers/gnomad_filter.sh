#!/bin/bash

sample=$1

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/filter_gnomad.py \
	-v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/annotation_files/gnomad/"$sample".gnomad_annotated.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".no_gnomad.tsv
