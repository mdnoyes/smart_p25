#!/bin/bash

sample=$1

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/remove_distorted_regions.py \
	-i /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".collapsed_duplicates.tsv  \
	-b1 /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/annotation_files/distortion_maps/beds_contigs/"$sample"_hap1.donor_strong.contig.bed \
	-b2 /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/annotation_files/distortion_maps/beds_contigs/"$sample"_hap2.donor_strong.contig.bed \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".no_distortions.tsv

