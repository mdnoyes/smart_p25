#!/bin/bash

sample=$1


python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/inter_snp_distance.py \
	-i /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".in_good_regions.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".intersnp_annotated.tsv

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/get_clusters.py  \
	-i /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".intersnp_annotated.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".clusters.tsv

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/remove_clusters.py \
	-v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".intersnp_annotated.tsv \
	-c /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".clusters.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".no_clustered_variants.tsv

rm /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/"$sample".intersnp_annotated.tsv  

