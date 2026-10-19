#!/bin/bash

sample=$1
dir=/mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants

awk '{if($18!="False") {print $0}}' $dir/"$sample".combined_tissues.tsv | awk '{if($19!="False") {print $0}}' > $dir/"$sample".passed_binomial_test.tsv
awk '{if($20!="FAIL") {print $0}}' $dir/"$sample".passed_binomial_test.tsv > $dir/"$sample".passed_strand_bias.tsv

