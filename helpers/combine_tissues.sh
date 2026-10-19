#!/bin/bash

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/combine_tissues.py \
        -v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT001.hprc_filtered.tsv \
	-t adrenal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.adrenal.tsv \
	-t aorta /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.aorta.tsv \
	-t blood /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.blood.tsv \
	-t brain-cerebellum /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.brain-cerebellum.tsv \
	-t brain-frontal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.brain-frontal.tsv \
	-t brain-hippocampus-L /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.brain-hippocampus-L.tsv \
	-t brain-hippocampus-R /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.brain-hippocampus-R.tsv \
	-t brain-temporal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.brain-temporal.tsv \
	-t colon-asc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.colon-asc.tsv \
	-t colon-desc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.colon-desc.tsv \
	-t esophagus /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.esophagus.tsv \
	-t fibroblast /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.fibroblast.tsv \
	-t heart /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.heart.tsv \
	-t liver /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.liver.tsv \
	-t lung /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.lung.tsv \
	-t muscle /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.muscle.tsv \
	-t skin-abd /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.skin-abd.tsv \
	-t skin-calf /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.skin-calf.tsv \
	-t testis /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT001.testis.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT001.combined_tissues.tsv


python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/combine_tissues.py \
	-v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT004.hprc_filtered.tsv \
	-t aorta /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.aorta.tsv \
	-t blood /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.blood.tsv \
	-t buccal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.buccal.tsv \
	-t colon-asc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.colon-asc.tsv \
	-t colon-desc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.colon-desc.tsv \
	-t esophagus /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.esophagus.tsv \
	-t fibroblast /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.fibroblast.tsv \
	-t heart /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.heart.tsv \
	-t liver /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.liver.tsv \
	-t lung /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.lung.tsv \
	-t muscle /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.muscle.tsv \
	-t skin-abd /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.skin-abd.tsv \
	-t skin-calf /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.skin-calf.tsv \
	-t testis /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT004.testis.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT004.combined_tissues.tsv

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/combine_tissues.py \
        -v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT005.hprc_filtered.tsv \
	-t adrenal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.adrenal.tsv \
	-t aorta /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.aorta.tsv \
	-t blood /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.blood.tsv \
	-t brain-cerebellum /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.brain-cerebellum.tsv \
	-t brain-frontal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.brain-frontal.tsv \
	-t brain-hippocampus-L /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.brain-hippocampus-L.tsv \
	-t brain-hippocampus-R /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.brain-hippocampus-R.tsv \
	-t brain-temporal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.brain-temporal.tsv \
	-t buccal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.buccal.tsv \
	-t colon-asc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.colon-asc.tsv \
	-t colon-desc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.colon-desc.tsv \
	-t esophagus /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.esophagus.tsv \
	-t fibroblast /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.fibroblast.tsv \
	-t heart /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.heart.tsv \
	-t liver /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.liver.tsv \
	-t lung /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.lung.tsv \
	-t muscle /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.muscle.tsv \
	-t skin-abd /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.skin-abd.tsv \
	-t skin-calf /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.skin-calf.tsv \
	-t testis /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT005.testis.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT005.combined_tissues.tsv

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/combine_tissues.py \
        -v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT012.hprc_filtered.tsv \
	-t adrenal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.adrenal.tsv \
	-t aorta /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.aorta.tsv \
	-t blood /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.blood.tsv \
	-t colon-asc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.colon-asc.tsv \
	-t esophagus /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.esophagus.tsv \
	-t fibroblast /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.fibroblast.tsv \
	-t heart /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.heart.tsv \
	-t liver /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.liver.tsv \
	-t lung /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.lung.tsv \
	-t muscle /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.muscle.tsv \
	-t skin-abd /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.skin-abd.tsv \
	-t skin-calf /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.skin-calf.tsv \
	-t testis /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT012.testis.tsv \
	-o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT012.combined_tissues.tsv

python /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/scripts/combine_tissues.py \
        -v /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT024.hprc_filtered.tsv \
	-t adrenal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.adrenal.tsv \
	-t blood /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.blood.tsv \
	-t brain-cerebellum /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.brain-cerebellum.tsv \
	-t brain-frontal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.brain-frontal.tsv \
	-t brain-hippocampus-L /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.brain-hippocampus-L.tsv \
	-t brain-hippocampus-R /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.brain-hippocampus-R.tsv \
	-t brain-temporal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.brain-temporal.tsv \
	-t buccal /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.buccal.tsv \
	-t colon-asc /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.colon-asc.tsv \
	-t colon-desc.tsvv /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.colon-desc.tsvv \
	-t esophagus /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.esophagus.tsv \
	-t fibroblast /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.fibroblast.tsv \
	-t heart /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.heart.tsv \
	-t liver /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.liver.tsv \
	-t lung /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.lung.tsv \
	-t muscle /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.muscle.tsv \
	-t skin-abd /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.skin-abd.tsv \
	-t skin-calf /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.skin-calf.tsv \
	-t testis /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/tissue_data/SMHT024.testis.tsv \
        -o /mmfs1/gscratch/smaht/mdnoyes/P25_first_five/filtered_variants/SMHT024.combined_tissues.tsv
