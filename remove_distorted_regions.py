from pybedtools import BedTool
from typing import Tuple, List, Dict, Iterable, Optional
from tqdm import tqdm
import argparse


def read_variant_file(variant_file: str):
	"""
	Read in variants, save them to a dictionary with bed intervals.
	:param variant_file: Path to file of de novo mutations.
	:return: Header of variant file as list of strings.
	:return: Dictionary of variant info mapped to bed interval.
	"""
	variants = {}
	with open(variant_file, 'r') as var_file:
		header = var_file.readline().rstrip().split('\t')
		for line in var_file:
			var_data = line.rstrip().split('\t')
			pos = int(var_data[header.index('pos')])
			bed_info = (var_data[header.index('chrom')], pos - 1, pos)
			variants[bed_info] = {k:v for k,v in zip(header, var_data)}

	return header, variants


def get_flagged_intervals(test_bed: BedTool, bed_file: str) -> List[Tuple]:
	"""
	Find the overlap between the input variants and the bed regions.
	:param test_bed: BedTool object of intervals for your de novo variants.
	:param bed_file: Path to bed file of regions of interest.
	:return: List of variant intervals that overlap with the bed regions.
	"""
	ref_bed = BedTool(bed_file)
	intersected_variants = test_bed.intersect(ref_bed, wo=True)
	
	distorted = []
	for interval in intersected_variants:
		distorted.append((interval.chrom, interval.start, interval.end))

	return distorted


def main(variant_file, h1_bed, h2_bed, output_file):
	header, variant_dict = read_variant_file(variant_file)
	variant_bed = BedTool(list(variant_dict.keys()))
	
	h1_distorted_intervals = get_flagged_intervals(variant_bed, h1_bed)
	h2_distorted_intervals = get_flagged_intervals(variant_bed, h2_bed)

	distorted_intervals = h1_distorted_intervals + h2_distorted_intervals

	with open(output_file, 'w') as outfile:
		print('\t'.join(header), file = outfile)

		for interval, var_info in tqdm(variant_dict.items(), total = len(variant_dict.items())):
			if interval in distorted_intervals:
				continue
		
			row_data = [var_info[col] for col in header]
			print('\t'.join(row_data), file = outfile)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Variant file")
	ap.add_argument("-b1", "--h1_distorted_bed", required=True, help="distortion map bed")
	ap.add_argument("-b2", "--h2_distorted_bed", required=True, help="distortion map bed")
	ap.add_argument("-o", "--output", required=True, help="output file")


	args = ap.parse_args()

	main(args.input, args.h1_distorted_bed, args.h2_distorted_bed, args.output)
