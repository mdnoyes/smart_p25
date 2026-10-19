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


def get_flagged_intervals(test_bed: BedTool, bed_file: str) -> Dict[Tuple, str]:
	"""
	Find the overlap between the input variants and the bed regions.
	:param test_bed: BedTool object of intervals for your de novo variants.
	:param bed_file: Path to bed file of regions of interest.
	:return: List of variant intervals that overlap with the bed regions.
	"""
	ref_bed = BedTool(bed_file)
	intersected_variants = test_bed.intersect(ref_bed, wo=True)
	
	flagged = {}
	for interval in intersected_variants:
		annotation = interval.fields[6]
		flagged[(interval.chrom, interval.start, interval.end)] = annotation

	return flagged


def main(variant_file, flagged_bed, output_file):
	header, variant_dict = read_variant_file(variant_file)
	out_header = header + ['NucFlag']
	variant_bed = BedTool(list(variant_dict.keys()))
	
	flagged_intervals = get_flagged_intervals(variant_bed, flagged_bed)


	with open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)

		for interval, var_info in tqdm(variant_dict.items(), total = len(variant_dict.items())):
			annotation = flagged_intervals.get(interval)
			if not annotation:
				annotation = 'NA'

			if annotation in ('collapse', 'misjoin', 'low_quality', 'false_dup', 'het_or_mismap'):
				continue
			
			row_data = [var_info[col] for col in header] + [annotation]
			print('\t'.join(row_data), file = outfile)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Variant file")
	ap.add_argument("-b", "--flagged_bed", required=True, help="bed file of flagged variants")
	ap.add_argument("-o", "--output", required=True, help="output file")


	args = ap.parse_args()

	main(args.input, args.flagged_bed, args.output)
