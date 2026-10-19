import argparse
import itertools

def import_snvs(snv_file):
	mut_dict = {}

	with open(snv_file, 'r') as infile:
		header = infile.readline().rstrip().split('\t')
		header_idc = {k: v for v,k in enumerate(header)}

		for line in infile:
			variant = line.rstrip().split('\t')

			ref, alt = [variant[header_idc[x]] for x in ('ref', 'alt')]
			mut_dict[variant[header_idc['id']]] = f'{ref}>{alt}'

	return mut_dict


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Variant file")
	ap.add_argument("-o", "--output", required=True, help="output file")

	args = ap.parse_args()

	mutations = import_snvs(args.input)

	cluster = []
	cluster_id = 0

	with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
		print('\t'.join(['id1', 'id2', 'mutation1', 'mutation2', 'distance', 'pair_id']), file = outfile)

		header = infile.readline().rstrip().split('\t')
		header_idc = {k: v for v,k in enumerate(header)}

		for line in infile:
			variant = line.rstrip().split('\t')

			if 0 < int(variant[header_idc['distance_to_closest_variant']]) <= 500:
				cluster.append(variant[header_idc['id']])

			else:
				if len(cluster) == 0:
					continue

				for pair in itertools.combinations(cluster, 2):

					pos_1, pos_2 = [int(x.split('_')[1]) for x in pair]
					distance = abs(pos_2 - pos_1)
					if distance > 500:
						continue

					out_data = list(pair) + [mutations[x] for x in pair] + [str(distance), str(cluster_id)]
					print('\t'.join(out_data), file = outfile)
				
				cluster_id += 1
				cluster = []

	print('dont')
