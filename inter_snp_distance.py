import argparse
from typing import Tuple, List, Dict, Iterable, Optional

class Neighbor:
	chrom: str
	pos: int
	row: List[str]
	previous_neighbor: 'Neighbor'
	next_neighbor: 'Neighbor'
	close_neighbors: List['Neighbor']

	def __init__(
		self,
		row: List[str],
		previous_neighbor: 'Neighbor' = None,
		next_neighbor: 'Neighbor' = None,
	):
		self.row = row
		self.chrom = self.row[0]
		self.pos = int(self.row[1])
		self.previous_neighbor = previous_neighbor
		self.next_neighbor = next_neighbor

	def dist_to(self, other: Optional['Neightbor']) -> int:
		if other is not None:
			if self.chrom == other.chrom:
				return abs(self.pos - other.pos)

		return -1

	def dist_to_prev(self) -> int:
		return self.dist_to(self.previous_neighbor)

	def dist_to_next(self) -> int:
		return self.dist_to(self.next_neighbor)

	def close_neighbors(self) -> Tuple[int, 'Neighbor']:
		"""
		:returns: (number of neighbors in cluster, last neighbor in cluster)
		"""
		current = self
		num_neighbors = 0

		while current is not None and current.dist_to_prev() in range(0, 1001):
			current = current.previous_neighbor
			num_neighbors += 1

		current = self
		while current is not None and current.dist_to_next() in range(0, 1001):
			current = current.next_neighbor
			num_neighbors += 1

		return num_neighbors, current

def main(input_file, output_file):
	with open(input_file, 'r') as snps, open(output_file, 'w') as output:
		header = snps.readline().rstrip().split('\t')
		header.append('distance_to_closest_variant')
		print('\t'.join(header), file = output)

		last = None
		for line in snps:
			row = line.rstrip().split('\t')
			last = Neighbor(row, last)

		first = last
		while first.previous_neighbor is not None:
			tmp, first = first, first.previous_neighbor
			if first is not None:
				first.next_neighbor = tmp

		current = first
		while current is not None:
			num_neighbors, last_neighbor = current.close_neighbors()
			neighbor_distance = [x for x in [current.dist_to_prev(), current.dist_to_next()] if x != -1]
			try:
				distance = min(neighbor_distance)
			except ValueError:
				distance = -1
			print('\t'.join([*current.row, str(distance)]), file=output)
			current = current.next_neighbor


if __name__ == "__main__":
	# to run from the command line, use commented out lines
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="vcf of all family variant calls")
	ap.add_argument("-o", "--output", required=True, help="name of output file")
	args = ap.parse_args()

	main(args.input, args.output)
