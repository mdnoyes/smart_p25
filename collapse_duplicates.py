from __future__ import annotations
from pysam import AlignmentFile, AlignedSegment, FastaFile
from typing import Dict, List, Tuple, Optional, TextIO, Set
from tqdm import tqdm
from collections import defaultdict
import argparse


class Variant:
    variant_id: str
    variant_info: List[str]
    chrom: str
    pos: int
    homolog_chrom: Optional[str]
    homolog_pos: Optional[int]
    homolog_ref_allele: Optional[str]
    alternate_dict: Dict[str, set[str]]
    alternate_alignments: set[str]
    alternate_calls: set[str]
    is_duplicate: bool

    def __init__(self, variant_info: Dict[str, str]) -> None:
        self.variant_info = variant_info
        self.variant_id = variant_info['id']
        self.chrom = variant_info['chrom']
        self.pos = int(variant_info['pos'])

        self.homolog_chrom = None
        self.homolog_pos = None
        self.homolog_ref_allele = None

        self.alternate_dict = {}
        self.alternate_alignments = set()

        self.is_duplicate = False
        self.alternate_calls = set()


    def add_homolog(self, chrom: str, pos: int, fasta: FastaFile):
        if self.chrom != chrom or self.pos != pos:
            self.homolog_chrom = chrom
            self.homolog_pos = pos
            self.homolog_ref_allele = fasta.fetch(chrom, pos - 1, pos).upper()


    def add_alternates(self, read: AlignedSegment) -> None:
        alt_alignments = read.get_tag('XA') if read.has_tag('XA') else None
        if not alt_alignments:
            return

        alt_alignments = alt_alignments.split(';')

        for alignment in alt_alignments[:-1]:
            chrom, pos = alignment.split(',')[:2]
            if pos.startswith('-'):
                pos = pos[1:]

            if chrom not in self.alternate_dict:
                self.alternate_dict[chrom] = set()

            self.alternate_dict[chrom].add(int(pos))


    def collect_alternates(self, threshold: int = 151):  
        """
        Define any potential paralogous sequences based on alternate alignments tag in cram. 
        """
        if not self.alternate_dict:
            self.alternate_alignments = None
            return 

        self.alternate_alignments = set()
        primary_alignments = [(self.chrom, self.pos)]
        if self.homolog_chrom and self.homolog_pos:
            primary_alignments.append((self.homolog_chrom, self.homolog_pos))
    
        for chrom, positions in self.alternate_dict.items():
            position_list = sorted(positions)
            start = position_list[0]
            current_end = start

            for position in position_list[1:]:
                if position - start <= threshold:
                    current_end = position
                else:
                    if not any(chrom == tc and start <= tp <= current_end for tc, tp in primary_alignments):
                        self.alternate_alignments.add((chrom, start, current_end))
                    start = position
                    current_end = position

            if not any(chrom == tc and start <= tp <= current_end for tc, tp in primary_alignments):
                self.alternate_alignments.add((chrom, start, current_end))


    def check_alternates(self, variant_index: VariantIndex) -> List[str]:
        if not self.alternate_alignments:
            return []
            
        variants = variant_index.variants
        alternate_calls = set()
        for chrom, start, end in self.alternate_alignments:
            alternate_calls.update(variant_index.find_variants_in_interval(
                chrom = chrom, 
                start = start, 
                end = end, 
                exclude_id = self.variant_id
            ))
           
        return alternate_calls

    def print_variant(self, outfile: TextIO) -> None:
        if self.alternate_calls:
            collapsed_calls = ';'.join(self.alternate_calls)
        else:
            collapsed_calls = 'NA'

        if self.homolog_chrom and self.homolog_pos:
            homolog_chrom = self.homolog_chrom
            homolog_pos = str(self.homolog_pos)
            homolog_ref_allele = self.homolog_ref_allele
        else:
            homolog_chrom = homolog_pos = homolog_ref_allele = 'NA'

        output = [self.chrom, str(self.pos), self.variant_id, self.variant_info['ref'], self.variant_info['alt'],
            homolog_chrom, homolog_pos, homolog_ref_allele, collapsed_calls]

        print('\t'.join(output), file = outfile)

    
class VariantIndex:
    variants: Dict[str, Variant]
    coordinate_lookup: Dict[Tuple(str, int), List[str]]
    pos_index: Dict[str, Dict[int, List[str]]]
    
    def __init__(self, variant_file: str) -> None:
        self.variants = {}

        with open(variant_file, 'r') as infile:
            header = infile.readline().rstrip().split('\t')

            for line in infile:
                variant = {k:v for k,v in zip(header, line.rstrip().split('\t'))}
                self.variants[variant['id']] = Variant(variant) 


    def import_bed(self, bed_path: str, fasta: FastaFile) -> None:
        coordinate_lookup = defaultdict(list)
        self.pos_index = defaultdict(lambda: defaultdict(set))
#        pos_index = defaultdict(lambda: defaultdict(list))

        with open(bed_path, 'r') as infile:
            header = infile.readline()
            for line in infile:
                variant = line.rstrip().split('\t')
                variant_id = variant[3]

                chrom = variant[0]
                pos = int(variant[2])

                if variant_id in self.variants:
                    self.variants[variant_id].add_homolog(chrom, pos, fasta)
                    self.pos_index[chrom][pos].add(variant_id)
                    coordinate_lookup[(chrom, pos)].append(variant_id)
#                if (chrom, pos) not in coordinate_lookup:
#                    coordinate_lookup[(chrom, pos)] = []
#                coordinate_lookup[(chrom, pos)].append(variant_id)
#                pos_index[chrom][pos].append(variant_id)
        
#        self.pos_index = pos_index
        self.coordinate_lookup = coordinate_lookup


    def import_alternates(self, file_path: str) -> None:
        with open(file_path, 'r') as infile:
            header = infile.readline().rstrip().split('\t')
            for line in infile:
                variant = {k:v for k,v in zip(header, line.rstrip().split('\t'))}
                variant_id = variant['id']
                # if variant_id.startswith('haplotype1'):
                #     variant_id = f'SMHT004#1#{variant_id}'
                # elif variant_id.startswith('haplotype2'):
                #     variant_id = f'SMHT004#2#{variant_id}'

                if variant_id in self.variants:
                    if variant['alternate_alignments'] == 'NA':
                        self.variants[variant_id].alternate_alignments = None
                        continue
                    
                    alignment_list =  variant['alternate_alignments'].split(';')
                    alternate_alignments = set()

                    for alignment in alignment_list:
                        chrom, interval = alignment.split(':')
                        # if chrom.startswith('haplotype1'):
                        #     chrom = f'SMHT004#1#{chrom}'
                        # elif chrom.startswith('haplotype2'):
                        #     chrom = f'SMHT004#2#{chrom}'
                        start, end = [int(x) for x in interval.split('-')]
                        alternate_alignments.add((chrom, start, end))

                    self.variants[variant_id].alternate_alignments = alternate_alignments


    def query_cram(self, variant: Variant, cram: AlignmentFile, chromosome: str, position: int) -> None:
        for read in cram.fetch(chromosome, position - 1, position):
            variant.add_alternates(read)

    def sort_variants(self, cram: AlignmentFile) -> List[Tuple[int, int, str, int, Variant]]:
        chrom_order = {name: i for i, name in enumerate(cram.references)}

        query_targets = []
        for variant in self.variants.values():
            if variant.chrom in chrom_order:
                query_targets.append((chrom_order[variant.chrom], variant.pos, variant.chrom, variant.pos, variant))
            
            if variant.homolog_chrom and variant.homolog_pos and variant.homolog_chrom in chrom_order:
                query_targets.append((chrom_order[variant.homolog_chrom], variant.homolog_pos, variant.homolog_chrom, variant.homolog_pos, variant))

        query_targets.sort(key=lambda x: (x[0], x[1]))

        return query_targets


    def import_cram(self, cram_path: str) -> None:
        cram = AlignmentFile(cram_path, "rc")

        sorted_variants = self.sort_variants(cram)

        for _ , _, chrom, pos, variant in tqdm(sorted_variants):
            self.query_cram(variant, cram, chrom, pos)

        for variant in self.variants.values():
            variant.collect_alternates()


    def check_homologous_coordinates(self) -> None:    
        coordinate_sets = defaultdict(set)
        for group in self.coordinate_lookup.values():
            for variant_id in group:
                coordinate_sets[variant_id].update(x for x in group if x != variant_id)         

        seen = set()
        for variant_id in list(self.variants.keys()):
            if variant_id in seen or variant_id not in coordinate_sets:
                continue

            # condense all matches into one big group 
            # for example, if both h1 variants A and B lift to the site of variant C on h2
            # then A, B, and C should all be grouped
            # this is important because not all lifts are reciprocal, so we have to essentially chain them
            matching_coordinates = set()
            match = [variant_id]
            while match:
                curr = match.pop()
                if curr not in matching_coordinates:
                    matching_coordinates.add(curr)
                    match.extend(list(coordinate_sets[curr] - matching_coordinates))

            seen.update(matching_coordinates)

            if len(matching_coordinates) > 1:
                sorted_matches = sorted(matching_coordinates)
                keep = sorted_matches[0]
                remove = sorted_matches[1:]

                self.variants[keep].alternate_calls.update(remove)

                for variant_id in remove:
                    if variant_id in self.variants:
                        self.variants[variant_id].is_duplicate = True


    def find_variants_in_interval(self, chrom: str, start: int, end: int, exclude_id: str) -> Set[str]:
        """
        Fast lookup: Returns all variant IDs on `chrom` in range [start, end].
        """
        matching_ids = set()
        
        if chrom not in self.pos_index:
            return matching_ids

        # Only inspect integer positions present in the index on this chromosome
        chrom_positions = self.pos_index[chrom]            
        
        # Iterate over positions within [start, end]
        # (For genomic windows ~150bp, this loops ~150 steps maximum)
        for p in range(start, end + 1):
            if p in chrom_positions:
                for variant_id in chrom_positions[p]:
                    if variant_id != exclude_id:
                        matching_ids.add(variant_id)

        return matching_ids


    def check_alternate_alignments(self) -> None:
        seen = set()

        for variant in self.variants.values():
            if variant.variant_id in seen:
                continue
            
            seen.add(variant.variant_id)
            
            alternate_calls = set(variant.check_alternates(self))
            if alternate_calls:
                seen.update(alternate_calls) 
                alternate_calls.add(variant.variant_id)
                sorted_matches = sorted(alternate_calls)
                
                keep = None
                remove = []

                for variant_id in sorted_matches:
                    if not self.variants[variant_id].is_duplicate and not keep:
                        keep = variant_id
                        continue
                    remove.append(variant_id)
                        
                if not keep:
                    try:
                        keep = next(k for k,v in self.variants.items() if remove[0] in v.alternate_calls)
                    except StopIteration:
                        continue

                self.variants[keep].alternate_calls.update(remove) 

                for variant_id in remove:
                    self.variants[variant_id].is_duplicate = True
                    if self.variants[variant_id].alternate_calls:
                        self.variants[keep].alternate_calls.update(self.variants[variant_id].alternate_calls)


    def print_variants(self, output_path) -> None:
        output_header = ('chrom', 'pos', 'id', 'ref', 'alt', 
            'homolog_chrom', 'homolog_pos', 'homolog_ref_allele', 'collapsed_calls')

        with open(output_path, 'w') as outfile:
            print('\t'.join(output_header), file = outfile)
            for variant in self.variants.values():
                if not variant.is_duplicate:
                    variant.print_variant(outfile)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=True, help="variant file")
    ap.add_argument("-b", "--bed", required=True, help="bed file with homolog positions")
    ap.add_argument("-f", "--fasta", required=True, help="fasta file")
    ap.add_argument("-a", "--alternates", required=False, help="debugging alternate file")
    ap.add_argument("-c", "--cram", required=True, help="fibroblast illumina cram")
    ap.add_argument("-o", "--output", required=True, help="Output")

    args = ap.parse_args()  

    print('importing fasta...')
    fasta = FastaFile(args.fasta)

    print('importing variants...')
    variants = VariantIndex(args.variants)

    # for debugging
    # print('importing alternates...')
    # variants.import_alternates(args.alternates)

    print('importing bed...')
    variants.import_bed(args.bed, fasta)

    print('importing cram...')
    variants.import_cram(args.cram)

    print('checking homologs...')
    variants.check_homologous_coordinates()

    print('checking alternates...')
    variants.check_alternate_alignments()

    print('pringing variants...')
    variants.print_variants(args.output)

    print('done')

   
