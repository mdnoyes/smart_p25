from __future__ import annotations
from pysam import AlignmentFile, AlignedSegment, FastaFile
from typing import Dict, List, Tuple, Optional, TextIO
from tqdm import tqdm
from collections import Counter, defaultdict
import argparse
import gzip

class PresenceFilter:
    """
    Determines whether a variant has sufficient read support based on filtered alt allele count.
    The threshold for support is predefined for diffrent sequencing platforms. 

    Attributes:
        read_threshold: Minimum number of reads to establish variant presence, based on sequencing platform.
    """
    read_threshold: int

    def __init__(self, read_threshold, depth_threshold):
        self.read_threshold = read_threshold
        self.depth_threshold = depth_threshold

    @staticmethod
    def for_type(data_type: str) -> PresenceFilter:
        """
        Define read threshold for provided sequencing platform.

        Args:
            data_type: Sequencing platform ('hifi' or 'illumina'). Expect 'ont' in future versions.

        Returns:
            An instance of PresenceFilter with platform-specific threshold.

        Raises:
            KeyError: If unsupported data_type is specified. 
        """
        # HiFi tissues were sequenced to lower coverage, so 1 supporting read has to be sufficient.
        # Illumina has higher coverage, so we expected to see more supporting reads for a true somatic event.
        return {
            'hifi': PresenceFilter(read_threshold=1, depth_threshold=5),
            'illumina': PresenceFilter(read_threshold=3, depth_threshold=20),
        }[data_type]


    def determine_presence(self, variant: Variant) -> str:
        total_passing_alt, h1_count, h2_count, ambiguous_count, failed_count = [variant.get_read_count(x, 2) for x in ('passed', 'unique_to_h1', 'unique_to_h2', 'ambiguous', 'failed')]
        total_pass, total_fail, h1_dp, h2_dp = [variant.get_read_count(x, 0) for x in ('passed', 'failed', 'unique_to_h1', 'unique_to_h2')]

        # require at least 5 passing hifi or 20 passing illumina reads to make this call
        if total_pass >= self.depth_threshold:
            if variant.flag == 'germline_heterozygote':
                if h1_count == h1_dp and h1_count > 0:
                    if h2_count == h2_dp and h1_count + h2_count > 0:
                        return 'germline'
                    elif h2_count > 0:
                        return 'on_h2'
                    else:
                        return 'germline_h1'
                if h2_count == h2_dp and h2_count > 0:
                    if h1_count > 0:
                        return 'on_h1'
                    else:
                        return 'germline_h2'
          
            if h1_count >= self.read_threshold and h2_count == 0:
                return 'on_h1'
            elif h2_count >= self.read_threshold and h1_count == 0:
                return 'on_h2'
            elif total_passing_alt >= self.read_threshold:
                return 'on_unknown_haplotype'
        
        if total_passing_alt > 0: 
            return 'below_read_threshold'
        elif total_passing_alt == 0 and failed_count > 0:
            return 'on_failed_reads'
        elif sum((h1_count, h2_count, ambiguous_count, failed_count)) == 0:
            return 'absent'
        elif total_fail > 0:
            return 'no_passing_reads'
        elif total_pass + total_fail == 0:
            return 'no_read_data'
        else:
            raise ValueError(f'missing case h1: {h1_count} h2:{h2_count} unphased:{ambiguous_count} failed:{failed_count} passing_dp:{total_pass} failing_dp:{total_fail}')


class ReadFilter:
    """
    Determines whether a sequencing read is of sufficient quality to confidently report allele at variant site.
    Quality filters are predefined for diffrent sequencing platforms.

    Attributes:
        min_basequal: Minimum phred-scaled base quality at variant site.
        min_length: Minimum total length of sequencing read in bp.
        min_gap_id: Minimum gap-adjusted sequence identity (0 - 1.0).
        min_aln_frac: Minimum fraction of read aligned to donor-specific reference (0-1.0). 
        min_dist_to_end: Minimum distance (in bp) between variant site and end of sequencing read. 
    """
    min_basequal: int 
    min_length: int
    min_gap_id: float
    min_aln_frac: float
    min_dist_to_end: int

    def __init__(self, min_basequal, min_length, min_gap_id, min_aln_frac, min_dist_to_end):
        self.min_basequal = min_basequal
        self.min_length = min_length
        self.min_gap_id = min_gap_id
        self.min_aln_frac = min_aln_frac
        self.min_dist_to_end = min_dist_to_end
    

    @staticmethod
    def for_type(data_type: str) -> ReadFilter:
        """
        Define quality thresholds for provided sequencing platform.

        Args:
            data_type: Sequencing platform ('hifi' or 'illumina'). Expect 'ont' in future versions.

        Returns:
            An instance of ReadFilter with platform-specific thresholds.

        Raises:
            KeyError: If unsupported data_type is specified. 
        """
        # For min_gap_id, a 151bp Illumina read with 3 SNVs would have id of 0.98, which could represent a true clustered event. 
        # Larger HiFi reads should not have that density of variation, so the min_gap_id threshold is higher. 
        # For min_aln_frac, large HiFi reads might span SVs, reducing aligned fraction, while 151bp Illumina reads are unlikely to extend that far.
        # For min_dist_to_end, HiFi reads are not expected to have quality drop off at end of read, while Illumina reads are.
        return {
            'hifi': ReadFilter(
                min_basequal=20,
                min_length=10000,
                min_gap_id=0.99,
                min_aln_frac=0.9,
                min_dist_to_end=0
            ),
            'illumina': ReadFilter(
                min_basequal=20,
                min_length=150,
                min_gap_id=0.98,
                min_aln_frac=1,
                min_dist_to_end=10
            ),
        }[data_type]


    def check_read(self, read_dict: Dict[str, str]) -> bool:
        """
        Checks if sequencing read meets minimum quality thresholds. 

        Args:
            read_dict: Dictionary of read quality metrics.

        Returns:
            True is read meets all thresholds, False if it fails any threshold.
        """
        try:
            if int(read_dict['read_length']) < self.min_length:
                return False
            if int(read_dict['base_quality']) < self.min_basequal:
                return False
            if float(read_dict['gap_compressed_identity']) < self.min_gap_id:
                return False
            if float(read_dict['aligned_fraction']) < self.min_aln_frac:
                return False
            if int(read_dict['dist_to_read_end']) < self.min_dist_to_end:
                return False
            return True
        except ValueError:
            return False


class Read:
    read_name: str
    read_data: Dict[str, Dict[str, str]]
    allele: str
    strand: str
    original_haplotype: Optional[str]
    alignment_scores: Dict[str, int]
    uniquely_mapped: str
    pass_filters: Dict[str, bool]
    alternate_alignments: List[str]

    def __init__(self, read_name) -> None:
        self.read_name = read_name
        self.read_data = {}
        self.alignment_scores = {}
        self.pass_filters = {}
        self.alternate_alignments = []
        self.original_haplotype = None


    def add_haplotype(self, read_info: Dict[str, str], read_filter: ReadFilter) -> None:
        haplotype = read_info['chromosome'].split('haplotype')[1][0]
        self.pass_filters[haplotype] = read_filter.check_read(read_info)
        self.read_data[haplotype] = read_info
        

    def format_alternates(self, read: AlignedSegment) -> None:
        alt_alignments = read.get_tag('XA') if read.has_tag('XA') else None
        if not alt_alignments:
            return

        alt_alignments = alt_alignments.split(';')

        for alignment in alt_alignments[:-1]:
            chrom, pos = alignment.split(',')[:2]
            if pos.startswith('-'):
                pos = pos[1:]

            coordinates = (chrom, int(pos))
            if coordinates not in self.alternate_alignments:
                self.alternate_alignments.append(coordinates)


    def annotate_cram(self, read: AlignedSegment) -> None:
        haplotype = read.reference_name.split('haplotype')[1][0]
        self.alignment_scores[haplotype] = read.get_tag('AS')
        self.original_haplotype = str(read.get_tag('oh'))
        self.format_alternates(read)


    def get_allele(self) -> None:
        hap_info = self.read_data.get(self.original_haplotype) or next(iter(self.read_data.values()))
        self.allele = hap_info['base_at_snv']
        self.strand = hap_info['strand']

    def recover_mappings(self, lost_reads: Dict[str, Dict[str, str]], read_filter: ReadFilter) -> None:
        recovered_read = lost_reads.get(self.read_name)

        if recovered_read:
            haplotype = recovered_read['chromosome'].split('haplotype')[1][0]
            if haplotype not in self.read_data:
                self.add_haplotype(recovered_read, read_filter)


    def check_mapping(self):
        # when read fails all quality filters, we don't care how it maps
        if False in self.pass_filters.values():
            self.uniquely_mapped = 'failed'
            return

        aligned_to_original = self.original_haplotype in self.read_data and self.original_haplotype in self.alignment_scores
        homolog_hap = '2' if self.original_haplotype == '1' else '1'
        aligned_to_homolog = homolog_hap in self.read_data and homolog_hap in self.alignment_scores

        # when read aligns to both haplotypes, we can check to see if it's uniquely mapped to one or the other
        if aligned_to_original and aligned_to_homolog:
            h1_score, h2_score = [self.alignment_scores.get(x) for x in ('1', '2')]
            h1_id, h2_id = [round(float(self.read_data.get(x)['gap_compressed_identity']), 3) for x in ('1', '2')]
            same_alignment_score = h1_score == h2_score
            same_gap_compressed_id = h1_id == h2_id

            if h1_score > h2_score:
                self.uniquely_mapped = 'unique_to_h1'
            elif h2_score > h1_score:
                self.uniquely_mapped = 'unique_to_h2'
            else:
                if h1_id > h2_id:
                    self.uniquely_mapped = 'unique_to_h1'
                elif h2_id > h1_id:
                    self.uniquely_mapped = 'unique_to_h2'
                else:
                    self.uniquely_mapped = 'ambiguous'

            return

        # if read only aligns to original haplotype (so no entries in dictionaries for other haplotype)
        if aligned_to_original and not aligned_to_homolog:
            self.uniquely_mapped = f'unique_to_h{self.original_haplotype}'
        # if read only aligns to homologous haplotype (no entries in dictionaries for orig. haplotype)
        elif aligned_to_homolog and not aligned_to_original:
            self.uniquely_mapped = f'unique_to_h{homolog_hap}'
        else:
            self.uniquely_mapped = 'error'


            
class Variant:
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    called_on_both_haps: str
    original_haplotype: int
    homolog_chrom: str = 'NA'
    homolog_pos : str = 'NA'
    homolog_ref_allele: str = 'NA'
    same_orientation: bool
    reads: Dict[str, Read]
    read_summary: Dict[str, Tuple[str, str, str, str, str]]
    collapsed_calls: str
    is_phased: str
    flag: str
    in_tissue: str

    def __init__(self, variant_info: Dict[str, str]) -> None:
        self.variant_id = variant_info['variant_id']
        self.original_haplotype = variant_info['chrom'].split('haplotype')[1][0]
        self.chrom = self.rename_contig(variant_info['chrom'], False)
        self.pos = int(variant_info['pos'])
        self.ref = variant_info['ref']
        self.alt = variant_info['alt']
        self.reads = {}
        self.called_on_both_haps = 'NA'
        self.collapsed_calls = variant_info['collapsed_calls']

        if variant_info['homolog_chrom'] != 'NA':
            self.homolog_chrom = variant_info['homolog_chrom']
            self.homolog_pos = int(variant_info['homolog_pos'])
            self.homolog_ref_allele = variant_info['homolog_ref_allele']


    def rename_contig(self, contig: str, homolog: bool) -> str:
        if homolog:
            haplotype = 2 if self.original_haplotype == 1 else 1
        else:
            haplotype = self.original_haplotype

        if '#' not in contig:
            return f'SMHT004#{haplotype}#{contig}' 
        else:
            return contig


    def check_reverse_complement(self, fasta: FastaFile) -> None:
        ref_allele = fasta.fetch(self.chrom, self.pos-1, self.pos).upper()
        if ref_allele != self.ref and ref_allele == rc(self.ref):
            self.ref = rc(self.ref)
            self.alt = rc(self.alt)


    def add_read(self, read: Dict[str, str], read_filter: ReadFilter):
        read_name = read['read_name']
        if read_name not in self.reads:
            self.reads[read_name] = Read(read_name)

        self.reads[read_name].add_haplotype(read, read_filter)


    def count_reads(self, read_list: List[Tuple[str, str]], allele: str):
        reads = [x for x in read_list if x[0] == allele]
        strand_count = Counter([x[1] for x in reads])

        return str(len(reads)), f"{strand_count['+']};{strand_count['-']}"

        
    def summarize_category(self, read_list: Optional[List[Tuple[str, str]]] = None) -> Tuple[str,str,str,str,str]:
        if not read_list:
            return ('NA', 'NA', 'NA', 'NA', 'NA')
        
        total_reads = str(len(read_list))
        ref_count, ref_strands = self.count_reads(read_list, self.ref)
        alt_count, alt_strands = self.count_reads(read_list, self.alt)

        return (total_reads, ref_count, alt_count, ref_strands, alt_strands)

        
    def summarize_reads(self):
        self.read_summary = {}
        
        alternates = {}
        read_dict = {x: [] for x in ('unique_to_h1', 'unique_to_h2', 'ambiguous', 'failed', 'error')}

        for read_name, read in self.reads.items():
            read_dict[read.uniquely_mapped].append((read.allele, read.strand))

        read_dict['passed'] = read_dict['unique_to_h1'] + read_dict['unique_to_h2'] + read_dict['ambiguous']  
        for category in ('passed', 'unique_to_h1', 'unique_to_h2', 'ambiguous', 'failed'):
             self.read_summary[category] = self.summarize_category(read_dict[category])


    def get_read_count(self, read_class: str, index: int) -> int:
        """
        Args:
            read_class: Any key to read_dict ('passed', 'unique_to_h1', 'unique_to_h2', 'ambiguous', 'failed').
            index: Index of read type in read_dict value (0 = total depth, 1 = ref count, 2 = alt count).

        Returns:
            Number of reads in category, or 0 if category is represented with 'NA'.
        """
        try:
            return int(self.read_summary[read_class][index])
        except ValueError:
            return 0

    
    def confirm_phasing(self):
        unphased_reads = self.get_read_count('ambiguous', 0)
        h1_reads = self.get_read_count('unique_to_h1', 0)
        h2_reads = self.get_read_count('unique_to_h2', 0)

        total_phased_reads = h1_reads + h2_reads
        total_passing_reads = total_phased_reads + unphased_reads
        

        if total_phased_reads == 0 and unphased_reads == 0:
            self.is_phased = 'NA'
        else:
            phased_faction = total_phased_reads / total_passing_reads
            if phased_faction >= 0.9:
                self.is_phased = 'True'
            elif 0 < phased_faction < 0.9:
                self.is_phased = 'partially_phased'
            elif total_phased_reads == 0 and unphased_reads > 0:
                self.is_phased = 'False'
        
    
    def flag(self) -> None:
        if self.alt == self.homolog_ref_allele:
            self.flag = 'germline_heterozygote'
            return

        if self.is_phased == 'False':
            self.flag = 'unknown_haplotype'
            return

        h1_alt_count = self.get_read_count('unique_to_h1', 2)
        h2_alt_count = self.get_read_count('unique_to_h2', 2)
        if self.is_phased == 'True' and h1_alt_count > 0 and h2_alt_count > 0:
            self.flag = 'recurrent_error'
            return 

        passing_read_count = self.get_read_count('passed', 0)
        if passing_read_count == 0:
            self.flag = 'no_passing_reads'
            return 

        self.flag = 'PASS'


    def print_output(self, outfile: TextIO) -> None:
        output_data = [self.chrom, str(self.pos), self.variant_id, self.ref, self.alt,
            self.homolog_chrom, str(self.homolog_pos), self.homolog_ref_allele,
            self.is_phased, self.flag, self.in_tissue]

        for read_class in ('passed', 'unique_to_h1', 'unique_to_h2', 'ambiguous', 'failed'):
            output_data += self.read_summary[read_class]

        print('\t'.join(output_data), file = outfile)

          
class VariantIndex:
    variants: Dict[str, Variant]
    lost_reads: Dict[str, Dict[str, str]]

    def __init__(self, variants: Dict[str, Variant]) -> None:
        self.variants = variants

    @staticmethod
    def import_file(variant_file: str) -> VariantIndex:
        print('importing variants...')
        variants = {}
        coordinate_lookup = defaultdict(list)

        with open(variant_file, 'r') as infile:
            header = infile.readline().rstrip().split('\t')

            for line in infile:
                variant = line.rstrip().split('\t')
                variant_info = {k: v for k, v in zip(header, variant)}

                if '#mt#' not in variant_info['chrom']:
                    var_object = Variant(variant_info) 
                    variants[variant_info['variant_id']] = var_object
                    coordinate_lookup[(var_object.chrom, str(var_object.pos))].append(var_object.variant_id)
                    if var_object.homolog_chrom != 'NA':
                        coordinate_lookup[(var_object.homolog_chrom, str(var_object.homolog_pos))].append(var_object.variant_id)

        return VariantIndex(variants), coordinate_lookup


    def import_read_data(self, read_file: str, coordinate_lookup:  Dict[Tuple[str, str], str], read_filter: ReadFilter) -> None:
        print('importing read data...')
        self.lost_reads = {}

        with gzip.open(read_file, 'rt') as infile:
            header = infile.readline().rstrip().split('\t')
            header_idx = {k:v for v,k in enumerate(header)}

            for line in infile:
                read = line.rstrip().split('\t')
                try:
                    variant_ids = coordinate_lookup[(read[header_idx['chromosome']], read[header_idx['snv_position']])]
                except KeyError:
                    self.lost_reads[read[header_idx['read_name']]] = dict(zip(header, read))
                    continue 
                
                for variant_id in variant_ids:
                    if variant_id not in self.variants:
                        continue
                    self.variants[variant_id].add_read(dict(zip(header, read)), read_filter)


    def query_cram(self, variant: Variant, cram: AlignmentFile, chromosome: str, position: int) -> None:
        for read in cram.fetch(chromosome, position - 1, position):
            read_name = read.query_name

            if read_name not in variant.reads:
                continue

            variant.reads[read_name].annotate_cram(read)


    def import_cram_data(self, cram_file: str, read_filter: ReadFilter) -> None:
        cram = AlignmentFile(cram_file, "rc")
        read_dict = {}

        print("importing cram data...")
        for variant_id, variant in tqdm(self.variants.items()):
            self.query_cram(variant, cram, variant.chrom, variant.pos)
            if variant.homolog_chrom != 'NA' and isinstance(variant.homolog_pos, int) :
                self.query_cram(variant, cram, variant.homolog_chrom, variant.homolog_pos)
            for read_name, read in variant.reads.items():
                read.get_allele()
                read.recover_mappings(self.lost_reads, read_filter)
                read.check_mapping()
                read_dict[read_name] = read.uniquely_mapped

        return read_dict

    
    def process_read_data(self, presence_filter: PresenceFilter, output_file: str) -> None:
        output_header = ('chrom', 'pos', 'id', 'ref', 'alt', 
            'homolog_chrom', 'homolog_pos', 'homolog_ref_allele',
            'phased', 'haplotype_flag', 'in_tissue',
            'passed_dp', 'passed_ref', 'passed_alt', 'passed_ref_strands', 'passed_alt_strands',
            'h1_dp', 'h1_ref', 'h1_alt', 'h1_ref_strands', 'h1_alt_strands',
            'h2_dp', 'h2_ref', 'h2_alt', 'h2_ref_strands', 'h2_alt_strands',
            'ambiguous_dp', 'ambiguous_ref', 'ambiguous_alt', 'ambiguous_ref_strands', 'ambiguous_alt_strands',
            'failed_dp', 'failed_ref', 'failed_alt', 'failed_ref_strands', 'failed_alt_strands'
        )

        with open(output_file, 'w') as outfile:
            print('processing read data...')
            print('\t'.join(output_header), file = outfile)
            for variant_id, variant in tqdm(self.variants.items()):
                variant.summarize_reads()
                variant.confirm_phasing()
                variant.flag()
                variant.in_tissue = presence_filter.determine_presence(variant)
                variant.print_output(outfile)
            

def rc(allele):
    if allele == 'A':
        return 'T'
    elif allele == 'C':
        return 'G'
    elif allele == 'G':
        return 'C'
    elif allele == 'T':
        return 'A'
    else:
        return None


if __name__ == "__main__":
    # to run from the command line, use commented out lines
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=True, help="Candidate mutation file")
    ap.add_argument("-r", "--read_data", required=True, help="File of read data with quality metrics")
    ap.add_argument("-p", "--platform", required=True, help="Sequencing platform")
    ap.add_argument("-c", "--cram", required=True, help="cram file")
    ap.add_argument("-f", "--fasta", required=True, help="fasta file")
    ap.add_argument("-o", "--output", required=True, help="Output file")
    
    args = ap.parse_args()

    read_filter = ReadFilter.for_type(args.platform)
    presence_filter = PresenceFilter.for_type(args.platform)


    variants, coordinate_lookup = VariantIndex.import_file(args.variants)

    variants.import_read_data(args.read_data, coordinate_lookup, read_filter)
    read_dict = variants.import_cram_data(args.cram, read_filter)

    variants.process_read_data(presence_filter, args.output)
    print('done!')
