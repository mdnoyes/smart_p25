from __future__ import annotations
from typing import Dict, List, Tuple, Optional, TextIO
from scipy.stats import binomtest, fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from collections import Counter
import argparse
from dataclasses import dataclass

@dataclass
class ReadData:
    haplotype: str
    dp: int
    ref: int
    alt: int
    ref_strands: str
    alt_strands: str

    @staticmethod
    def from_read_counts(variant_info: Dict[str, str], haplotype: str) -> ReadData:
        dp, ref, alt = [
            ReadData.to_int(variant_info[f'{haplotype}_{x}']) 
            for x in ('dp', 'ref', 'alt')
        ]
        ref_strands, alt_strands = [
            variant_info[f'{haplotype}_{x}'] 
            for x in ('ref_strands', 'alt_strands')
        ]
        return ReadData(haplotype, dp, ref, alt, ref_strands, alt_strands)
    

    @staticmethod
    def to_int(val):
        return int(val) if val != 'NA' else 0


class Tissue:
    name: str
    platforms: List[str]
    read_data: Dict[str, ReadData]
    in_tissue: str
    phased_by: str
    haplotype: str
    concordant_vaf: str
    is_recurrent_error: str

    def __init__(self, tissue_name: str, platforms: List[str], variant_info: Dict[str, str]) -> None:
        self.name = tissue_name
        self.platforms = platforms
        self.in_tissue = variant_info['in_tissue']
        self.phased_by = variant_info['phased_by']
        self.haplotype = variant_info['haplotype']
        self.concordant_vaf = variant_info['concordant_vaf']
        self.is_recurrent_error = variant_info['is_recurrent_error'] == 'True'
        self.read_data = {}

        if 'hifi' in platforms:
            for haplotype in ('hifi_passed', 'hifi_h1', 'hifi_h2', 'hifi_ambiguous'):
                self.read_data[haplotype] = ReadData.from_read_counts(variant_info, haplotype)
        if 'illumina' in platforms:
            for haplotype in ('illumina_passed', 'illumina_h1', 'illumina_h2', 'illumina_ambiguous'):
                self.read_data[haplotype] = ReadData.from_read_counts(variant_info, haplotype)


    def format_read_data(self, platform: str):
        try:
            h1_counts = self.read_data[f'{platform}_h1']
            h2_counts = self.read_data[f'{platform}_h2']
            ambiguous_counts = self.read_data[f'{platform}_ambiguous']
        except KeyError:
            return 'NA'

        output = []

        for hap_counts in (h1_counts, h2_counts, ambiguous_counts):
            read_counts = ';'.join([str(x) for x in (hap_counts.dp, hap_counts.ref, hap_counts.alt)])
            hap_name = hap_counts.haplotype.split('_')[1]
            output.append(f'{hap_name}:{read_counts}')

        return '|'.join(output)


class Variant:
    variant_id: str
    variant_info: Dict[str, str]
    tissues: Dict[str, Tissue]
    is_germline: bool
    is_recurrent_error: bool
    haplotype: str
    unadjusted_binom_pvals: Dict[str, float]
    phased_by: str
    discordant_tissues: List[str]
    tissues_with_variant: List[str]
    low_confidence_tissues: List[str]
    pass_hifi_binom_test: str
    pass_ilm_binom_test: str

    def __init__(self, variant_info: Dict[str, str]) -> None:
        self.variant_id = variant_info['variant_id']
        self.variant_info = variant_info
        self.tissues = {}
        self.is_germline = False
        self.unadjusted_pvals = {}
        self.adjusted_pvals = {}
        self.haplotype = 'NA'
        self.discordant_tissues = []
        self.pass_hifi_binom_test = 'NA'
        self.pass_ilm_binom_test = 'NA'


    def add_tissue(self, tissue_name: str, tissue_platforms: List[str], tissue_data: Dict[str, str]):
        if tissue_data['is_germline'] == 'True':
            self.is_germline = True

        self.tissues[tissue_name] = Tissue(tissue_name, tissue_platforms, tissue_data)

    
    def in_haplotype(self, haplotype_counts: Dict[str, int], haplotype: str, threshold: int) -> bool:
        homolog = 'h2' if haplotype == 'h1' else 'h1'

        hap_count = haplotype_counts.get(haplotype, 0)
        homolog_count = haplotype_counts.get(homolog, 0)

        if hap_count >= 3 and homolog_count <= threshold:
            return True
        
        return False
       

    def examine_mismatched_haplotypes(self, haplotype_counts: Dict[str, int]) -> None:
        if self.in_haplotype(haplotype_counts, 'h1', 1):
            assumed_hap = 'h1'
        elif self.in_haplotype(haplotype_counts, 'h2', 1):
            assumed_hap = 'h2'
        else:
            self.is_recurrent_error = True
            self.haplotype = 'NA'
            return

        discordant_tissue = next(iter([
            tissue.name 
            for tissue in self.tissues.values() 
            if tissue.haplotype not in (assumed_hap, 'NA')
        ]), None)

        if discordant_tissue:
            self.discordant_tissues.append(discordant_tissue)

        self.haplotype = assumed_hap
        self.is_recurrent_error = False
        

    def examine_recurrent_error(self, haplotype_counts: Dict[str, int]) -> None:
        discordant_tissue = next(iter([
            tissue.name 
            for tissue in self.tissues.values() 
            if tissue.is_recurrent_error
        ]), None)
            
        if self.in_haplotype(haplotype_counts, 'h1', 0):
            self.is_recurrent_error = False
            self.haplotype = 'h1'
            if discordant_tissue:
                self.discordant_tissues.append(discordant_tissue)
        elif self.in_haplotype(haplotype_counts, 'h2', 0):
            self.is_recurrent_error = False
            self.haplotype = 'h2'
            if discordant_tissue:
                self.discordant_tissues.append(discordant_tissue)
        else:
            self.is_recurrent_error = True
            self.haplotype = 'NA'


    def compare_haplotypes(self):
        # possible tissue.phased_by values: 'hifi', 'illumina', 'hifi_illumina', 'NA'
        phased_in_tissue = [
            tissue.phased_by for tissue in self.tissues.values() 
            if tissue.phased_by != 'NA'
        ]

        phasing_platforms = sorted(list(set([
            platform for platforms in phased_in_tissue for platform in platforms.split('_')
        ])))

        self.phased_by = '_'.join(phasing_platforms) if phasing_platforms else 'NA'


        # possible tissue.is_recurrent_error values: True, False
        recurrent_in_tissue = [
            tissue.is_recurrent_error for tissue in self.tissues.values() 
            if tissue.phased_by != 'NA'
        ]

        # possible tissue.haplotype values: 'h1', 'h2', 'NA'
        haplotype_counts = Counter([
            tissue.haplotype for tissue in self.tissues.values() 
            if tissue.haplotype != 'NA'
        ])

        total_phased_tissues = sum(haplotype_counts.values())
        recurrent_error_count = sum(recurrent_in_tissue)

        # not phased in any tissue, so there's nothing we can do
        if recurrent_error_count == 0 and total_phased_tissues == 0:
            self.is_recurrent_error = False
            self.haplotype = 'NA'
            return

        # if it's a recurrent error in only one tissue, we might be able to rescue it
        if recurrent_error_count == 1:
            self.examine_recurrent_error(haplotype_counts)
            return
        # if it's recurrent on more than one tissue, it's an error
        elif recurrent_error_count > 1:
            self.is_recurrent_error = True
            self.haplotype = 'NA'
            return

        # if phased to both haplotypes across tissues, we might be able to rescue it
        if 'h1' in haplotype_counts and 'h2' in haplotype_counts:
            self.examine_mismatched_haplotypes(haplotype_counts)
            return 

        # otherwise, it should be a clean assignment to h1 or h2
        if 'h1' in haplotype_counts and 'h2' not in haplotype_counts:
            self.is_recurrent_error = False
            self.haplotype = 'h1'
        elif 'h2' in haplotype_counts and 'h1' not in haplotype_counts:
            self.is_recurrent_error = False
            self.haplotype = 'h2'
        else:
            raise ValueError('you have not handled this phasing case')


    def get_read_counts(self, platform: str, haplotype: str, allele: str, in_tissue: bool = False) -> int:
            """
            Get the total number of passing reads that support an allele across all tissues for a given platform.

            Args:
                platform: Sequencing platform, either 'hifi' or 'illumina'.
                haplotype: Haplotype, 'passed', 'h1', or 'h2'.
                allele: Allele of interest, one of 'ref', 'alt', or 'dp' for all passing reads.
            """
            read_type = f'{platform}_{haplotype}'

            if in_tissue:
                tissue_dict = {k:v for k,v in  self.tissues.items() if k in self.tissues_with_variant}
            else:
                tissue_dict = self.tissues

            tissue_reads = [tissue.read_data.get(read_type) for tissue in tissue_dict.values() if tissue.name not in self.discordant_tissues]

            return sum(getattr(tissue, allele) for tissue in tissue_reads if tissue is not None)


    def binomial_test(self, allele_count: int, total_count: int, expected_VAF: float, alternative: str) -> Optional[float]:
        if allele_count >= 0 and total_count > 0: 
            return binomtest(allele_count, total_count, expected_VAF, alternative=alternative).pvalue
        return None


    def run_binomial_tests(self, platform: str):
        # do not run on recurrent errors, this should help with multiple testing 
        if self.is_recurrent_error:
            return

        total_h1, total_h2 = [self.get_read_counts(platform, x, 'dp') for x in ('h1', 'h2')]

        if total_h1 > 0 and total_h2 > 0 and self.haplotype in ('h1', 'h2') and platform == 'hifi':
            is_phased = True
            mutant_hap = self.haplotype
        else:
            is_phased = False
            mutant_hap = 'passed'
                
        ref, alt, total = [self.get_read_counts(platform, mutant_hap, x) for x in ('ref', 'alt', 'dp')]
        if total == 0:
            return

        if is_phased and not self.is_germline and platform == 'hifi':
            # looking only at reads from mutant haplotype
            # germline event (expected alt VAF = 1 - error rate, expected ref VAF = error rate) 
            # somatic event (expected alt VAF < 1, expected ref VAF < 1)
            mutant_alt_pval = self.binomial_test(alt, total, 0.9, 'less')
            if mutant_alt_pval  is not None:
                self.unadjusted_pvals[f'{platform}_mutant_alt_pval'] = mutant_alt_pval
        else:
            # looking at reads from both haplotypes
            # germline event (expected alt VAF = 0.5, expected ref VAF = 0.5) 
            # somatic event (expected alt VAF < 1, expected ref VAF > 1)
            pooled_alt_pval = self.binomial_test(alt, total, 0.5, 'two-sided')
            if pooled_alt_pval is not None:
                self.unadjusted_pvals[f'{platform}_pooled_alt_pval'] = pooled_alt_pval
                

    def interpret_binomial_test(self, platform: str):
        pass_test = None

        # if we have phased reads
        if f'{platform}_mutant_alt_pval' in self.adjusted_pvals:
            pval = self.adjusted_pvals[f'{platform}_mutant_alt_pval']
            pass_test = str(pval < 0.01)

        # if we have unphased reads
        if f'{platform}_pooled_alt_pval' in self.adjusted_pvals:
            pval = self.adjusted_pvals[f'{platform}_pooled_alt_pval']
            
            if self.is_germline:      
                # heterozygous gentoype, so we have to be more creative here - we're trying to figure out mutation direction          
                if pval < 0.01:
                    pass_test = 'True'

                    h = self.haplotype if self.haplotype != 'NA' else 'passed'

                    alt, dp = [self.get_read_counts(platform, h, x, True) for x in ('alt', 'dp')]
                    alt_vaf = alt / dp if dp > 0 else 1 
                else:
                    # vaf is not significantly different from 0.5, looks germline
                    pass_test = 'False' 
            else:
                # this is just the normal case standard binomial test 
                # homozygous ref genotype, unphased reads, is VAF significantly different from 0.5?
                pass_test = str(pval < 0.01)

        if pass_test is None:
            return

        if platform == 'hifi':
            self.pass_hifi_binom_test = pass_test
        else:
            self.pass_ilm_binom_test = pass_test


    def get_strand_counts(self, platform: str) -> int:
        tissue_reads = [tissue.read_data.get(f'{platform}_passed') for tissue in self.tissues.values() if tissue.name in self.tissues_with_variant]

        ref_strands = [tissue.ref_strands.split(';') for tissue in tissue_reads if tissue is not None]
        ref_strands = [x for x in ref_strands if x != ['NA']]
        alt_strands = [tissue.alt_strands.split(';') for tissue in tissue_reads if tissue is not None]
        alt_strands = [x for x in alt_strands if x != ['NA']]

        ref_plus, ref_minus = [
            sum(int(x) for x in col if x != 'NA') 
            for col in zip(*ref_strands)
        ] if ref_strands else (0, 0)

        alt_plus, alt_minus = [
            sum(int(x) for x in col if x != 'NA') 
            for col in zip(*alt_strands)
        ] if alt_strands else (0, 0)

        contingency_table = [
            [ref_plus, ref_minus],
            [alt_plus, alt_minus]
        ]

        return contingency_table


    def run_strand_bias_test(self, platform: str):
        # do not run on recurrent errors
        if self.is_recurrent_error:
            return

        contingency_table = self.get_strand_counts(platform)

        # confirm that the alternate allele is observed in passing reads
        if sum(contingency_table[1]) == 0:
            return

        odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')

        self.unadjusted_pvals[f'{platform}_strand_bias'] = p_value


    def interpret_strand_bias_test(self):
        illumina_pval = self.adjusted_pvals.get('illumina_strand_bias')
        self.strand_bias = []

        if illumina_pval is None:
            self.strand_bias.append('NA')
            return

        if illumina_pval is not None and illumina_pval < 0.01:
            self.strand_bias.append('FAIL')
            
        if not self.strand_bias:
            self.strand_bias.append('PASS')

    
    def get_confident_tissues(self, confidence_level: str) -> List[str]:
        tissue_list = [
            tissue.name for tissue in self.tissues.values() 
            if tissue.name not in self.discordant_tissues and
            tissue.in_tissue == confidence_level
        ]

        return tissue_list


    def single_platform_support(self, tissue: str, platform: str):
        threshold = 3 if platform == 'hifi' else 5
        read_info = self.tissues[tissue].read_data.get(f'{platform}_passed')
    
        if read_info is None:
            return False
        
        return read_info.alt >= threshold


    def define_tissues(self):
        self.tissues_with_variant = []
        self.low_confidence_tissues = []

        # tissue.in_tissue possible values: 
        # 'high_confidence', 'medium_confidence', 'low_confidence', 'absent'
        
        confidence_counts = Counter([
            tissue.in_tissue for tissue in self.tissues.values() 
            if tissue.name not in self.discordant_tissues
        ])

        self.low_confidence_tissues += self.get_confident_tissues('low_confidence')
        
        # if it is only low confidence or absent, then we will consider it absent
        if 'high_confidence' not in confidence_counts and 'medium_confidence' not in confidence_counts:
            return

        # if there are no high confidence tissues (requires HiFi and Illumina support)
        if confidence_counts['high_confidence'] == 0:
            # if there are are at least two medium confidence tissues
            if confidence_counts['medium_confidence'] >= 2:
                # if they have sufficient read depth in the supporting platforms
                self.tissues_with_variant += self.get_confident_tissues('medium_confidence')
                return
            # if there is only one medium confidence tissues, throw them in with the low confidence
            # this is mostly just for bookkeeping, sites with no tissues_with_variant will be excluded later
            else:
                single_tissue = self.get_confident_tissues('medium_confidence')[0]
                # if there's only 1 medium confidence tissue, but it has strong hifi support, let's rescue it
                if self.single_platform_support(single_tissue, 'hifi') or self.single_platform_support(single_tissue, 'illumina'):
                    self.tissues_with_variant.append(single_tissue)
                else:
                    self.low_confidence_tissues.append(single_tissue)
                return

        # if we see at least one high confidence tissue, all high and medium confidence tissues should be retained
        elif confidence_counts['high_confidence'] > 0:
            self.tissues_with_variant += self.get_confident_tissues('high_confidence')
            for tissue in self.get_confident_tissues('medium_confidence'):
                if self.single_platform_support(tissue, 'hifi') or self.single_platform_support(tissue, 'illumina'):
                    self.tissues_with_variant.append(tissue)
            return

        raise ValueError('there is a confidence combination you did not account for')


    def check_platform_support(self, platform):
        supporting_reads = 0

        if self.haplotype not in ('h1', 'h2'):
            hap = 'passed'
        else:
            hap = self.haplotype

        for tissue in self.tissues_with_variant:
            if f'{platform}_passed' not in self.tissues[tissue].read_data:
                continue

            if hap != 'passed':
                alt_count = self.tissues[tissue].read_data[f'{platform}_{hap}'].alt
                if alt_count != 0:
                    supporting_reads += alt_count
                    continue
            
            supporting_reads += self.tissues[tissue].read_data[f'{platform}_ambiguous'].alt

        return str(supporting_reads > 0)
        

    def print(self, output_header: List[str], outfile: TextIO):
        output = [
            self.variant_info[x] 
            for x in output_header
            if x in self.variant_info
        ]

        output += [str(self.is_germline), self.phased_by, self.haplotype, 
             str(self.pass_hifi_binom_test), str(self.pass_ilm_binom_test), ';'.join(self.strand_bias)]

        output += [self.check_platform_support('hifi'), self.check_platform_support('illumina')]
        
        for tissue_list in (self.tissues_with_variant, self.low_confidence_tissues, self.discordant_tissues):
            if tissue_list:
                tissue_list.sort()
                output.append(';'.join(tissue_list))
            else:
                output.append('NA')

        for tissue, platform in [x.split('_') for x in output_header[25:]]:
            if tissue in self.tissues:
                output.append(self.tissues[tissue].format_read_data(platform))
            else:
                output.append('NA')

        print('\t'.join(output), file = outfile)

        return


class VariantIndex:
    variants: Dict[str, Variant]

    def __init__(self, variant_file: str) -> None:
        print('importing variants...')
        self.variants = {}

        with open(variant_file, 'r') as infile:
            header = infile.readline().rstrip().split('\t')

            for line in infile:
                variant = line.rstrip().split('\t')
                variant_info = {k:v for k,v in zip(header, variant)}

                self.variants[variant_info['variant_id']] = Variant(variant_info) 


    def import_tissue(self, tissue_name: str, tissue_file: str):
        with open(tissue_file, 'r') as infile:
            header = infile.readline().rstrip().split('\t')
            platforms = [p for p in ('hifi', 'illumina') if f'{p}_in_tissue' in header]

            for line in infile:
                variant = line.rstrip().split('\t')
                variant_info = {k: v for k, v in zip(header, variant)}

                try:
                    self.variants[variant_info['variant_id']].add_tissue(tissue_name, platforms, variant_info)
                except KeyError:
                    continue


    def multiple_test_correction(self, test_type: str):
            variant_ids = list([k for k, v in self.variants.items() if test_type in v.unadjusted_pvals.keys()])
            unadjusted_pvalues = [self.variants[key].unadjusted_pvals[test_type] for key in variant_ids]
            
            adjusted_pvalues = multipletests(unadjusted_pvalues, method='fdr_bh')[1]

            for variant_id, pval in zip(variant_ids, adjusted_pvalues):
                self.variants[variant_id].adjusted_pvals[test_type] = pval


    def combine_tissue_data(self):
        for variant_id, variant in self.variants.items():
            variant.compare_haplotypes()
            variant.define_tissues()
            variant.run_binomial_tests('hifi')
            variant.run_binomial_tests('illumina')

        self.multiple_test_correction('hifi_mutant_alt_pval')
        self.multiple_test_correction('hifi_pooled_alt_pval')
        self.multiple_test_correction('illumina_pooled_alt_pval')

        for variant_id, variant in self.variants.items():
            variant.interpret_binomial_test('hifi')
            variant.interpret_binomial_test('illumina')
            variant.run_strand_bias_test('illumina')

        self.multiple_test_correction('illumina_strand_bias')

        for variant_id, variant in self.variants.items():
            variant.interpret_strand_bias_test()

    def print_variants(self, output_header: List[str], output_file: str):
        with open(output_file, 'w') as outfile:
            print('\t'.join(output_header), file = outfile)
            for variant_id, variant in self.variants.items():
                if not variant.is_recurrent_error and len(variant.tissues_with_variant) > 0:
                # if len(variant.tissues_with_variant) > 0:
                    variant.print(output_header, outfile)



if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=False, help="variant file")
    ap.add_argument("-t", "--tissue", action='append', nargs=2, required=True, metavar=('NAME', 'FILE'), help="Tissue name and read data file")
    ap.add_argument("-o", "--output", required=True, help="Output")

    args = ap.parse_args()  

    variants = VariantIndex(args.variants)

    output_header = [
        'chrom', 'pos', 'variant_id', 'ref', 'alt', 
        'homolog_chrom', 'homolog_pos', 'homolog_ref_allele', 
        'collapsed_calls', 'flipped_alleles', 'unique', 'satellite', 'segdup', 'repeat',
        'is_germline', 'phased_by', 'haplotype', 
        'hifi_binom_test', 'illumina_binom_test', 'strand_bias',
        'supported_in_hifi', 'supported_in_illumina',
        'tissues_with_variant', 'low_confidence_tissues', 'discordant_tissues'
    ]


    for tissue_name, tissue_file in args.tissue:
        print(f'importing {tissue_name}...')
        variants.import_tissue(tissue_name, tissue_file)
        output_header += [f'{tissue_name}_hifi', f'{tissue_name}_illumina']

    variants.combine_tissue_data()

    print('printing variants...')
    variants.print_variants(output_header, args.output)

    print('done')

