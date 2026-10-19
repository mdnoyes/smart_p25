from __future__ import annotations
from typing import Dict, List, Tuple, Optional
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import argparse

class Read_Data:
    haplotype: str
    dp: int
    ref: int
    alt: int
    ref_strands: str
    alt_strands: str

    def __init__(self, variant_info: Dict[str, str], haplotype: str):
        self.haplotype = haplotype
        self.dp, self.ref, self.alt = [
            Read_Data.to_int(variant_info[f'{haplotype}_{x}']) 
            for x in ('dp', 'ref', 'alt')
        ]
        self.ref_strands, self.alt_strands = [
            variant_info[f'{haplotype}_{x}'] 
            for x in ('ref_strands', 'alt_strands')
        ]


    @staticmethod
    def to_int(val):
        return int(val) if val != 'NA' else 0


class Platform_Data:
    phased: str
    haplotype_flag: str
    in_tissue: str
    read_data: Dict[str, Read_Data]

    def __init__(self, platform_info: Dict[str, str]) -> None:
        self.phased, self.haplotype_flag, self.in_tissue = [
            platform_info[x] 
            for x in ('phased', 'haplotype_flag', 'in_tissue')
        ]

        self.read_data = {}
        for haplotype in ('passed', 'h1', 'h2', 'ambiguous', 'failed'):
            self.read_data[haplotype] = Read_Data(platform_info, haplotype)


    def summarize_reads(self):
        output = []
        for haplotype in ('passed', 'h1', 'h2', 'ambiguous'):
            output += [str(x) for x in 
                (self.read_data[haplotype].dp, 
                self.read_data[haplotype].ref, 
                self.read_data[haplotype].alt, 
                self.read_data[haplotype].ref_strands, 
                self.read_data[haplotype].alt_strands)]

        return output


class Variant:
    variant_id: str
    variant_info: List[str]
    platforms: Dict[str, Platform_Data]
    is_germline: bool
    is_recurrent_error: bool
    phased_by: str
    haplotype: str
    in_tissue: str
    vaf_pvalue: Optional[str]
    concordant_vaf: str

    def __init__(self, variant_info: Dict[str, str]) -> None:
        self.variant_id = variant_info['variant_id']
        self.variant_info = variant_info
        self.platforms = {}
        self.is_germline = False
        self.is_recurrent_error = False
        self.phased_by = 'NA'
        self.haplotype = 'NA'
        self.vaf_pvalue = None
        self.concordant_vaf = 'NA'
        self.in_tissue = 'absent'


    def add_platform(self, platform: str, platform_info: Dict[str, str]):
        self.platforms[platform] = Platform_Data(platform_info)


    def process_single_platform(self):
        platform = next(iter(self.platforms))

        if self.platforms[platform].haplotype_flag == 'germline_heterozygote':
            self.is_germline = True
            self.process_germline_variant()
        elif self.platforms[platform].haplotype_flag == 'recurrent_error':
            self.in_tissue = 'low_confidence'
            self.is_recurrent_error = True
            self.phased_by = platform

            return

        if self.platforms[platform].phased in ('True', 'partially_phased'):
            self.phased_by = platform

        if self.platforms[platform].in_tissue == 'on_h1':
            self.haplotype = 'h1'
            self.in_tissue = 'medium_confidence'
        elif self.platforms[platform].in_tissue == 'on_h2':
            self.haplotype = 'h2'
            self.in_tissue = 'medium_confidence'
        elif self.platforms[platform].in_tissue == 'on_unknown_haplotype':
            self.in_tissue = 'medium_confidence'
        elif self.platforms[platform].in_tissue == 'below_read_threshold':
            self.in_tissue = 'low_confidence'
        elif self.platforms[platform].in_tissue == 'on_failed_reads':
            self.in_tissue = 'low_confidence'
        elif self.platforms[platform].in_tissue in ('germline', 'h1_germline', 'h2_germline'):
            self.in_tissue = 'absent'
        else:
            self.in_tissue = 'absent'
        

    def resolve_phase(self, hifi_haplotype: str, illumina_haplotype: str) -> None:
        unphased_options = ('absent', 'on_failed_reads', 'on_unknown_haplotype', 'below_read_threshold', 'germline')
        phased_options = ('on_h1', 'on_h2', 'germline_h1', 'germline_h2')

        # if it is a recurrent error, it is useful to know which tissue(s) told us that
        if self.is_recurrent_error:
            hifi_active = hifi_haplotype in ('on_h1', 'on_h2', 'on_unknown_haplotype')
            ilm_active = illumina_haplotype in ('on_h1', 'on_h2', 'on_unknown_haplotype')

            if hifi_active and ilm_active:
                self.phased_by = 'hifi_illumina'
            elif hifi_active:
                self.phased_by = 'hifi'
            elif ilm_active:
                self.phased_by = 'illumina'
            else:
                self.phased_by = 'NA'
            return

        # now, our only options are either 2x 'on_h1' or 2x 'on_h2' 
        # or 'on_h*' with one of: 'absent', 'on_failed_reads', 'on_unknown_haplotype'
        if 'on_h1' in (hifi_haplotype, illumina_haplotype):
            self.haplotype = 'h1'
        elif 'on_h2' in (hifi_haplotype, illumina_haplotype):
            self.haplotype = 'h2'
        # in case there's a mistake and it hasn't been assigned to a haplotype, 
        # quit before assigning phasing platform
        else:
            return

        if hifi_haplotype in phased_options and illumina_haplotype in phased_options:
            self.phased_by = 'hifi_illumina' 
        elif illumina_haplotype in unphased_options and hifi_haplotype in phased_options:
            self.phased_by = 'hifi'
        elif hifi_haplotype in unphased_options and illumina_haplotype in phased_options:
            self.phased_by = 'illumina'

    
    def process_germline_variant(self) -> None:
        haplotype = 'h1' if 'haplotype1' in self.variant_info['chrom'] else 'h2'
        homolog_haplotype = 'h2' if 'haplotype1' in self.variant_info['chrom'] else 'h1'

        for platform in self.platforms:
            in_tissue = self.platforms[platform].in_tissue
            if in_tissue == f'on_{homolog_haplotype}':
                self.platforms[platform].in_tissue = 'absent'
            elif 'germline' in in_tissue:
                if self.platforms[platform].read_data['ambiguous'].alt > 0:
                    self.platforms[platform].in_tissue = 'on_unknown_haplotype'
                else:
                    self.platforms[platform].in_tissue = 'germline'
            elif in_tissue == 'on_unknown_haplotype':
                self.platforms[platform].in_tissue = f'on_{haplotype}'


    def compare_phase(self) -> None:
        # possible haplotype_flag values:
        # germline_heterozygote, no_passing_reads, PASS, recurrent_error, unknown_haplotype
        hifi_flag, ilm_flag = [self.platforms[x].haplotype_flag for x in ('hifi', 'illumina')]

        if 'germline_heterozygote' in (hifi_flag, ilm_flag):
            self.is_germline = True
            self.process_germline_variant()
    
        # possible phased values:
        # False, NA, partially_phased, True
        hifi_is_phased, ilm_is_phased = [self.platforms[x].phased for x in ('hifi', 'illumina')]

        # possible in_tissue values:
        # absent, on_failed_reads, on_h1, on_h2, on_unknown_haplotype, below_read_threshold
        hifi_in_tissue, ilm_in_tissue = [self.platforms[x].in_tissue for x in ('hifi', 'illumina')]

        if 'recurrent_error' in (hifi_flag, ilm_flag) or {'on_h1', 'on_h2'} <= {hifi_in_tissue, ilm_in_tissue}:
            self.is_recurrent_error = True

        if any(p in ('True', 'partially_phased') for p in (hifi_is_phased, ilm_is_phased)):
            self.resolve_phase(hifi_in_tissue, ilm_in_tissue)


    def establish_presence(self) -> None:
        # in the case where it's a recurrent error on both tissues, it could still end up looking high confidence
        # to get ahead of that, let's call all recurrent errors low confidence and move on
        if self.is_recurrent_error:
            self.in_tissue = 'low_confidence'
            return

        # possible status values:
        # absent, on_failed_reads, on_h1, on_h2, on_unknown_haplotype, below_read_threshold
        hifi_status, ilm_status = [self.platforms[x].in_tissue for x in ('hifi', 'illumina')]

        confirmed_in_tissue = ('on_h1', 'on_h2', 'on_unknown_haplotype')
        hifi_confirmed = hifi_status in confirmed_in_tissue
        ilm_confirmed = ilm_status in confirmed_in_tissue

        if hifi_confirmed:
            if ilm_confirmed or ilm_status == 'below_read_threshold':
                self.in_tissue = 'high_confidence'
            else:
                self.in_tissue = 'medium_confidence'
            return
        
        if ilm_confirmed:
            if hifi_status == 'below_read_threshold':
                self.in_tissue = 'high_confidence'
            else:
                self.in_tissue = 'medium_confidence'
            return

        if ilm_status == 'below_read_threshold' or hifi_status == 'below_read_threshold':
            self.in_tissue = 'low_confidence'
            return

        self.in_tissue = 'absent'


    def compare_vaf(self):
        hifi_alt_count, ilm_alt_count = [self.platforms[x].read_data['passed'].alt for x in ('hifi', 'illumina')]
        hifi_dp, ilm_dp = [self.platforms[x].read_data['passed'].dp for x in ('hifi', 'illumina')]

        if hifi_dp == 0 or ilm_dp == 0:
            return

        table = [
            [hifi_alt_count, hifi_dp - hifi_alt_count],
            [ilm_alt_count, ilm_dp - ilm_alt_count]
        ]

        self.vaf_pvalue = fisher_exact(table, alternative='two-sided')[1]


    def combine_platforms(self):
        if len(self.platforms) == 0:
            pass
        elif len(self.platforms) == 1:
            self.process_single_platform()
        else:
            self.compare_phase()
            self.establish_presence()
            self.compare_vaf()
        

            
class VariantIndex:
    variants: Dict[str, Variant]
    header: List[str]
    
    def __init__(self, variant_file: str) -> None:
        print('importing variants...')
        self.variants = {}

        with open(variant_file, 'r') as infile:
           self.header = infile.readline().rstrip().split('\t')

           for line in infile:
                variant = line.rstrip().split('\t')
                variant_info = {k:v for k,v in zip(self.header, variant)}

                self.variants[variant_info['variant_id']] = Variant(variant_info) 


    def import_platform_counts(self, platform: str, platform_file: str) -> None:
        with open(platform_file, 'r') as infile:
            header = infile.readline().rstrip().split('\t')

            for line in infile:
                variant = line.rstrip().split('\t')
                variant_info = {k: v for k, v in zip(header, variant)}

                try:
                    self.variants[variant_info['id']].add_platform(platform, variant_info)
                except KeyError:
                    continue


    def multiple_test_correction(self) -> None:
        variant_ids = list([k for k, v in self.variants.items() if v.vaf_pvalue])
        unadjusted_pvalues = [self.variants[key].vaf_pvalue for key in variant_ids]
        
        adjusted_pvalues = multipletests(unadjusted_pvalues, method='fdr_bh')[1]

        for variant_id, pval in zip(variant_ids, adjusted_pvalues):
            if pval < 0.05:
                self.variants[variant_id].concordant_vaf = 'False'
            else:
                self.variants[variant_id].concordant_vaf = 'True'


    def process_variants(self, read_files: List[str]) -> None:
        for variant_id, variant in tqdm(self.variants.items()):
            variant.combine_platforms()

        if None not in read_files:
            self.multiple_test_correction()


    def print_variants(self, output_header: List[str], output_file: str, platforms: List[str]) -> None:
        with open(output_file, 'w') as outfile:
            print('\t'.join(output_header), file = outfile)

            for variant_id, variant in self.variants.items():
                if len(variant.platforms) == 0:
                    continue

                output = [variant.variant_info[x] for x in output_header if x in variant.variant_info]
                output += [variant.in_tissue, variant.phased_by, variant.haplotype, variant.concordant_vaf]
                output += [str(x) for x in (variant.is_germline, variant.is_recurrent_error)]

                for platform in platforms:
                    if platform in variant.platforms:
                        output.append(variant.platforms[platform].in_tissue) 
                        output +=  variant.platforms[platform].summarize_reads()
                    else:
                        output += ['NA'] * 21
                
                print('\t'.join(output), file = outfile)


def make_header(platform: str) -> List[str]:
    header = [
        f'{platform}_in_tissue', f'{platform}_passed_dp',f'{platform}_passed_ref', f'{platform}_passed_alt', f'{platform}_passed_ref_strands', f'{platform}_passed_alt_strands',
        f'{platform}_h1_dp', f'{platform}_h1_ref', f'{platform}_h1_alt', f'{platform}_h1_ref_strands', f'{platform}_h1_alt_strands',
        f'{platform}_h2_dp', f'{platform}_h2_ref', f'{platform}_h2_alt', f'{platform}_h2_ref_strands', f'{platform}_h2_alt_strands',
        f'{platform}_ambiguous_dp', f'{platform}_ambiguous_ref', f'{platform}_ambiguous_alt', f'{platform}_ambiguous_ref_strands', f'{platform}_ambiguous_alt_strands'
    ]

    return header


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=False, help="variant file")
    ap.add_argument("-hifi", "--hifi_counts", required=False, help="HiFi read counts")
    ap.add_argument("-ilm", "--illumina_counts", required=False, help="Illumina read counts")
    ap.add_argument("-o", "--output", required=True, help="Output")

    args = ap.parse_args()  

    variants = VariantIndex(args.variants)

    platforms = []
    output_header = variants.header + ['in_tissue', 'phased_by', 'haplotype', 
        'concordant_vaf', 'is_germline', 'is_recurrent_error']

    if args.hifi_counts:
        print('importing hifi read data...')
        output_header += make_header('hifi')
        platforms.append('hifi')
        variants.import_platform_counts('hifi', args.hifi_counts)

    if args.illumina_counts:
        print('importing illumina read data...')
        output_header += make_header('illumina')
        platforms.append('illumina')
        variants.import_platform_counts('illumina', args.illumina_counts)

    print('processing variants...')
    variants.process_variants([args.hifi_counts, args.illumina_counts])

    print('printing variants...')
    variants.print_variants(output_header, args.output, platforms)

    print('done')