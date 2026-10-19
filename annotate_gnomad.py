import hail as hl
from gnomad_toolbox.load_data import get_gnomad_release
import argparse
from tqdm import tqdm


def import_bed(bed_file):
    bed_dict = {}
    with open(bed_file, 'r') as infile:
        header = infile.readline().rstrip().split('\t')
        for line in infile:
            variant = line.rstrip().split('\t')
            bed_dict[variant[3]] = (variant[0], int(variant[2]))

    return bed_dict


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


def get_reverse_complement(allele_list):
    return [rc(allele_list[0]), rc(allele_list[1])]


def get_gnomad_data(gnomad_table, chrom, pos, allele_list):
    target_ht = ht.filter(ht.locus ==  hl.locus(chrom, pos, reference_genome='GRCh38'))
    final_info = target_ht.select(
        gene = target_ht.vep.transcript_consequences[0].gene_symbol,
        consequence = target_ht.vep.most_severe_consequence,
        af = target_ht.freq[0].AF,
        rsid = hl.delimit(target_ht.rsid, ",")
    )

    variant_data = final_info.collect()
    for variant_info in variant_data:
        if variant_info['alleles'] == allele_list:
            return variant_info
        elif variant_info['alleles'] == get_reverse_complement(allele_list):
            return variant_info

    return None

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=True, help="variant file")
    ap.add_argument("-b", "--bed", required=True, help="bed file")
    ap.add_argument("-o", "--output", required=True, help="output file")
    args = ap.parse_args()

    hg38_positions = import_bed(args.bed)

    hl.init(
        backend="local",
        master="local[16]"
    )
    ht = get_gnomad_release(data_type='genomes', version='4.1')
    print(f'''
    ''')

    with open(args.variants, 'r') as infile, open(args.output, 'w') as outfile:
        header = infile.readline().rstrip().split('\t')
        header += ['hg38_chrom', 'hg38_pos', 'rsid', 'gnomad4.1_frequency', 'gene', 'vep_most_severe']
        print('\t'.join(header), file = outfile)

        for line in tqdm(infile, total = 15873):
            variant = line.rstrip().split('\t')
            var_id = variant[2]
            alleles = [variant[3], variant[4]]
            if var_id not in hg38_positions:
                variant += ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
                print('\t'.join(variant), file = outfile)
                continue

            hg38_chrom, hg38_pos = hg38_positions[var_id]
            variant += [hg38_chrom, str(hg38_pos)]

            gnomad_info = get_gnomad_data(ht, hg38_chrom, hg38_pos, alleles)

            if not gnomad_info:
                variant += ['NA', 'NA', 'NA', 'NA']
                print('\t'.join(variant), file = outfile)
                continue

            variant += [str(gnomad_info[x]) for x in ('rsid', 'af', 'gene', 'consequence')]
            print('\t'.join(variant), file = outfile)

