import argparse
import gzip
from pybedtools import BedTool

def format_data(variant_line, tissue_name):
    info_index = {value: index for index, value in enumerate(variant_line[8].split(':'))}
    var_info = variant_line[9].split(':')

    var_data = { f'{tissue_name}_germline': str(variant_line[6] == 'GERMLINE'), 
            f'{tissue_name}_qual': variant_line[5]}

    for item in ('GQ', 'DP', 'AD', 'VAF'):
        var_data[f'{tissue_name}_{item}'] = var_info[info_index[item]]

    return var_data



def import_variants(tissue_name, vcf, var_dict):
    with gzip.open(vcf, 'rt') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            else:
                variant = line.rstrip().split('\t')
                if variant[6] not in ('GERMLINE', 'PASS'):
                    continue
                genotype = variant[9].split(':')[0]
                if genotype in ('./.'):
                    continue
                ref, alt = variant[3:5]
                if len(ref) != 1 or len(alt) != 1:
                    continue

                chrom, pos = variant[:2]
                qual = variant[5]
                info = variant[9]

                id = f'{chrom}_{pos}_{ref}_{alt}'
                if id not in var_dict:
                    var_dict[id] = {'bed': (chrom, int(pos)-1, int(pos))}
                
                var_dict[id][tissue_name] = format_data(variant, tissue_name)
                if tissue_name == "merged_tissues":
                    var_dict[id][tissue_name]["annotation"] = variant[7]

    return var_dict


def find_variants_in_region(var_dict, bed_file):
    variant_bed = BedTool([v['bed'] for k,v in variant_dict.items()])
    flagged_bed = BedTool(bed_file)

    intersected_variants = variant_bed.intersect(flagged_bed)
    variants_in_region = intersected_variants.to_dataframe()

    return list(variants_in_region.itertuples(index=False, name=None))



def flag_variants(var_dict, bed_list):
    flagged_variants = []

    for bed_file in bed_list:
        flagged_variants += find_variants_in_region(var_dict, bed_file)
        print(len(flagged_variants))


    return set(flagged_variants)


#def filter_callsets(var_dict, flagged):
def filter_callsets(var_dict):
    filtered_calls = []

    for id, info in var_dict.items():
        #if info['bed'] in flagged:
        #    continue
        
        variant_data = dict(zip(['chrom', 'pos', 'ref', 'alt'], id.split('_')))
        variant_data['id'] = id

        callsets = [x for x in info.keys() if x != 'bed']
        variant_data['callsets'] = '_'.join(callsets)

        for callset in callsets:
            variant_data.update(info[callset])

        
        filtered_calls.append(variant_data)

    return filtered_calls



def print_callset(final_calls, output_file):
    output_header = ('chrom', 'pos', 'id', 'ref', 'alt', 'callsets', 'annotation',
            'fibroblast_germline',  'fibroblast_qual', 'fibroblast_GQ', 'fibroblast_DP', 'fibroblast_AD', 'fibroblast_VAF',
            'merged_tissues_germline','merged_tissues_qual', 'merged_tissues_GQ',  'merged_tissues_DP', 'merged_tissues_AD', 'merged_tissues_VAF')

    with open(output_file, 'w') as outfile:
        print('\t'.join(output_header), file = outfile)

        for variant in final_calls:
            for column in output_header:
                variant.setdefault(column, 'NA')

            out_string = '\t'.join([variant[column] for column in output_header])
            print(out_string, file = outfile)



if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fibroblast", required=True, help="fibroblast vcf")
    ap.add_argument("-t", "--tissue", required=True, help="tissue vcf")
#    ap.add_argument("-b", "--exclusion_beds", required=True, help="flagger bed", nargs = "+")
    ap.add_argument("-o", "--output", required=True, help="output")
    args = ap.parse_args()

    variant_dict = {}

    variant_dict = import_variants('fibroblast', args.fibroblast, variant_dict)
    variant_dict = import_variants('merged_tissues', args.tissue, variant_dict)
    print(len(variant_dict))

#    variants_in_excluded_regions = flag_variants(variant_dict, args.exclusion_beds) 
#    print(len(variants_in_excluded_regions))

#    merged_callset = filter_callsets(variant_dict, variants_in_excluded_regions)
    merged_callset = filter_callsets(variant_dict)
    print(len(merged_callset))

    print_callset(merged_callset, args.output)

