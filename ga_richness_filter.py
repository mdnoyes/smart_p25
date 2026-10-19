import argparse
from collections import Counter
from pysam import FastaFile

def check_homopolymer(sequence):
    candidate_homopolymers = []

    last_base = sequence[0]
    homopolymer_length = 1

    for base in sequence[1:]:
        if base == last_base:
            homopolymer_length += 1
        elif base != last_base:
            candidate_homopolymers.append((homopolymer_length, last_base))
            homopolymer_length = 1
        last_base = base

    candidate_homopolymers.append((homopolymer_length, last_base))

    homopolymers = [f'{x[0]}{x[1]}' for x in candidate_homopolymers if x[0] > 4]

    if homopolymers:    
        return (';'.join(homopolymers))
    else:
        return 'NA'


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


def calculate_percent(base_counts, sequence_length, base_1, base_2):
    if base_1 in base_counts and base_2 in base_counts:
        return (base_counts[base_1] + base_counts[base_2]) / sequence_length
    else:
        return 0


def get_base_content(base_counts, sequence_length, base_1, base_2):
    forward_percent = calculate_percent(base_counts, sequence_length, base_1, base_2)
    reverse_percent = calculate_percent(base_counts, sequence_length, rc(base_1), rc(base_2))

    if forward_percent >= reverse_percent:
        return forward_percent, 'forward'
    else:
        return reverse_percent, 'reverse'



def check_sequence_context(chrom, pos, fasta):
    local_sequence = fasta.fetch(chrom, pos-6, pos+5).upper()
    sequence_length = len(local_sequence)
    base_counts = Counter(local_sequence)

    ga_content, ga_direction = get_base_content(base_counts, sequence_length, 'G', 'A')
    gc_content, gc_direction = get_base_content(base_counts, sequence_length, 'G', 'C')
    homopolymers = check_homopolymer(local_sequence)

    return ga_content, ga_direction, gc_content, homopolymers, local_sequence


def new_homopolymer_filter(homopolymers, ref, alt, local_sequence):
    if homopolymers == 'NA':
        return True

    if ref != local_sequence[5] and alt != local_sequence[5]:
        if rc(ref) != local_sequence[5] and rc(alt) != local_sequence[5]:
            raise ValueError('Your sequence may be incorrectly defined! Check indexing.')

    left_flank = local_sequence[0:5] 
    right_flank = local_sequence[6:11]

    homopolymers = homopolymers.split(';')
    for h in homopolymers:
        count = int(h[:-1]) 
        base = h[-1]
        homo_seq = base * count

        is_adjacent = (left_flank == homo_seq) or (right_flank == homo_seq)
        flank_len = min(count - 1, 5)

        left  = left_flank[-flank_len:]
        right = right_flank[:flank_len]
        target = base * (flank_len + 1)

        part_of_five_prime = (left + ref == target) or (left + alt == target)
        part_of_three_prime = (ref + right == target) or (alt + right == target)


        if is_adjacent and base in (ref, alt):
            return False
        if part_of_five_prime or part_of_three_prime:
            return False
        
    return True
        



def ga_filter(ga_content, ga_direction, ref, alt, supported_by_illumina, tissues_with_variant):
    # I don't think I'm worried about transitions anymore, but leaving this in case I change my mind
    # if ga_direction == 'forward':
    #     transitions = {('G', 'A'), ('A', 'G')}
    # elif ga_direction == 'reverse':
    #     transitions = {('C', 'T'), ('T', 'C')}

    # Illumina shouldn't have the same GA-richness problems as HiFi, so if we see Illumina support, it's probably ok.
    if ga_content > 0.7 and supported_by_illumina == 'False':
        # if shared_tissues in ('NA', 'fibroblast') and (ref, alt) in transitions:
        #         return False
        if len(tissues_with_variant) == 1:
            return False

    # In really GA-rich regions, the probability of these errors being incorporated into the DSA
    # is higher, so let's be really strict here. 
    # if ga_content > 0.9 and (ref, alt) in transitions:
    # if ga_content > 0.9:
    #     return False
       
    return True


def apply_filters(variant_dict, fasta):
    chrom, pos, ref, alt = [variant_dict[x] for x in ('chrom', 'pos', 'ref', 'alt')]
    # chrom, pos, ref, alt = var_dict['id'].split('_')
    ga, ga_direction, gc, homopolymers, sequence = check_sequence_context(chrom, int(pos), fasta)

    variant_tissues = var_dict['tissues_with_variant'].split(';')
    supported_by_illumina = var_dict['supported_in_illumina']

    # total_illumina_read_support = 0

    # for key, value in var_dict.items():
    #     if 'illumina' not in key:
    #         continue
    #     if key.split('_')[0] not in variant_tissues:
    #         continue

    #     read_counts = value.split('|')
    #     total_illumina_read_support += sum([int(x.split(';')[2]) for x in read_counts])



    # Filter for GA-richness and homopolymers
    # No filter yet for GC-richness, but we can revisit if necessary
    pass_ga_filter = ga_filter(ga, ga_direction, ref, alt, supported_by_illumina, variant_tissues)
    pass_homopolymer_filter = new_homopolymer_filter(homopolymers, ref, alt, sequence)

    if pass_ga_filter and pass_homopolymer_filter:
        return ga, gc, homopolymers, True
    else:
        return ga, gc, homopolymers, False

    

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=True, help="variant file")
    ap.add_argument("-f", "--fasta", required=True, help="fasta file")
    ap.add_argument("-o", "--output", required=True, help="output file")
    args = ap.parse_args()

    fasta = FastaFile(args.fasta)

    with open(args.variants, 'r') as infile, open(args.output, 'w') as outfile:
        header = infile.readline().rstrip().split('\t')
        header += ['ga_content_10bp', 'gc_content_10bp', 'homopolymers_10bp']
        print('\t'.join(header), file = outfile)

        for line in infile:
            variant = line.rstrip().split('\t')
            var_dict = {k:v for k,v in zip(header, variant)}
            ga, gc, homopolymers, pass_filters = apply_filters(var_dict, fasta)

            if pass_filters:
                variant += [str(x) for x in (ga, gc, homopolymers)]
                print('\t'.join(variant), file = outfile)

    print('done')