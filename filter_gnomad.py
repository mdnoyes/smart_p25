import argparse


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=True, help="variant file")
    ap.add_argument("-o", "--output", required=True, help="output file")
    args = ap.parse_args()

    with open(args.variants, 'r') as infile, open(args.output, 'w') as outfile:
        header = infile.readline().rstrip().split('\t')
        header_idx = {k:v for v,k in enumerate(header)}
        print('\t'.join(header), file = outfile)

        for line in infile:
            variant = line.rstrip().split('\t')
            gnomad_frequency = variant[header_idx['gnomad4.1_frequency']]
            if gnomad_frequency != 'NA' and gnomad_frequency != 'None':
                gnomad_frequency = float(gnomad_frequency)
                if gnomad_frequency >= 10 ** -4:
                    continue
                
            print('\t'.join(variant), file = outfile)


    print('done')