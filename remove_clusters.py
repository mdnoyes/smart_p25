import argparse
from collections import defaultdict


def import_clusters(cluster_file):
    clusters = defaultdict(set)

    with open(cluster_file, 'r') as infile:
        header = infile.readline().rstrip().split('\t')

        for line in infile:
            pair = line.rstrip().split('\t')
            cluster_id = int(pair[-1])

            clusters[cluster_id].add(pair[0])
            clusters[cluster_id].add(pair[1])

    return clusters


def summarize_clusters(cluster_file):
    clusters = import_clusters(cluster_file)

    blacklist = set()
    whitelist = {}

    for clustered_variants in clusters.values():
        size = len(clustered_variants)
        clustered_ids = ";".join(clustered_variants)
        
        if size >= 4:
            blacklist.update(clustered_variants)
        else:
            for variant in clustered_variants:
                whitelist[variant] = ['True', str(size), clustered_ids]

    return blacklist, whitelist

    
def filter_variants(variant_file, blacklist, whitelist, output_file):
    with open(variant_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = infile.readline().rstrip().split('\t')
        header_idx = {k:v for v,k in enumerate(header)}
        header += ['in_cluster', 'cluster_size', 'clustered_variants']
        print('\t'.join(header), file = outfile)

        for line in infile:
            variant_info = line.rstrip().split('\t')
            variant_id = variant_info[header_idx['id']]

            if variant_id in blacklist:
                continue

            output_info = variant_info + whitelist.get(variant_id, ['False', '0', 'NA'])

            print('\t'.join(output_info), file = outfile)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--variants", required=True, help="variant file")
    ap.add_argument("-c", "--clusters", required=True, help="cluster file")
    ap.add_argument("-o", "--output", required=True, help="output file")
    args = ap.parse_args()

    failing_clusters, passing_clusters = summarize_clusters(args.clusters)
    filter_variants(args.variants, failing_clusters, passing_clusters, args.output)

    print('done')
