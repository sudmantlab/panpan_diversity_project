import sys
import argparse
import csv
from Bio import SeqIO
from bisect import bisect_left

def get_auN(contig_lens, g_sum):
    """Calculates the area under the Nx curve (auN)."""
    if g_sum == 0:
        return 0.0
    sq_l_sum = sum(l * l for l in contig_lens)
    return sq_l_sum / g_sum

def get_nx_stats(contig_lens, total_sum):
    """
    Calculates Nx and Lx stats for 1-100%.
    contig_lens is sorted descending.
    """
    # Create cumulative sum list
    cumsum = 0
    nx_ratios = []
    for l in contig_lens:
        cumsum += l
        nx_ratios.append(cumsum / total_sum)

    nx_stats = {}
    for x in range(1, 101):
        x_pct = x / 100.0
        # Find the first index where cumulative sum ratio >= x_pct
        loc = bisect_left(nx_ratios, x_pct)
        
        if loc >= len(contig_lens):
            nx_stats[x] = {"n": 0, "l": len(contig_lens)}
        else:
            nx_stats[x] = {"n": contig_lens[loc], "l": loc + 1}
            
    return nx_stats

def main():
    parser = argparse.ArgumentParser(description="Calculates basic assembly stats (Python version)")
    parser.add_argument("fn_fa", help="input fasta")
    parser.add_argument("fn_out", help="output csv (tsv)")
    parser.add_argument("--genomesize", type=int, help="length for NG50, defaults to assembly length")
    parser.add_argument("--genomename", default="NA", help="name of genome for output")
    
    args = parser.parse_args()

    # Read contig lengths
    contig_lens = []
    try:
        for record in SeqIO.parse(args.fn_fa, "fasta"):
            contig_lens.append(len(record.seq))
    except Exception as e:
        print(f"Error reading fasta: {e}")
        sys.exit(1)

    if not contig_lens:
        print("No sequences found in fasta.")
        sys.exit(1)

    # Basic stats
    g_sum = sum(contig_lens)
    g_size = args.genomesize if args.genomesize else g_sum
    
    # Sort descending for Nx calculations
    contig_lens.sort(reverse=True)

    nx_stats = get_nx_stats(contig_lens, g_sum)
    ngx_stats = get_nx_stats(contig_lens, g_size)

    auN = get_auN(contig_lens, g_sum)
    auNG = get_auN(contig_lens, g_size)

    # Print to console
    print(f"Total contigs: {len(contig_lens)}")
    print(f"Total assembly length: {g_sum}")
    print(f"max contig length: {max(contig_lens)}")
    print(f"min contig length: {min(contig_lens)}")
    print(f"auN: {auN:.2f} auNG: {auNG:.2f}")

    # Mb thresholds (1, 5, 10, 50, 100 Mb)
    # Sort ascending for binary search of specific lengths
    contig_lens_asc = sorted(contig_lens)
    for mb in [1, 5, 10, 50, 100]:
        size = mb * 1_000_000
        # Find index of first element >= size
        loc = len(contig_lens_asc) - bisect_left(contig_lens_asc, size)
        print(f"\t{loc} contigs >= {mb}Mb")

    # Output to TSV
    with open(args.fn_out, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["genomeName", "stat", "v1", "v2"])
        
        writer.writerow([args.genomename, "auN", g_sum, auN])
        writer.writerow([args.genomename, "auNG", g_size, auNG])

        for x in sorted(nx_stats.keys()):
            writer.writerow([args.genomename, "Nx", x, nx_stats[x]["n"]])
        for x in sorted(nx_stats.keys()):
            writer.writerow([args.genomename, "Lx", x, nx_stats[x]["l"]])
        for x in sorted(ngx_stats.keys()):
            writer.writerow([args.genomename, "NGx", x, ngx_stats[x]["n"]])
        for x in sorted(ngx_stats.keys()):
            writer.writerow([args.genomename, "LGx", x, ngx_stats[x]["l"]])

if __name__ == "__main__":
    main()
