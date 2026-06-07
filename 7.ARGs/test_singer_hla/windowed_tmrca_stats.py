import tskit
import numpy as np
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compute coalescence metrics from ARG samples.')
    parser.add_argument('--prefix', type=str, required=True,
                        help='Filename prefix for ARG samples (e.g., "hprc_ht2t_chr6_")')
    parser.add_argument('--offset', type=int, default=26390105,
                        help='T2T start coordinate (default: 26390105)')
    parser.add_argument('--chrom-end', type=int, default=33135780,
                        help='T2T end coordinate (default: 33135780)')
    parser.add_argument('--num-args', type=int, default=99,
                        help='Number of ARG samples (default: 99)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output CSV filename')
    
    args = parser.parse_args()

    # Define fixed 1kb windows in T2T coordinates
    t2t_windows = np.arange(args.offset, args.chrom_end + 1000, 1000)
    t2t_windows = np.clip(t2t_windows, args.offset, args.chrom_end)
    num_windows = len(t2t_windows) - 1

    # Initialize window data storage
    window_data = {
        "start": t2t_windows[:-1],
        "end": t2t_windows[1:],
        "tmrca_sum": np.zeros(num_windows),
        "tmrca_count": np.zeros(num_windows),
        "pairwise_sum": np.zeros(num_windows),
        "pairwise_count": np.zeros(num_windows),
    }

    # Process each ARG
    for i in range(args.num_args):
        ts = tskit.load(f"{args.prefix}{i}.trees")
        
        # Convert T2T windows to local coordinates for this ARG
        local_windows = t2t_windows - args.offset
        local_windows = np.clip(local_windows, 0, ts.sequence_length)

        # Compute TMRCA metrics
        tree_tmrcas = np.array([tree.time(tree.root) for tree in ts.trees()])
        tree_lefts = np.array([tree.interval.left for tree in ts.trees()])
        tree_rights = np.array([tree.interval.right for tree in ts.trees()])
        
        for window_idx in range(num_windows):
            w_start = local_windows[window_idx]
            w_end = local_windows[window_idx + 1]
            
            overlaps = (tree_lefts < w_end) & (tree_rights > w_start)
            if np.any(overlaps):
                overlap_lengths = (
                    np.minimum(tree_rights[overlaps], w_end) 
                    - np.maximum(tree_lefts[overlaps], w_start)
                )
                avg_tmrca = np.average(tree_tmrcas[overlaps], weights=overlap_lengths)
                window_data["tmrca_sum"][window_idx] += avg_tmrca
                window_data["tmrca_count"][window_idx] += 1

        # Compute pairwise coalescence times
        diversity = ts.diversity(
            sample_sets=[ts.samples()],
            windows=local_windows,
            mode="branch"
        )
        avg_pairwise_time = diversity / 2
        
        for window_idx in range(num_windows):
            if not np.isnan(avg_pairwise_time[window_idx]):
                window_data["pairwise_sum"][window_idx] += avg_pairwise_time[window_idx]
                window_data["pairwise_count"][window_idx] += 1

    # Calculate final averages
    df_avg = pd.DataFrame({
        "start": window_data["start"],
        "end": window_data["end"],
        "avg_tmrca": window_data["tmrca_sum"] / window_data["tmrca_count"],
        "avg_pairwise_coalescence_time": window_data["pairwise_sum"] / window_data["pairwise_count"]
    })

    df_avg.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()