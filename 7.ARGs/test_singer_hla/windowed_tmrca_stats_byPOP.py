import tskit
import numpy as np
import pandas as pd
from collections import defaultdict
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compute population-specific coalescence metrics from ARG samples.')
    parser.add_argument('--prefix', type=str, required=True,
                        help='Filename prefix for ARG samples (e.g., "hprc_ht2t_chr6_")')
    parser.add_argument('--pops', type=str, required=True,
                        help='Path to populations CSV file')
    parser.add_argument('--samples', type=str, required=True,
                        help='Path to samples text file')
    parser.add_argument('--output', type=str, required=True,
                        help='Output CSV filename')
    parser.add_argument('--offset', type=int, default=26390105,
                        help='T2T start coordinate (default: 26390105)')
    parser.add_argument('--chrom-end', type=int, default=33135780,
                        help='T2T end coordinate (default: 33135780)')
    parser.add_argument('--num-args', type=int, default=99,
                        help='Number of ARG samples (default: 99)')
    
    args = parser.parse_args()

    # Load population data and sample list
    pop_df = pd.read_csv(args.pops)
    with open(args.samples, "r") as f:
        sample_order = [line.strip() for line in f]
    sample_to_pop = pop_df.set_index("sample")["population"].to_dict()

    # Define fixed 1kb windows in T2T coordinates
    t2t_windows = np.arange(args.offset, args.chrom_end + 1000, 1000)
    t2t_windows = np.clip(t2t_windows, args.offset, args.chrom_end)
    num_windows = len(t2t_windows) - 1

    # Initialize storage
    window_data = defaultdict(lambda: {
        "diversity_sum": np.zeros(num_windows),
        "diversity_count": np.zeros(num_windows),
        "tmrca_sum": np.zeros(num_windows),
        "tmrca_count": np.zeros(num_windows),
    })

    # Process each ARG
    for i in range(args.num_args):
        ts = tskit.load(f"{args.prefix}{i}.trees")
        
        # Map sample names to node IDs
        samples_in_ts = ts.samples()
        sample_node_map = {name: node_id for node_id, name in zip(samples_in_ts, sample_order)}
        
        # Group nodes by population
        pop_nodes = defaultdict(list)
        for name, node_id in sample_node_map.items():
            if name in sample_to_pop:
                pop = sample_to_pop[name]
                pop_nodes[pop].append(node_id)
        
        # Convert windows to local coordinates
        local_windows = (t2t_windows - args.offset).clip(0, ts.sequence_length)
        
        # Process each population
        for pop, nodes in pop_nodes.items():
            if len(nodes) < 2:
                continue
            
            # Pairwise Diversity
            diversity = ts.diversity(sample_sets=[nodes], windows=local_windows, mode="branch")
            avg_pairwise_time = diversity / 2
            
            for w_idx in range(num_windows):
                if not np.isnan(avg_pairwise_time[w_idx]):
                    window_data[pop]["diversity_sum"][w_idx] += avg_pairwise_time[w_idx]
                    window_data[pop]["diversity_count"][w_idx] += 1

            # Compute Windowed TMRCA
            tree_tmrcas = []
            tree_lefts = []
            tree_rights = []
            for tree in ts.trees():
                root = tree.root
                if any(tree.is_descendant(node, root) for node in nodes):
                    tree_tmrcas.append(tree.time(root))
                    tree_lefts.append(tree.interval.left)
                    tree_rights.append(tree.interval.right)
            
            tree_tmrcas = np.array(tree_tmrcas)
            tree_lefts = np.array(tree_lefts)
            tree_rights = np.array(tree_rights)
            
            for w_idx in range(num_windows):
                w_start = local_windows[w_idx]
                w_end = local_windows[w_idx + 1]
                
                mask = (tree_lefts < w_end) & (tree_rights > w_start)
                if np.any(mask):
                    overlap_lengths = (
                        np.minimum(tree_rights[mask], w_end) 
                        - np.maximum(tree_lefts[mask], w_start)
                    )
                    weighted_tmrca = np.sum(tree_tmrcas[mask] * overlap_lengths)
                    total_overlap = np.sum(overlap_lengths)
                    
                    window_data[pop]["tmrca_sum"][w_idx] += weighted_tmrca
                    window_data[pop]["tmrca_count"][w_idx] += total_overlap

    # Compile results
    results = []
    for pop, data in window_data.items():
        for w_idx in range(num_windows):
            if data["diversity_count"][w_idx] > 0 or data["tmrca_count"][w_idx] > 0:
                diversity_avg = (
                    data["diversity_sum"][w_idx] / data["diversity_count"][w_idx] 
                    if data["diversity_count"][w_idx] > 0 else np.nan
                )
                tmrca_avg = (
                    data["tmrca_sum"][w_idx] / data["tmrca_count"][w_idx] 
                    if data["tmrca_count"][w_idx] > 0 else np.nan
                )
                results.append({
                    "population": pop,
                    "start": t2t_windows[w_idx],
                    "end": t2t_windows[w_idx + 1],
                    "avg_pairwise_coalescence_time": diversity_avg,
                    "avg_tmrca": tmrca_avg
                })

    pd.DataFrame(results).to_csv(args.output, index=False)

if __name__ == "__main__":
    main()