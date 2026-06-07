import tskit
import numpy as np
import pandas as pd
import argparse
import os
import sys
import traceback
from collections import defaultdict
import warnings

def validate_tree(tree, ts):
    """Check for common ARG generation errors"""
    issues = []
    if tree.interval.right <= tree.interval.left:
        issues.append(f"Invalid span {tree.interval}")
    return issues

def main():
    try:
        # Parse arguments
        parser = argparse.ArgumentParser(description='Compute coalescence metrics per population')
        parser.add_argument('--chrom', required=True, help='Chromosome name')
        parser.add_argument('--chrom-sizes', required=True, help='Chromosome sizes file')
        parser.add_argument('--pops', required=True, help='Population definitions CSV')
        parser.add_argument('--samples', required=True, help='Sample list text file')
        parser.add_argument('--output', required=True, help='Output CSV path')
        parser.add_argument('--num-args', type=int, default=99, help='Number of ARG samples (0-based)')
        parser.add_argument('--diag-output', default="pop_void_diagnostics", help='Diagnostic output directory')
        args = parser.parse_args()

        # Input validation
        print(f"Processing chromosome: {args.chrom}", file=sys.stderr)
        for f in [args.chrom_sizes, args.pops, args.samples]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"File missing: {f}")

        # Load chromosome size
        chrom_sizes = pd.read_csv(args.chrom_sizes, sep='\t', names=["chrom", "length"])
        chrom_row = chrom_sizes[chrom_sizes["chrom"] == args.chrom]
        if chrom_row.empty:
            raise ValueError(f"Chromosome {args.chrom} not found in {args.chrom_sizes}")
        chrom_length = chrom_row["length"].values[0]

        # Create 1kb windows
        windows = np.arange(0, chrom_length + 1000, 1000)
        windows = np.clip(windows, 0, chrom_length)
        num_windows = len(windows) - 1

        # Initialize storage
        pooled_metrics = defaultdict(lambda: {
            'pairwise_sum': 0.0,
            'pairwise_count': 0
        })
        
        pop_metrics = defaultdict(lambda: defaultdict(lambda: {
            'tmrca_sum': 0.0,
            'tmrca_weight': 0.0,
            'pairwise_sum': 0.0,
            'pairwise_count': 0
        }))

        # Load population data
        pop_df = pd.read_csv(args.pops)
        with open(args.samples) as f:
            sample_order = [line.strip() for line in f]
        sample_to_pop = pop_df.set_index("sample")["population"].to_dict()
        populations = pop_df["population"].unique()

        processed_files = 0
        total_args = args.num_args + 1

        for i in range(total_args):
            ts_path = f"singer/{args.chrom}/{args.chrom}_{i}.trees"
            if not os.path.exists(ts_path) or os.path.getsize(ts_path) == 0:
                print(f"Warning: Missing or empty {ts_path}", file=sys.stderr)
                continue

            try:
                ts = tskit.load(ts_path)
                print(f"Processing {ts_path} ({processed_files+1}/{total_args})", file=sys.stderr)
                
                # Create population-node mapping
                pop_nodes = defaultdict(list)
                for idx, sample in enumerate(sample_order):
                    if idx < ts.num_samples and sample in sample_to_pop:
                        pop = sample_to_pop[sample]
                        pop_nodes[pop].append(ts.samples()[idx])

                # Process pooled metrics (all samples)
                all_samples = list(range(ts.num_samples))
                if len(all_samples) >= 2:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        diversity = ts.diversity(sample_sets=[all_samples], windows=windows, mode="branch")
                    # Convert to pairwise times
                    avg_pairwise = diversity / 2
                    
                    for w_idx in range(num_windows):
                        # Get scalar value for this window
                        pairwise_val = avg_pairwise[w_idx].item()
                        # Filter negative values
                        if not np.isnan(pairwise_val) and pairwise_val >= 0:
                            pooled_metrics[w_idx]['pairwise_sum'] += pairwise_val
                            pooled_metrics[w_idx]['pairwise_count'] += 1

                # Precompute tree root times and positions
                tree_root_times = []
                tree_lefts = []
                tree_rights = []
                
                for tree in ts.trees():
                    issues = validate_tree(tree, ts)
                    if issues:
                        continue
                    
                    roots = tree.roots
                    root_time = np.mean([tree.time(r) for r in roots])
                    
                    tree_root_times.append(root_time)
                    tree_lefts.append(tree.interval.left)
                    tree_rights.append(tree.interval.right)
                
                # Convert to arrays for vectorized processing
                tree_root_times = np.array(tree_root_times)
                tree_lefts = np.array(tree_lefts)
                tree_rights = np.array(tree_rights)

                # Process each population
                for pop, nodes in pop_nodes.items():
                    num_nodes = len(nodes)
                    if num_nodes < 2:
                        continue
                    
                    # Pairwise diversity for population
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        diversity_pop = ts.diversity(sample_sets=[nodes], windows=windows, mode="branch")
                    avg_pairwise_pop = diversity_pop / 2
                    
                    # Process windows
                    for w_idx in range(num_windows):
                        w_start = windows[w_idx]
                        w_end = windows[w_idx+1]
                        
                        # Pairwise (within population)
                        pairwise_val = avg_pairwise_pop[w_idx].item()
                        # Filter negative values
                        if not np.isnan(pairwise_val) and pairwise_val >= 0:
                            pop_metrics[pop][w_idx]['pairwise_sum'] += pairwise_val
                            pop_metrics[pop][w_idx]['pairwise_count'] += 1
                        
                        # TMRCA
                        mask = (tree_lefts < w_end) & (tree_rights > w_start)
                        if np.any(mask):
                            overlaps = (
                                np.minimum(tree_rights[mask], w_end) -
                                np.maximum(tree_lefts[mask], w_start)
                            )
                            valid = overlaps > 0
                            if np.any(valid):
                                weighted_tmrca = float(np.sum(tree_root_times[mask][valid] * overlaps[valid]))
                                total_overlap = float(np.sum(overlaps[valid]))
                                
                                pop_metrics[pop][w_idx]['tmrca_sum'] += weighted_tmrca
                                pop_metrics[pop][w_idx]['tmrca_weight'] += total_overlap

                processed_files += 1

            except Exception as e:
                print(f"Error processing {ts_path}: {str(e)}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

        # Compile results
        results = []
        void_data = []
        for w_idx in range(num_windows):
            start = int(windows[w_idx])
            end = int(windows[w_idx+1])
            
            # Pooled metrics
            pooled_avg = (pooled_metrics[w_idx]['pairwise_sum'] / pooled_metrics[w_idx]['pairwise_count']
                          if pooled_metrics[w_idx]['pairwise_count'] > 0 else np.nan)
            
            for pop in populations:
                if pop not in pop_metrics or w_idx not in pop_metrics[pop]:
                    # Population not present in this window
                    tmrca_avg = np.nan
                    pairwise_avg = np.nan
                else:
                    metrics = pop_metrics[pop][w_idx]
                    
                    # Within-population pairwise
                    pairwise_avg = (metrics['pairwise_sum'] / metrics['pairwise_count']
                                   if metrics['pairwise_count'] > 0 else np.nan)
                    
                    # TMRCA
                    tmrca_avg = (metrics['tmrca_sum'] / metrics['tmrca_weight']
                                if metrics['tmrca_weight'] > 0 else np.nan)
                
                # Ratio calculation - only when both components are valid
                if np.isnan(pooled_avg) or np.isnan(pairwise_avg) or pairwise_avg == 0:
                    ratio = np.nan
                else:
                    ratio = pooled_avg / pairwise_avg
                
                # Record void status
                is_void = np.isnan(tmrca_avg) or np.isnan(pairwise_avg)
                void_data.append({
                    "chromosome": args.chrom,
                    "start": start,
                    "end": end,
                    "population": pop,
                    "is_void": is_void
                })
                
                # Add results with both components and ratio
                results.append({
                    "chromosome": args.chrom,
                    "start": start,
                    "end": end,
                    "population": pop,
                    "avg_tmrca": float(tmrca_avg) if not np.isnan(tmrca_avg) else np.nan,
                    "T_pooled": float(pooled_avg) if not np.isnan(pooled_avg) else np.nan,
                    "T_within": float(pairwise_avg) if not np.isnan(pairwise_avg) else np.nan,
                    "Tpooled_Twithin_ratio": float(ratio) if not np.isnan(ratio) else np.nan
                })

        # Save results
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        result_df = pd.DataFrame(results)
        result_df.to_csv(args.output, index=False)
        print(f"Processed {processed_files} ARG samples", file=sys.stderr)

        # VOID DIAGNOSTICS
        diag_dir = args.diag_output
        os.makedirs(diag_dir, exist_ok=True)
        
        # 1. Void summary statistics
        void_df = pd.DataFrame(void_data)
        void_summary = void_df.groupby("population")["is_void"].agg(
            void_count="sum", total="count"
        ).reset_index()
        void_summary["void_percent"] = void_summary["void_count"] / void_summary["total"] * 100
        
        # 2. Per-population void BED files
        for pop in populations:
            pop_voids = void_df[(void_df["population"] == pop) & void_df["is_void"]]
            if not pop_voids.empty:
                bed_path = f"{diag_dir}/{args.chrom}_{pop}_voids.bed"
                pop_voids[["chromosome", "start", "end"]].to_csv(
                    bed_path, sep='\t', index=False, header=False
                )
        
        # 3. Void positions text file
        void_positions_path = f"{diag_dir}/{args.chrom}_void_positions.txt"
        with open(void_positions_path, "w") as f:
            f.write("chromosome\tstart\tend\tpopulation\n")
            for _, row in void_df[void_df["is_void"]].iterrows():
                f.write(f"{row['chromosome']}\t{row['start']}\t{row['end']}\t{row['population']}\n")
        
        # 4. Void summary report
        summary_path = f"{diag_dir}/{args.chrom}_void_summary.txt"
        with open(summary_path, "w") as f:
            f.write(f"Chromosome: {args.chrom}\n")
            f.write(f"Total windows: {num_windows}\n")
            f.write(f"ARG samples processed: {processed_files}/{total_args}\n\n")
            f.write("Population void statistics:\n")
            f.write(void_summary.to_string(index=False))
            
        print(f"Saved void diagnostics to {diag_dir}", file=sys.stderr)

    except Exception as e:
        print(f"Fatal error: {str(e)}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
