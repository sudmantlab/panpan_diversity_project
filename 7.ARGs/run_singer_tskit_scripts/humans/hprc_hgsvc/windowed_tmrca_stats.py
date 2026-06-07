import tskit
import numpy as np
import pandas as pd
import argparse
import os
import sys
import traceback

def validate_tree(tree, ts):
    """Check for common ARG generation errors"""
    issues = []
    if tree.interval.right <= tree.interval.left:
        issues.append(f"Invalid span {tree.interval}")
    if len(tree.roots) == ts.num_samples:
        issues.append(f"Full non-coalescence ({len(tree.roots)} roots)")
    return issues

def main():
    try:
        # Parse arguments with validation
        parser = argparse.ArgumentParser(description='Compute genome-wide coalescence metrics')
        parser.add_argument('--chrom', required=True, help='Chromosome name')
        parser.add_argument('--chrom-sizes', required=True, help='Chromosome sizes file')
        parser.add_argument('--output', required=True, help='Output CSV path')
        parser.add_argument('--num-args', type=int, default=99, help='Number of ARG samples (0-based)')
        parser.add_argument('--diag-output', default="void_diagnostics", help='Diagnostic output directory')
        args = parser.parse_args()

        # Input validation
        if not os.path.exists(args.chrom_sizes):
            raise FileNotFoundError(f"Chrom sizes file missing: {args.chrom_sizes}")
            
        if not os.path.exists(f"singer/{args.chrom}"):
            raise FileNotFoundError(f"Chromosome directory missing: singer/{args.chrom}")

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

        # Initialize results with sum/count for accurate averaging
        results = []
        for w_idx in range(num_windows):
            results.append({
                "chromosome": args.chrom,
                "start": int(windows[w_idx]),
                "end": int(windows[w_idx+1]),
                "tmrca_sum": 0.0,
                "tmrca_count": 0,
                "pairwise_sum": 0.0,
                "pairwise_count": 0
            })

        processed_files = 0

        # Process each ARG sample
        for i in range(args.num_args + 1):
            ts_path = f"singer/{args.chrom}/{args.chrom}_{i}.trees"
            if not os.path.exists(ts_path):
                print(f"Warning: Missing ARG file {ts_path}", file=sys.stderr)
                continue

            try:
                ts = tskit.load(ts_path)
                
                # Skip if insufficient samples
                if ts.num_samples < 2:
                    print(f"Skipping {ts_path} - only {ts.num_samples} samples", file=sys.stderr)
                    continue

                # Process trees with validation
                tree_tmrcas, tree_lefts, tree_rights = [], [], []
                invalid_trees = 0
                
                for tree in ts.trees():
                    issues = validate_tree(tree, ts)
                    if issues:
                        print(f"Invalid tree in {ts_path}: {'; '.join(issues)}", file=sys.stderr)
                        invalid_trees += 1
                        continue
                    
                    # Calculate average root time
                    roots = tree.roots
                    root_time = np.mean([tree.time(r) for r in roots])
                    
                    tree_tmrcas.append(root_time)
                    tree_lefts.append(tree.interval.left)
                    tree_rights.append(tree.interval.right)

                # Convert to arrays for vector operations
                tree_tmrcas = np.array(tree_tmrcas)
                tree_lefts = np.array(tree_lefts)
                tree_rights = np.array(tree_rights)

                # Calculate pairwise times
                diversity = ts.diversity(windows=windows, mode="branch")
                avg_pairwise_time = diversity / 2

                # Window processing
                for w_idx in range(num_windows):
                    w_start = windows[w_idx]
                    w_end = windows[w_idx+1]

                    # TMRCA calculation
                    mask = (tree_lefts < w_end) & (tree_rights > w_start)
                    if np.any(mask):
                        overlaps = (
                            np.minimum(tree_rights[mask], w_end) -
                            np.maximum(tree_lefts[mask], w_start)
                        )
                        valid_overlaps = overlaps > 0
                        
                        if np.any(valid_overlaps):
                            weighted_tmrca = np.sum(tree_tmrcas[mask][valid_overlaps] * overlaps[valid_overlaps])
                            total_overlap = np.sum(overlaps[valid_overlaps])
                            current_tmrca = weighted_tmrca / total_overlap
                            
                            results[w_idx]["tmrca_sum"] += current_tmrca
                            results[w_idx]["tmrca_count"] += 1

                    # Pairwise time calculation
                    val = avg_pairwise_time[w_idx]
                    if not np.isnan(val) and val >= 0:
                        results[w_idx]["pairwise_sum"] += val
                        results[w_idx]["pairwise_count"] += 1

                processed_files += 1
                print(f"Processed {ts_path} ({processed_files}/{args.num_args+1})", file=sys.stderr)

            except Exception as e:
                print(f"Error processing {ts_path}: {str(e)}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

        # Final calculations and cleanup
        final_results = []
        void_windows = []
        for res in results:
            tmrca_val = (res["tmrca_sum"] / res["tmrca_count"] if res["tmrca_count"] > 0 else np.nan)
            pairwise_val = (res["pairwise_sum"] / res["pairwise_count"] if res["pairwise_count"] > 0 else np.nan)
            
            final_entry = {
                "chromosome": res["chromosome"],
                "start": res["start"],
                "end": res["end"],
                "avg_tmrca": tmrca_val,
                "avg_pairwise_coalescence_time": pairwise_val
            }
            final_results.append(final_entry)
            
            # Record void windows
            if np.isnan(tmrca_val) or np.isnan(pairwise_val):
                void_windows.append({
                    "chromosome": res["chromosome"],
                    "start": res["start"],
                    "end": res["end"]
                })

        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        # Save results
        result_df = pd.DataFrame(final_results)
        result_df.to_csv(args.output, index=False)
        print(f"Successfully processed {processed_files} ARG samples for {args.chrom}", file=sys.stderr)

        # ===== DIAGNOSTIC OUTPUTS (NO MATPLOTLIB) =====
        diag_dir = args.diag_output
        os.makedirs(diag_dir, exist_ok=True)
        
        # 1. Chromosome-level void summary
        void_count = len(void_windows)
        void_percent = void_count / len(results) * 100 if len(results) > 0 else 0
        
        with open(f"{diag_dir}/{args.chrom}_void_summary.txt", "w") as f:
            f.write(f"Chromosome: {args.chrom}\n")
            f.write(f"Total windows: {len(results)}\n")
            f.write(f"Void windows: {void_count} ({void_percent:.2f}%)\n")
            
        # 2. BED file of void regions
        if void_windows:
            void_df = pd.DataFrame(void_windows)
            void_df.to_csv(f"{diag_dir}/{args.chrom}_void_regions.bed", 
                          sep="\t", index=False, header=False)
        
        # 3. Void positions text file
        with open(f"{diag_dir}/{args.chrom}_void_positions.txt", "w") as f:
            f.write("chromosome\tstart\tend\n")
            for window in void_windows:
                f.write(f"{window['chromosome']}\t{window['start']}\t{window['end']}\n")

    except Exception as e:
        print(f"FATAL ERROR: {str(e)}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()