# Save as: extract_regional_trees.py

import tskit
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Extract and trim a genomic region from a .trees file.")
    parser.add_argument("--trees-in", required=True, help="Path to input tree sequence (.trees).")
    parser.add_argument("--start", required=True, type=int, help="Start coordinate of the region.")
    parser.add_argument("--end", required=True, type=int, help="End coordinate of the region.")
    parser.add_argument("--trees-out", required=True, help="Path for the output (subsetted) .trees file.")
    args = parser.parse_args()

    try:
        if not os.path.exists(args.trees_in) or os.path.getsize(args.trees_in) == 0:
            print(f"Warning: Input file is missing or empty: {args.trees_in}. Creating empty output.", file=sys.stderr)
            open(args.trees_out, 'a').close()
            sys.exit(0)
            
        ts = tskit.load(args.trees_in)
        
        # Step 1: Extract the gappy interval
        intervals = [[args.start, args.end]]
        regional_ts_gappy = ts.keep_intervals(intervals)

        # Step 2: Use trim() to squash the sequence and re-zero coordinates
        trimmed_ts = regional_ts_gappy.trim()

        # --- DEBUGGING on the FINAL object ---
        print(f"DEBUG: Original sequence_length was: {ts.sequence_length:,.0f}", file=sys.stderr)
        print(f"DEBUG: Final trimmed sequence_length is now: {trimmed_ts.sequence_length:,.0f}", file=sys.stderr)

        if trimmed_ts.num_trees == 0:
            print(f"Warning: No trees found in the interval {args.start}-{args.end}.", file=sys.stderr)

        # Step 3: Save the final, trimmed tree sequence
        trimmed_ts.dump(args.trees_out)
        print(f"Successfully saved trimmed regional tree sequence to: {args.trees_out}", file=sys.stderr)

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()