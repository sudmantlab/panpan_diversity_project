# Save as: ts_to_newick.py

import tskit
import argparse
import sys
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Extracts a single tree at a specific absolute coordinate from a trimmed regional tree sequence.")
    parser.add_argument("--ts-file", required=True, help="Path to the input regional (trimmed) tree sequence (.trees) file.")
    parser.add_argument("--newick-file", required=True, help="Path for the output Newick (.nwk) file.")
    # --- THESE ARGUMENTS MUST BE PRESENT ---
    parser.add_argument("--absolute-position", type=float, required=True, help="Absolute genomic coordinate to extract the tree from.")
    parser.add_argument("--region-start", type=float, required=True, help="The original start coordinate of the region before trimming.")
    # ------------------------------------
    args = parser.parse_args()

    try:
        if not os.path.exists(args.ts_file) or os.path.getsize(args.ts_file) == 0:
            print(f"Warning: Input file is missing or empty: {args.ts_file}. Cannot generate Newick tree.", file=sys.stderr)
            with open(args.newick_file, 'w') as f:
                pass # Create empty file
            sys.exit(0)

        # Load the trimmed regional tree sequence
        ts = tskit.load(args.ts_file)

        # Calculate the position relative to the start of the trimmed sequence
        relative_position = args.absolute_position - args.region_start

        # Check if the relative position is valid within the trimmed sequence
        if not (0 <= relative_position < ts.sequence_length):
            print(f"Error: Calculated relative position {relative_position} "
                  f"(from absolute {args.absolute_position}) is outside the trimmed sequence bounds "
                  f"(0-{ts.sequence_length}). Check your peak coordinate.", file=sys.stderr)
            sys.exit(1)

        print(f"Using absolute position {args.absolute_position}, which corresponds to "
              f"relative position {relative_position:.2f} in the trimmed sequence.", file=sys.stderr)

        # Get the tree at the calculated relative position
        tree = ts.at(relative_position)

        # Handle multiple roots
        target_root = None # Define target_root
        if tree.num_roots > 1:
            print(f"Warning: Tree at relative position {relative_position:.2f} has {tree.num_roots} roots. Writing the largest subtree.", file=sys.stderr)
            root_sizes = [tree.num_samples(r) for r in tree.roots]
            largest_root_index = np.argmax(root_sizes)
            target_root = tree.roots[largest_root_index]
            newick_string = tree.newick(root=target_root)
        else:
            newick_string = tree.newick() # No root needed if single root

        with open(args.newick_file, "w") as f:
            f.write(newick_string)

        print(f"Successfully wrote Newick tree from relative position {relative_position:.2f} to {args.newick_file}")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()