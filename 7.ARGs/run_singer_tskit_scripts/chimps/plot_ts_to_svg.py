# Save as: plot_ts_to_svg.py (version for older tskit)

import tskit
import argparse
import sys
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Generate an SVG plot of a single tree from a tree sequence.")
    parser.add_argument("--ts-file", required=True, help="Path to the input tree sequence (.trees) file.")
    parser.add_argument("--svg-file", required=True, help="Path for the output SVG file.")
    args = parser.parse_args()

    try:
        if not os.path.exists(args.ts_file) or os.path.getsize(args.ts_file) == 0:
            print(f"Warning: Input file is missing or empty: {args.ts_file}", file=sys.stderr)
            with open(args.svg_file, 'w') as f:
                f.write('<svg></svg>')
            sys.exit(0)

        ts = tskit.load(args.ts_file)
        
        midpoint = ts.sequence_length / 2
        tree = ts.at(midpoint)
        
        # --- WORKAROUND FOR OLDER TSKIT ---
        if tree.num_roots > 1:
            print(f"Warning: Tree at position {midpoint:.2f} has {tree.num_roots} roots. "
                  f"Plotting the largest subtree.", file=sys.stderr)
            
            # Find the root of the largest subtree
            root_sizes = [tree.num_samples(r) for r in tree.roots]
            largest_root_index = np.argmax(root_sizes)
            target_root = tree.roots[largest_root_index]
            
            # Get all samples under that root
            samples_under_root = tree.samples(target_root)
            
            # Create a new, simplified tree sequence containing only that subtree
            subtree_ts = ts.simplify(samples=samples_under_root)
            
            # Draw the first (and only) tree from this new tree sequence
            svg_content = subtree_ts.first().draw_svg(size=(1000, 600), node_labels={})
        else:
            # Original behavior for single-root trees
            svg_content = tree.draw_svg(size=(1000, 600), node_labels={})
        
        with open(args.svg_file, "w") as f:
            f.write(svg_content)
            
        print(f"Successfully wrote SVG plot to {args.svg_file}")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()