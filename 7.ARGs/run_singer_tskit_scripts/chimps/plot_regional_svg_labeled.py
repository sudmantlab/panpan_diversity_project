# Save as: plot_regional_svg_labeled.py

import tskit
import argparse
import sys
import os
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Generate a labeled SVG plot of a tree at a specific coordinate +/- offset.")
    parser.add_argument("--ts-file", required=True, help="Path to input regional (trimmed) tree sequence.")
    parser.add_argument("--svg-file", required=True, help="Path for the output SVG file.")
    parser.add_argument("--samples-file", required=True, help="Path to samples.txt.")
    parser.add_argument("--pop-file", required=True, help="Path to population CSV.")
    parser.add_argument("--absolute-position", type=float, required=True, help="Absolute genomic coordinate of the central peak.")
    parser.add_argument("--region-start", type=float, required=True, help="Original start coordinate of the region.")
    # --- THIS ARGUMENT MUST BE PRESENT ---
    parser.add_argument("--offset", type=int, default=0, help="Offset (bp) from the absolute position (e.g., -5, 0, 5).")
    # ------------------------------------
    args = parser.parse_args()

    # --- Load Sample Names & Pop Data ---
    try:
        with open(args.samples_file) as f: sample_names = [ln.strip() for ln in f if ln.strip()]
        num_diploid_samples = len(sample_names)
        num_expected_tips = num_diploid_samples * 2
        print(f"Loaded {num_diploid_samples} sample names.", file=sys.stderr)
        pop_df = pd.read_csv(args.pop_file)
        sample_to_pop = pop_df.set_index("sample")["population"].to_dict()
        print(f"Loaded pop data for {len(sample_to_pop)} samples.", file=sys.stderr)
    except Exception as e:
        print(f"Error loading samples/pop file: {e}", file=sys.stderr); sys.exit(1)

    # --- Load Tree Sequence ---
    try:
        if not os.path.exists(args.ts_file) or os.path.getsize(args.ts_file) == 0:
            print(f"Warning: Input empty: {args.ts_file}", file=sys.stderr)
            with open(args.svg_file, 'w') as f: f.write('<svg></svg>'); sys.exit(0)
        ts = tskit.load(args.ts_file)
        print(f"Loaded trimmed ts length {ts.sequence_length:.0f} bp.", file=sys.stderr)
    except Exception as e:
        print(f"Error loading ts: {e}", file=sys.stderr); sys.exit(1)

    # --- Create Label Map ---
    node_labels = {}
    for i, sample_name in enumerate(sample_names):
        population = sample_to_pop.get(sample_name, 'Unknown')
        node_labels[2 * i] = f"{sample_name}_hap1_{population}"
        node_labels[2 * i + 1] = f"{sample_name}_hap2_{population}"

    # --- Select Tree at Specific Coordinate + Offset ---
    tree_to_draw = None
    try:
        target_absolute_pos = args.absolute_position + args.offset
        relative_position = target_absolute_pos - args.region_start

        if not (0 <= relative_position < ts.sequence_length):
            print(f"Error: Target relative position {relative_position:.2f} (abs={target_absolute_pos}, offset={args.offset}) is out of bounds (0-{ts.sequence_length}).", file=sys.stderr)
            with open(args.svg_file, 'w') as f: f.write(f'<svg><text x="10" y="50">Error: Target position {target_absolute_pos:.0f} out of bounds.</text></svg>')
            sys.exit(0) # Exit gracefully

        print(f"Extracting tree at abs pos {target_absolute_pos:.2f} (rel: {relative_position:.2f}, offset: {args.offset})", file=sys.stderr)
        original_tree = ts.at(relative_position)

        # Handle multi-root using simplify
        if original_tree.num_roots > 1:
            print(f"Warning: Tree has {original_tree.num_roots} roots. Plotting largest.", file=sys.stderr)
            root_sizes = [original_tree.num_samples(r) for r in original_tree.roots]
            target_root = original_tree.roots[np.argmax(root_sizes)]
            samples_under_root = original_tree.samples(target_root)
            subtree_ts = ts.simplify(samples=samples_under_root, map_nodes=False)
            if subtree_ts.num_trees > 0:
                 tree_to_draw = subtree_ts.first()
                 new_node_labels = {nid: lbl for nid, lbl in node_labels.items() if nid in samples_under_root and nid < subtree_ts.num_nodes}
                 node_labels = new_node_labels
            else: print(f"Warning: Could not create subtree.", file=sys.stderr)
        elif original_tree.num_roots == 1:
            tree_to_draw = original_tree
        else: print(f"Warning: Tree has no roots.", file=sys.stderr)

        if tree_to_draw is None:
             print("Error: No valid tree object available for drawing.", file=sys.stderr)
             with open(args.svg_file, 'w') as f: f.write('<svg><text x="50" y="50">Error: No valid tree.</text></svg>')
             sys.exit(0)

    except Exception as e:
        print(f"Error selecting/simplifying tree: {e}", file=sys.stderr); sys.exit(1)

    # --- Generate SVG with Customizations ---
    try:
        style = """
            .node > .lab { font-size: 10px; }
            .leaf > .lab { font-size: 6px; text-anchor: start; transform: translateX(5px); }
        """
        svg_content = tree_to_draw.draw_svg(
            size=(3000, 1000), # Wider, less tall
            node_labels=node_labels,
            y_axis=True, y_label="Time (Generations Ago)",
            x_axis=True, style=style
        )
        with open(args.svg_file, "w") as f: f.write(svg_content)
        print(f"Successfully wrote labeled SVG plot to {args.svg_file}")
    except Exception as e:
        print(f"Error generating SVG: {e}", file=sys.stderr); sys.exit(1)

# --- Standard execution block ---
if __name__ == "__main__":
    main()