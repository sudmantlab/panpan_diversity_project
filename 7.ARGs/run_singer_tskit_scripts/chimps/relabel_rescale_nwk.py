# Save as: relabel_rescale_nwk.py

import dendropy
import argparse
import sys
import os
import re
import pandas as pd # Need pandas to read the population CSV

def main():
    parser = argparse.ArgumentParser(description="Relabel tips (haplotype & population aware, 1-based input) and rescale branch lengths.")
    parser.add_argument("--nwk-in", required=True, help="Path to the input Newick (.nwk) file (expected labels 1-48).")
    parser.add_argument("--nwk-out", required=True, help="Path for the output modified Newick (.nwk) file.")
    parser.add_argument("--samples-file", required=True, help="Path to the samples.txt file (one diploid sample name per line, in order).")
    # --- ADDED ARGUMENT for population file ---
    parser.add_argument("--pop-file", required=True, help="Path to the population CSV file (e.g., populations_v2.csv).")
    parser.add_argument(
        "-g", "--gen-time", type=float, default=25.0,
        help="Generation time in years (default: 25)."
    )
    args = parser.parse_args()

    # --- 1. Load Diploid Sample Names ---
    try:
        with open(args.samples_file) as f:
            sample_names = [line.strip() for line in f if line.strip()]
        num_diploid_samples = len(sample_names)
        num_expected_tips = num_diploid_samples * 2
        print(f"Loaded {num_diploid_samples} diploid sample names from {args.samples_file}.", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Samples file not found at {args.samples_file}", file=sys.stderr)
        sys.exit(1)

    # --- 1b. Load Population Data ---
    try:
        pop_df = pd.read_csv(args.pop_file)
        # Create a dictionary mapping sample name to population
        sample_to_pop = pop_df.set_index("sample")["population"].to_dict()
        print(f"Loaded population data for {len(sample_to_pop)} samples from {args.pop_file}.", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Population file not found at {args.pop_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading population file: {e}", file=sys.stderr)
        sys.exit(1)


    # --- 2. Load the Tree ---
    try:
        tns = dendropy.TaxonNamespace()
        tree = dendropy.Tree.get(
            path=args.nwk_in,
            schema="newick",
            taxon_namespace=tns,
            suppress_internal_node_taxa=False,
            suppress_leaf_node_taxa=False
        )
        print(f"Loaded tree with {len(tree.leaf_nodes())} tips from {args.nwk_in}.", file=sys.stderr)
        if len(tree.leaf_nodes()) != num_expected_tips:
             print(f"Warning: Tree has {len(tree.leaf_nodes())} tips, but expected {num_expected_tips}.", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input Newick file not found at {args.nwk_in}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading tree: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Relabel Tips (Haplotype & Population Aware, 1-Based Input) ---
    tip_map = {}

    # --- Build mapping FROM 1-based label ('1'-'48') TO SampleName_hapX_Population ---
    for i, sample_name in enumerate(sample_names):
        label_in_file_hap1 = str(2 * i + 1) # Expect '1', '3', '5', ...
        label_in_file_hap2 = str(2 * i + 2) # Expect '2', '4', '6', ...

        # Look up population, default to 'Unknown' if sample not in pop file
        population = sample_to_pop.get(sample_name, 'Unknown')
        if population == 'Unknown':
             print(f"Warning: Sample '{sample_name}' not found in population file.", file=sys.stderr)

        tip_map[label_in_file_hap1] = f"{sample_name}_hap1_{population}"
        tip_map[label_in_file_hap2] = f"{sample_name}_hap2_{population}"
    # --------------------------------------------------------------------------------

    processed_labels_count = 0
    missing_labels_in_tree = []
    processed_tip_objects = set()

    for tip in tree.leaf_node_iter():
        if tip.taxon and tip.taxon not in processed_tip_objects:
            current_label = tip.taxon.label
            cleaned_label = current_label.strip("'\"")

            if cleaned_label in tip_map:
                new_label = tip_map[cleaned_label]
                tip.taxon.label = new_label
                processed_labels_count += 1
                processed_tip_objects.add(tip.taxon)
            else:
                 try:
                     label_int = int(cleaned_label)
                     if not (1 <= label_int <= num_expected_tips):
                         print(f"Warning: Ignoring unexpected tip label '{cleaned_label}'.", file=sys.stderr)
                         processed_tip_objects.add(tip.taxon)
                     else:
                          missing_labels_in_tree.append(current_label)
                 except ValueError:
                     missing_labels_in_tree.append(current_label)

    if missing_labels_in_tree:
        unique_missing = sorted(list(set(missing_labels_in_tree)))
        print(f"Warning: Could not map expected labels: {', '.join(unique_missing)}", file=sys.stderr)

    if processed_labels_count != num_expected_tips:
         print(f"Warning: Relabeled {processed_labels_count} tips. Expected {num_expected_tips}.", file=sys.stderr)


    # --- 4. Rescale Branch Lengths ---
    scaling_factor = args.gen_time / 1_000_000.0
    for edge in tree.edges():
        if edge.length is not None:
            if edge.length < 0:
                 print(f"Warning: Negative edge length {edge.length} encountered. Setting to 0.", file=sys.stderr)
                 edge.length = 0.0
            edge.length *= scaling_factor

    # --- 5. Save the Modified Tree ---
    try:
        tree_string = tree.as_string("newick", suppress_rooting=True)
        with open(args.nwk_out, "w") as f:
            f.write(tree_string)
        print(f"Successfully relabeled tips (hap+pop), rescaled branches (factor={scaling_factor:.2e}), saved to {args.nwk_out}")

    except Exception as e:
        print(f"Error writing tree: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()