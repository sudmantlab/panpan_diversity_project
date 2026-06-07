import argparse
import numpy as np
import tskit
import os

def read_long_ARG(node_files, branch_files, mutation_files, block_coordinates, chrom_length):
    if len(node_files) != len(branch_files) or len(node_files) != len(block_coordinates):
        raise ValueError("Mismatched input file lists")

    tables = tskit.TableCollection(sequence_length=chrom_length)
    node_table = tables.nodes
    branch_table = tables.edges

    sample_num = 0

    for node_file_index, (node_file, branch_file, mutation_file) in enumerate(zip(node_files, branch_files, mutation_files)):
        block_start = block_coordinates[node_file_index]
        print(f"Processing block {node_file_index} (start={block_start})")

        # Load node data. np.loadtxt will raise FileNotFoundError if missing.
        node_time = np.loadtxt(node_file)

        # Robustly load branch file, skipping malformed lines.
        # The `open()` command will raise FileNotFoundError if the file is missing.
        valid_edges = []
        with open(branch_file) as f:
            for i, line in enumerate(f):
                parts = line.strip().split()
                if len(parts) == 4:
                    valid_edges.append([float(p) for p in parts])
                else:
                    print(f"WARNING: Skipping malformed line {i + 1} in {branch_file}: {line.strip()}")
        edge_span = np.array(valid_edges) if valid_edges else np.empty((0, 4))
        
        # Robustly load mutation file, skipping malformed lines.
        # `open()` will raise an error if the file is missing.
        valid_mutations = []
        if os.path.getsize(mutation_file) > 0:
            with open(mutation_file) as f:
                for i, line in enumerate(f):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        valid_mutations.append([float(p) for p in parts[:4]])
                    else:
                        print(f"WARNING: Skipping malformed line {i + 1} in {mutation_file}: {line.strip()}")
        mutations = np.array(valid_mutations) if valid_mutations else np.empty((0, 4))
        
        edge_span = edge_span[edge_span[:, 2] >= 0, :]

        # Add nodes
        node_num = node_table.num_rows - sample_num
        min_time = 0
        for t in node_time:
            if t == 0:
                if node_file_index == 0:
                    node_table.add_row(flags=tskit.NODE_IS_SAMPLE)
                    sample_num += 1
            else:
                t = max(min_time + 1e-7, t)
                node_table.add_row(time=t)
                min_time = t

        # Add edges
        if edge_span.shape[0] > 0:
            parent_indices = np.array(edge_span[:, 2], dtype=np.int32)
            child_indices = np.array(edge_span[:, 3], dtype=np.int32)
            parent_indices[parent_indices >= sample_num] += node_num
            child_indices[child_indices >= sample_num] += node_num

            branch_table.append_columns(
                left=edge_span[:, 0] + block_start,
                right=edge_span[:, 1] + block_start,
                parent=parent_indices,
                child=child_indices
            )

        # Add mutations
        if mutations.shape[0] > 0:
            mut_pos = 0
            for mut in mutations:
                pos = mut[0] + block_start
                if pos >= chrom_length:
                    continue
                if pos != mut_pos:
                    tables.sites.add_row(position=pos, ancestral_state='0')
                    mut_pos = pos
                site_id = tables.sites.num_rows - 1
                mut_node = int(mut[1])
                derived_state = str(int(mut[3]))
                if mut_node < sample_num:
                    tables.mutations.add_row(site=site_id, node=mut_node, derived_state=derived_state)
                else:
                    tables.mutations.add_row(site=site_id, node=mut_node + node_num, derived_state=derived_state)

    tables.sort()
    return tables.tree_sequence()

def main():
    parser = argparse.ArgumentParser(description="Merge ARG outputs")
    parser.add_argument("--file_table", required=True)
    parser.add_argument("--chrom_sizes", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    chrom_sizes = {}
    with open(args.chrom_sizes) as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)

    if not os.path.exists(args.file_table) or os.path.getsize(args.file_table) == 0:
        print(f"WARNING: Input file table is missing or empty: {args.file_table}. Creating empty output.")
        # Create an empty tree sequence if the file list is empty
        try:
            # Heuristic to get chrom name from output file path for setting sequence_length
            output_basename = os.path.basename(args.output)
            chrom_name = '_'.join(output_basename.split('_')[:-1])
            seq_len = chrom_sizes[chrom_name]
        except KeyError:
             raise ValueError(f"Cannot determine chromosome from output file name {args.output} to create empty tree.")
        tskit.TableCollection(sequence_length=seq_len).tree_sequence().dump(args.output)
        return

    node_files, branch_files, mutation_files, block_coordinates = [], [], [], []
    with open(args.file_table) as f:
        for line in f:
            node, branch, mut, start = line.strip().split()
            node_files.append(node)
            branch_files.append(branch)
            mutation_files.append(mut)
            block_coordinates.append(int(start))

    chrom_from_file = '_'.join(os.path.basename(node_files[0]).split('_')[:3])
    chrom_length = chrom_sizes.get(chrom_from_file)
    if chrom_length is None:
        raise ValueError(f"Chromosome '{chrom_from_file}' not found in sizes file.")

    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates, chrom_length)
    ts.dump(args.output)

if __name__ == "__main__":
    main()