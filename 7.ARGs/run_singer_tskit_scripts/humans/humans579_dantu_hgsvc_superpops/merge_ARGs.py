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

    time_zero_nodes_added = False
    sample_num = 0

    for node_file_index, (node_file, branch_file, mutation_file) in enumerate(zip(node_files, branch_files, mutation_files)):
        block_start = block_coordinates[node_file_index]
        print(f"Processing block {node_file_index} (start={block_start})")

        # Load data
        node_time = np.loadtxt(node_file)
        edge_span = np.loadtxt(branch_file)
        edge_span = edge_span[edge_span[:, 2] >= 0, :]  # Filter invalid branches
        mutations = np.loadtxt(mutation_file)

        # Add nodes
        node_num = node_table.num_rows - sample_num
        min_time = 0
        for t in node_time:
            if t == 0:
                if node_file_index == 0:  # Add samples only from first block
                    node_table.add_row(flags=tskit.NODE_IS_SAMPLE)
                    sample_num += 1
            else:
                t = max(min_time + 1e-7, t)
                node_table.add_row(time=t)
                min_time = t

        # Add edges
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

        # Add mutations (skip those beyond chrom_length)
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
            if mut_node < sample_num:
                tables.mutations.add_row(site=site_id, node=mut_node, derived_state=str(int(mut[3])))
            else:
                tables.mutations.add_row(site=site_id, node=mut_node + node_num, derived_state=str(int(mut[3])))

    tables.sort()
    return tables.tree_sequence()

def main():
    parser = argparse.ArgumentParser(description="Merge ARG outputs")
    parser.add_argument("--file_table", required=True)
    parser.add_argument("--chrom_sizes", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load chromosome sizes
    chrom_sizes = {}
    with open(args.chrom_sizes) as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)

    # Load input files
    node_files, branch_files, mutation_files, block_coordinates = [], [], [], []
    with open(args.file_table) as f:
        for line in f:
            node, branch, mut, start = line.strip().split()
            node_files.append(node)
            branch_files.append(branch)
            mutation_files.append(mut)
            block_coordinates.append(int(start))

    # Get chromosome name and length
    #chrom = '_'.join(os.path.basename(node_files[0]).split('_')[:3])  #line for chimps
    chrom = os.path.basename(node_files[0]).split('_')[0]  # Just takes "chr1" (humans have different names) 
    chrom_length = chrom_sizes.get(chrom, 0)
    if chrom_length == 0:
        raise ValueError(f"Chromosome {chrom} not found in {args.chrom_sizes}")

    # Merge and save
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates, chrom_length)
    ts.dump(args.output)

if __name__ == "__main__":
    main()