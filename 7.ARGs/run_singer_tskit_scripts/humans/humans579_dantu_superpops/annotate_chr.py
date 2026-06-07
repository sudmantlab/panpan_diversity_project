import os
import subprocess
import pandas as pd

# Configuration
chrom = "chr21"
metrics_file = f"results/{chrom}_metricsbyPop.csv"  # Changed to metricsbyPop
gene_annotation = "ht2t_CHM13_renamed_chr_gene.bed.gz"
output_file = f"results/{chrom}_metricsbyPop_annotated.csv"  # Changed output name

# Create unique temporary file names
temp_files = [
    f"temp_{chrom}_metrics.bed",
    f"temp_{chrom}_metrics_sorted.bed",
    f"temp_{chrom}_genes_sorted.bed",
    f"temp_{chrom}_annotated.bed"
]

try:
    # Step 1: Convert metrics CSV to BED format
    print("Converting metrics to BED format...")
    print(f"Input file: {metrics_file}")
    
    # Read the input file with correct column names
    df = pd.read_csv(metrics_file)
    
    # Write to BED format (just chromosome, start, end)
    df[["chromosome", "start", "end"]].to_csv(
        temp_files[0], 
        sep='\t', 
        index=False, 
        header=False
    )

    # Step 2: Sort metrics BED file
    print("Sorting metrics BED file...")
    with open(temp_files[1], "w") as sorted_out:
        subprocess.run(
            f"bedtools sort -i {temp_files[0]}",
            shell=True, check=True, stdout=sorted_out
        )

    # Step 3: Sort gene annotation BED file
    print("Sorting gene annotations...")
    with open(temp_files[2], "w") as genes_out:
        subprocess.run(
            f"zcat {gene_annotation} | bedtools sort",
            shell=True, check=True, stdout=genes_out
        )

    # Step 4: Annotate metrics with genes using bedtools map
    print("Annotating metrics with genes...")
    with open(temp_files[3], "w") as annotated_out:
        subprocess.run(
            f"bedtools map -a {temp_files[1]} -b {temp_files[2]} -c 4 -o collapse",
            shell=True, check=True, stdout=annotated_out
        )

    # Step 5: Merge gene annotations back with original data
    print("Creating final CSV output...")
    
    # Read the gene annotations
    gene_annotations = pd.read_csv(
        temp_files[3], 
        sep='\t', 
        header=None,
        names=["chromosome", "start", "end", "genes"]
    )
    
    # Merge with original dataframe
    annotated_df = pd.merge(
        df,
        gene_annotations,
        on=["chromosome", "start", "end"],
        how="left"
    )
    
    # Fill NaN genes with empty string
    annotated_df["genes"] = annotated_df["genes"].fillna("")
    
    # Reorder columns to match input with genes added at the end
    column_order = [
        "chromosome", "start", "end", "population", 
        "avg_tmrca", "T_pooled", "T_within", 
        "Tpooled_Twithin_ratio", "genes"
    ]
    annotated_df = annotated_df[column_order]
    
    # Save as CSV
    annotated_df.to_csv(output_file, index=False)
    print(f"Successfully created annotated file: {output_file}")

except Exception as e:
    print(f"Error during processing: {str(e)}")
    import traceback
    traceback.print_exc()

finally:
    # Cleanup temporary files
    print("Cleaning up temporary files...")
    for f in temp_files:
        if os.path.exists(f):
            os.remove(f)
    print("Done!")
