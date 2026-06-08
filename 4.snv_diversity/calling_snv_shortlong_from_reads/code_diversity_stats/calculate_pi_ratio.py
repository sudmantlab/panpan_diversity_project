import pandas as pd
import argparse
import sys

def calculate_ratios(within_file, pooled_file, output_file):
    """
    Merges pixy results and calculates pooled/within pi ratios.
    """
    try:
        # Load the within-population pi results
        within_df = pd.read_csv(within_file, sep='\t')
        print(f"Loaded {len(within_df)} rows from within-population results.")

        # Load the pooled pi results
        pooled_df = pd.read_csv(pooled_file, sep='\t')
        # We only need the window info and the pi value
        pooled_df = pooled_df[['chromosome', 'window_pos_1', 'window_pos_2', 'avg_pi']]
        pooled_df.rename(columns={'avg_pi': 'pi_pooled'}, inplace=True)
        print(f"Loaded {len(pooled_df)} windows from pooled results.")

        # Merge the two dataframes based on the window coordinates
        merged_df = pd.merge(
            within_df,
            pooled_df,
            on=['chromosome', 'window_pos_1', 'window_pos_2'],
            how='left'
        )
        print(f"Merged dataframes. Resulting shape: {merged_df.shape}")

        # Calculate the final ratio
        # The within-pi column from pixy is already named 'avg_pi'
        merged_df['pi_ratio'] = merged_df['pi_pooled'] / merged_df['avg_pi']
        
        # Save the final summary file
        merged_df.to_csv(output_file, index=False, na_rep='NA')
        print(f"Successfully calculated ratios and saved to {output_file}")

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge pixy results and calculate Pi_pooled/Pi_within.")
    parser.add_argument('--within', required=True, help='Path to the within-population pi file.')
    parser.add_argument('--pooled', required=True, help='Path to the pooled pi file.')
    parser.add_argument('--output', required=True, help='Path for the output summary CSV file.')
    
    args = parser.parse_args()
    calculate_ratios(args.within, args.pooled, args.output)
