import tskit
import pandas as pd
# Load all 99 trees (indices 0–98)
ts_list = [tskit.load(f"hprc_ht2t_chr6_{i}.trees") for i in range(99)]
# Convert relative coordinates to T2T positions (offset = 26,390,105)
offset = 26390105
data = []
for ts in ts_list:
    for tree in ts.trees():
        data.append({
            "start": tree.interval.left + offset,
            "end": tree.interval.right + offset,
            "tmrca": tree.time(tree.root)
        })
df = pd.DataFrame(data)
df.to_csv("tmrca_data_t2t.csv", index=False)
print(f"Last position: {df['end'].max()}")
#Last position: 33135780.0
print(f"Genomic range: {df['start'].min()} – {df['end'].max()}")
# Expected output: 26390105.0 – 33135780.0
pop_df = pd.read_csv("populations.csv")  # columns: sample, population
pop_map = pop_df.set_index("sample")["population"].to_dict()


ts_list = [tskit.load(f"hprc_ht2t_chr6_{i}.trees") for i in range(99)]
windowed_diversity_list = []
for ts in ts_list:
    windows = np.arange(0, ts.sequence_length, 1000)
    windows = np.append(windows, ts.sequence_length)
    windowed_diversity = ts.diversity(mode='branch', windows=windows)
    windowed_diversity_list.append(windowed_diversity_list)
