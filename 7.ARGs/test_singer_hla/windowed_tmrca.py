import tskit
import numpy as np
import pandas as pd
offset = 26390105  # T2T start coordinate
ts_list = [tskit.load(f"hprc_ht2t_chr6_{i}.trees") for i in range(99)]
data = []

for ts in ts_list:
    windows = np.arange(0, ts.sequence_length, 1000)
    windows = np.append(windows, ts.sequence_length)  
    tree_tmrcas = np.array([tree.time(tree.root) for tree in ts.trees()])
    tree_lefts = np.array([tree.interval.left for tree in ts.trees()])
    tree_rights = np.array([tree.interval.right for tree in ts.trees()])
    
    for i in range(len(windows) - 1):
        window_start = windows[i]
        window_end = windows[i + 1]
        
        overlaps = (tree_lefts < window_end) & (tree_rights > window_start)
        if np.any(overlaps):
            overlap_lengths = np.minimum(tree_rights[overlaps], window_end) - np.maximum(tree_lefts[overlaps], window_start)
            avg_tmrca = np.average(tree_tmrcas[overlaps], weights=overlap_lengths)
        else:
            avg_tmrca = np.nan
        
        data.append({
            "start": window_start + offset,
            "end": window_end + offset,
            "tmrca": avg_tmrca
        })

df = pd.DataFrame(data)
df.to_csv("windowed_tmrca_t2t.csv", index=False)
df_avg = df.groupby(["start", "end"], as_index=False)["tmrca"].mean()
df_avg.to_csv("windowed_tmrca_t2t_avg.csv", index=False)