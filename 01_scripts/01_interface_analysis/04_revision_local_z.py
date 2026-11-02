# Same qcut local test as 03_scripts/03_interface/02_local_parasite_convergence.py
# Reads the no-1.25 / term>=4 rerun RWR scores.
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
rerun = ROOT / "03_algorithm_tables"
src = rerun / "rerun_rwr_nodes.tsv"
out = rerun / "rerun_local_convergence.tsv"
frame = pd.read_csv(src, sep="\t")
frame = frame[frame["layer"].eq("pathogen") & frame["accessible_annot"].astype(str).str.upper().eq("TRUE")].copy()
frame["strict_member"] = frame["strict_member"].astype(str).str.upper().eq("TRUE")
frame = frame.dropna(subset=["degree", "same_side_diffusion_percentile"])
selected = frame["strict_member"]
frame["degree_bin"] = pd.qcut(frame["degree"], q=4, duplicates="drop")
counts = frame.loc[selected, "degree_bin"].value_counts().sort_index()
pools = {key: group.index.to_numpy() for key, group in frame.groupby("degree_bin", observed=True)}
observed = frame.loc[selected, "same_side_diffusion_percentile"].mean()
rng = np.random.default_rng(20260902)
null = np.empty(10000, dtype=float)
for i in range(10000):
    indices = []
    for key, count in counts.items():
        indices.extend(rng.choice(pools[key], size=int(count), replace=False))
    null[i] = frame.loc[indices, "same_side_diffusion_percentile"].mean()
null_mean = null.mean()
null_sd = null.std(ddof=1)
summary = pd.DataFrame([{
    "selected_n": int(selected.sum()),
    "accessible_network_n": int(len(frame)),
    "observed_mean": float(observed),
    "null_mean": float(null_mean),
    "null_sd": float(null_sd),
    "z": float((observed - null_mean) / null_sd),
    "empirical_p_one_sided": float((1 + np.sum(null >= observed - 1e-12)) / 10001),
}])
summary.to_csv(out, sep="\t", index=False)
print(summary.to_string(index=False))
