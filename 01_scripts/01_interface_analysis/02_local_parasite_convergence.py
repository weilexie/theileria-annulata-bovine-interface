"""Primary local parasite convergence test."""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ALG = os.path.join(ROOT, "03_algorithm_tables")
N_PERM = 10000
SEED = 20260902


def main():
    frame = pd.read_csv(os.path.join(ALG, "36_rwr_cross_layer_gain.tsv"), sep="\t")
    frame = frame[
        frame["layer"].eq("pathogen")
        & frame["accessible_annot"].astype(str).str.upper().eq("TRUE")
    ].copy()
    frame["strict_member"] = frame["strict_member"].astype(str).str.upper().eq("TRUE")
    frame = frame.dropna(subset=["degree", "same_side_diffusion_percentile"])
    selected = frame["strict_member"]
    if selected.sum() == 0 or selected.sum() >= len(frame):
        raise ValueError("local accessible layer must contain selected and background nodes")
    frame["degree_bin"] = pd.qcut(frame["degree"], q=4, duplicates="drop")
    counts = frame.loc[selected, "degree_bin"].value_counts().sort_index()
    pools = {key: group.index.to_numpy() for key, group in frame.groupby("degree_bin", observed=True)}
    if any(len(pools[key]) < count for key, count in counts.items()):
        raise ValueError("degree-bin pool is smaller than selected local layer")
    observed = frame.loc[selected, "same_side_diffusion_percentile"].mean()
    rng = np.random.default_rng(SEED)
    null = np.empty(N_PERM, dtype=float)
    for i in range(N_PERM):
        indices = []
        for key, count in counts.items():
            indices.extend(rng.choice(pools[key], size=int(count), replace=False))
        null[i] = frame.loc[indices, "same_side_diffusion_percentile"].mean()
    null_mean = null.mean()
    null_sd = null.std(ddof=1)
    summary = pd.DataFrame([{
        "analysis": "accessible_recurrent_core_local_convergence",
        "selected_n": int(selected.sum()),
        "accessible_network_n": int(len(frame)),
        "selected_fraction": float(selected.mean()),
        "observed_mean_same_side_diffusion": float(observed),
        "null_mean": float(null_mean),
        "null_sd": float(null_sd),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "z": float((observed - null_mean) / null_sd),
        "empirical_p_one_sided": float((1 + np.sum(null >= observed - 1e-12)) / (N_PERM + 1)),
        "permutations": N_PERM,
        "null_scheme": "accessible-layer degree-bin matched without replacement",
        "primary_layer": "accessible_annot=TRUE, pathogen network nodes",
        "selected_layer": "strict_member=TRUE within primary layer",
    }])
    summary.to_csv(os.path.join(ALG, "26_parasite_local_convergence.tsv"), sep="\t", index=False)
    pd.DataFrame({"iteration": np.arange(1, N_PERM + 1), "null_statistic": null}).to_csv(
        os.path.join(ALG, "27_parasite_local_convergence_null.tsv"), sep="\t", index=False
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
