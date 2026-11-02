# Theileria annulata–bovine interface scripts and panel data

Folder order is sequential at every level. Journal-page compose scripts are not included.

```
00_README.md
01_scripts/
  00_PIPELINE.md
  01_interface_analysis/     # 104/95, z = 4.66, 18 PCSF runs
  02_build_figure_data.R     # algorithm tables -> plot TSV
  03_panels/                 # 01_figure_1 ... 05_figure_5
02_processed_inputs/
  01_networks/
  02_evidence/
  03_id_map/
  04_inventory/
03_algorithm_tables/         # locked outputs (prefixed 01_, 02_, ...)
04_plot_data/                # 01_figure_1 ... 05_figure_5
05_extra_tables/             # module names, BPQ tables
06_panel_output/             # created when a panel script is run
```

```text
Rscript --vanilla 01_scripts/01_interface_analysis/01_aligned_interface.R
python 01_scripts/01_interface_analysis/02_local_parasite_convergence.py
Rscript --vanilla 01_scripts/02_build_figure_data.R
Rscript --vanilla 01_scripts/03_panels/01_figure_1/01_panel_a_included_studies.R
```

Locked checks after the interface scripts: `03_algorithm_tables/` files matching `pcsf_parameter_stability` (18/18), `pcsf_node_selection_frequency` / `pcsf_edge_selection_frequency` (≥0.50 → 104 / 95), `parasite_local_convergence` (z = 4.66, P = 0.0001).

Software versions, random seeds and STRING retrieval dates: `00_REPRODUCIBILITY.md`, `sessionInfo.txt`, `requirements.txt`.
