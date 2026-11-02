# Package pipeline

Run from the zip root (`theileria_interface_scripts_panels/`).

1. `01_scripts/01_interface_analysis/01_aligned_interface.R` — RWR, proximity, PCSF; writes `03_algorithm_tables/`
2. `01_scripts/01_interface_analysis/02_local_parasite_convergence.py` — local-core z = 4.66
3. `01_scripts/02_build_figure_data.R` — rebuilds `04_plot_data/`
4. `01_scripts/03_panels/0N_figure_N/0M_panel_*.R` — draws panels into `06_panel_output/`

Processed inputs are in `02_processed_inputs/`. Extra BPQ/module-name tables are in `05_extra_tables/`. Compose scripts are not in this deposit.
