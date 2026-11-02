# Reproducibility notes

## Software

See `sessionInfo.txt` (R 4.5.1 and attached packages) and `requirements.txt` (Python).

Main R packages: igraph 2.2.2, Matrix 1.7.3, matrixStats 1.5.0, PCSF 0.99.1, dplyr 1.2.0, readr 2.2.0, stringr 1.6.0, tidyr 1.3.2, magick 2.9.1, ggplot2 4.0.2, ggrepel 0.9.8, patchwork 1.3.2, cowplot 1.2.0, circlize 0.4.17, ragg 1.5.2.

Python 3.14.3: pandas 2.3.3, numpy 2.4.3.

## STRING graphs

STRING v12 (Szklarczyk et al., 2023), REST API `https://string-db.org/api`.

- Bos taurus (taxon 9913), combined score ≥700: retrieved 1 September 2026 (`02_processed_inputs/01_networks/01_host_string_network_700.tsv`).
- T. annulata (taxon 5874), combined score ≥400: retrieved 2 September 2026 (`02_processed_inputs/01_networks/02_parasite_string_network_400.tsv`).

## Random seeds

| Step | Seed | Permutations |
|---|---|---|
| Two-sided convergence (host strict, host relaxed, pathogen strict, pathogen relaxed) | 20260903, 20260904, 20260905, 20260906 | 10,000 |
| Module proximity | 20260902 | 1,000 |
| Local-core test (`02_local_parasite_convergence.py`) | 20260902 | 10,000 |
| Label-permutation null | 20260903 | 10,000 |
| Degree-preserving rewire null | 20260904 | 5,000 |
| Leiden (host / parasite) | 20260905 / 20260906 | — |
| Figure-data network layouts (host / parasite) | 20260903 / 20260904 | — |
| Figure-data PCSF / atlas layouts | 20260905, 20260906 | — |
