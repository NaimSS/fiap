# Exact dynamic programming for the Fault Indicator Allocation Problem

Companion code and data for the paper.

## `src/`

The solvers. 

`linear_dp/` holds the six exact DPs for path networks (main trunk): the naive $O(n^2m)$ DP (`dp_n2k`), CHT (`dp_nk`), SMAWK (`dp_nk_smawk`), LARSCH (`dp_nk_larsch`) and their Lagrangian versions (`dp_nlogW`, `dp_nlogW_larsch`). 

`general_dp/` holds the tree DP: `lagrange_solver` (Lagrangian search, the main solver), `default_solver` (all budgets, no Lagrangian search), and `lagrange_solver_without_hld`.

## `data/`

`instances/` has the 26 power distribution networks used in the paper: `<name>.fi` is the network and `<name>.fi.diameter` adds the main trunk. Edge lengths are stored as `floor(km * 10^5)`. `results/` has the CSVs behind the paper's tables and figures.

## `scripts/`

`bench/` holds the benchmark runners that produce the CSVs in `data/results/`; each script documents its options and output at the top. `viz/` draws a network with its optimal indicator placement.

## AI Use Disclosure

AI-assisted tools were used to write parts of the scripts, visualization tools, documentation and algorithm implementations. All algorithm code was reviewed and tested by the authors.
