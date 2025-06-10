# BalancedMinimumEvolution

This repository contains Julia code for the Balanced Minimum Evolution (BME) problem, including beam search algorithms and experimental scripts.

## Getting Started

### 1. Install Julia
Download and install Julia (version 1.11 or later recommended) from [https://julialang.org/downloads/](https://julialang.org/downloads/).

### 2. Clone the Repository
Clone this repository to your local machine:
```sh
git clone <repository-url>
cd BalancedMinimumEvolution
```

### 3. Instantiate the Julia Environment
Open Julia in the project directory and activate the environment:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
This will install all required dependencies as specified in `Project.toml` and `Manifest.toml`.

### 4. Load the Package
To use the package in Julia, run:
```julia
using BalancedMinimumEvolution
```

## Script Explanations

The `scripts/` directory contains the following main scripts:

- **beam_experiments.jl**: Runs the main set of experiments for the Balanced Minimum Evolution (BME) problem. It benchmarks different algorithms (beam search, beam search with LP, and FastME variants) on a variety of datasets, records solution quality and runtime, and saves results to `data/experimental_results.csv` and lower bounds to `data/lower_bounds.csv`. This is the main script to reproduce the experimental results in the paper.

- **generate_RDSM.jl**: Generates synthetic datasets of symmetric doubly stochastic matrices (RDSM instances) of varying sizes. The matrices are saved as text files in the `data/RDSM*/` directories. Use this script to create the random RDSM datasets used in the experiments.

- **generate_RIM.jl**: Generates synthetic datasets of symmetric integer matrices (RIM instances) of varying sizes. The matrices are saved as text files in the `data/RIM*/` directories. Use this script to create the random RIM datasets used in the experiments.

- **latex_tables.jl**: Processes the experimental results and lower bounds, then generates LaTeX tables summarizing the performance of the algorithms (e.g., solution quality, gaps, and runtimes). The tables are output as `.tex` files (such as `beam_tables.tex`, `beam_tables_gap_A.tex`, and `beam_tables_time.tex`) in the project root, ready for inclusion in publications.

---

- Use `generate_RDSM.jl` and `generate_RIM.jl` to generate synthetic datasets.
- Use `beam_experiments.jl` to run all experiments and collect results.
- Use `latex_tables.jl` to generate publication-ready tables from the results.

## Results

- Results are stored in the `data/` directory.
- Output files include CSVs with experimental results (e.g., `experimental_results.csv`, `lower_bounds.csv`) and subdirectories for each dataset.
- LaTeX tables summarizing results are generated as `beam_tables.tex`, `beam_tables_gap_A.tex`, `beam_tables_time.tex`, etc., in the project root.

## Notes
- Ensure you have write permissions for the `data/` directory to allow scripts to save results.
- For reproducibility, use the provided `Manifest.toml` to guarantee the same package versions.

## Support
For questions or issues, please open an issue in this repository.
