# Replication Package: "Spatial Price Competition and Buyer Power in the U.S. Beef Packing Industry"

This repository contains the code and data needed to replicate the tables and figures in Moschini and Smith (forthcoming, *American Journal of Agricultural Economics*). The replication package includes raw output from the fixed-point algorithm used to compute the Bertrand-Nash equilibrium. For details on the model and algorithmic procedure, please refer to the manuscript.

## Directory Structure

### Replication Inputs
- `data/`: Input datasets, including:
  - County-level fed cattle sales volumes
  - Plant-level data on beef packers
- `helpers/`: Helper functions used by the MATLAB code
- `raw_output/`: Raw output from the equilibrium algorithm (Bertrand-Nash oligopsony), used to generate results
- `replicate_figures.py`: Python script to generate manuscript figures
- `replicate_results_tables.m`: MATLAB script to generate manuscript tables

### Replication Outputs
- `images/`: Figures generated by `replicate_figures.py`
- `results/`: Tables generated by `replicate_results_tables.m`

## Replication Instructions

1. **Run `replicate_results_tables.m`**
   - Output:
     - Tables 3–6 will be saved in `results/`
     - `county_results.csv` will be saved in `data/`, used as input for `replicate_figures.py`

2. **Set up your Python environment**
   - Required libraries:
     - `pandas`
     - `plotly`
     - `kaleido`
   - Install using pip if needed:
     ```bash
     pip install pandas plotly kaleido
     ```

3. **Run `replicate_figures.py`**
   - Output: Figures saved to `images/`:
     - **Figure 1**: `LogFedCattle.png`
     - **Figure 2**: `plants.png`
     - **Figure 4**:
       - Panel (a): `Price.png`
       - Panel (b): `Markdown.png`
       - Panel (c): `MarkdownSpatial.png`
       - Panel (d): `MarkdownMulti.png`
       - Panel (e): `MarkdownContract.png`
       - Panel (f): `CapEffect.png`

## Citation

Moschini, G. and T.J. Smith. "*Spatial Competition and Buyer Power in the U.S. Beef Packing Industry*." *American Journal of Agricultural Economics*, forthcoming.
