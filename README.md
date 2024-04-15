[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10976282.svg)](https://doi.org/10.5281/zenodo.10976282)
# ifnl4-structure-paper

This is the repo for the MD data in the paper ....

The trajectories used for the analysis can be found here: https://zenodo.org/doi/10.5281/zenodo.10975902

The analysis.ipynb contains convergence analysis, RMSF, and ion analysis.
The angle_analysis_IL10RB.ipynb contains the analysis for calculating the angle between the three domains.

## Packages
Requires:
- MDAnalysis
- Numpy
- Matplotlib

## Install
Create conda environment from YAML:
```sh
conda create -f environment.yml
```
Activate the conda environment: 
```sh
conda activate ifnl4-analysis
```

### How to use
1) Clone this repo
   ```sh
   git clone https://github.com/bylehn/ifnl4-structure-paper
   ```
2) Follow the installation steps above for the conda environment
3) Download the trajectory files at https://zenodo.org/doi/10.5281/zenodo.10975902

4) Open notebooks and run them
