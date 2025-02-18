# MaxTDA: Robust Statistical Inference for Maximal Persistence in Topological Data Analysis

## Description
This repository contains code and data used in the paper:
**MaxTDA: Robust Statistical Inference for Maximal Persistence in Topological Data Analysis**.

## Directory Structure
- **`code/`**: Contains Python and R scripts used to generate various results in the manuscript.
- **`data/`**: Contains pre-generated data used for some simulations.
- **`images/`**: Stores figures from simulations and real data applications.

## Dependencies
### Python
The required Python libraries are listed in `code/requirements.txt`. Install them using:
```bash
pip install -r requirements.txt
```
The computations were performed using **Python 3.11.10**.

### R
The following R packages are required:
- **TDA** (version 1.9) - For various TDA computations.
- **ripserr** (version 0.1.1) - For persistence computation.
- **tidyverse** (version 2.0.0) - For data manipulation and plotting.
- **foreach** - For parallel computations.
- **doParallel** (version 1.0.17) - For parallel computations.
- **scatterplot3d** (version 0.3-44) - For 3D plotting.
- **nonlinearTseries** (version 0.2.12) - For estimating TDE parameters.

Install missing R dependencies using:
```r
install.packages(c("TDA", "ripserr", "tidyverse", "foreach", "doParallel", "scatterplot3d", "nonlinearTseries"))
```

## Usage
Each script contains inline instructions for execution. Most scripts have both `.R` and `.ipynb` versions:
- **`.ipynb` files**: Generate data and perform density thresholding.
- **`.R` scripts**: Perform TDA computations.

### Key Scripts
#### **Topology Recovery of Annulus**
- `annulus_topology_recovery.R`: Analyzes topology recovery of an annulus with noise.
- `annulus_topology_recovery.ipynb`: Generates data and performs density thresholding.

#### **Varying Data Distribution**
- `varying_data_distribution.R`: Simulates scenarios where data distribution is non-uniform.
- `varying_data_distribution.ipynb`: Performs density thresholding.

#### **Graphical Abstract**
- `graphical_abstract.R`: Plots density surfaces illustrating the smoothing effect.
- `graphical_abstract.ipynb`: Generates data and performs density thresholding.

#### **Other Scripts**
- `vr_illustration.R`: Provides an illustration of Vietoris-Rips (VR) filtration.
- `application.R`: Contains code for the Application section.
- `application.ipynb`: Generates data and performs density thresholding.
- `Utils.R`: Contains utility functions for routine operations.
- `kde_functions.py`: Implements **Algorithm 1**.
- `pds_functions.py`: Computes routine persistence diagrams.
- `DTM_filtrations.py`: Performs **DTM filtration** (adapted from [GUDHI's TDA-tutorial](https://github.com/GUDHI/TDA-tutorial)).

## Contact
- **Sixtus Dakurah & Jessi Cisewski-Kehe**
- **Email**: [sdakurah, jjkehe]@wisc.edu