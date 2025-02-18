=============================================================================================
**MaxTDA: Robust Statistical Inference for Maximal Persistence in Topological Data Analysis**                              
=============================================================================================

Description:
--------------------------------------------------------------------------------
This repository contains code and data used in the paper: 
**MaxTDA: Robust Statistical Inference for Maximal Persistence in Topological 
Data Analysis**

Directories:
--------------------------------------------------------------------------------
1. code: This contains Python and R scripts used to generate various results in 
         the manuscript
2. data: This contains pre-generated data used for some simulations
3. images: Figures from simulations and real data applications are saved here


Dependencies:
--------------------------------------------------------------------------------
The required python libraries are in the file: code/requirements.txt
- This can be installed with: %pip install -r requirements.txt
- The computations were done with Python version 3.11.10

The following R packages are required for various sections of the codes

- TDA: For various TDA computations (version 1.9)
- ripserr: For persistence computation (version 0.1.1)
- tidyverse: Miscellaneous data manipulation and making plots (version 2.0.0)
- foreach: For pararellel computations
- doParallel: For pararellel computations (version 1.0.17)
- scatterplot3d: For 3D plotting (version 0.3-44)
- nonlinearTseries: For estimating TDE parameters (version 0.2.12)


Usage:
--------------------------------------------------------------------------------
Instructions for running the codes are provided within each script. Most scripts
are in a pair of .R and .ipynb. The .ipynb is used to generate data whiles the
.R script is used for most TDA computations.

- annulus_topology_recovery.R: This looks at the topology recovery of the annulus
                               corrupted with noise 
- annulus_topology_recovery.ipynb: This generates the data together with the 
                                   debsity thresholding

- varying_data_distribution.R: This performs simulations where the data 
                               distribution is not uniform on all geometries
- varying_data_distribution.R: This contains code for the debsity thresholding

- graphical_abstract.R: This contains code that plots the density surface 
                        illustrating the smoothing effect of the method
- graphical_abstract.ipynb: This generates the data together with the 
                                   debsity thresholding

- vr_illustration.R: This contains code for an example illustration VR filtration
                                        
- application.R: This contains code for results in the Application section
- application.ipynb: This contains code for the debsity thresholding

- Utils.R: Functions for routine operations

- kde_functions.py: This contains python scripts that implements Algorithm 1

- pds_functions.py: Contains code for routine persistence diagrams computations 

- DTM_filtrations.py: This code performs DTM filtration. The script was taken 
                      from: https://github.com/GUDHI/TDA-tutorial

Contact:
--------------------------------------------------------------------------------
- Name: Sixtus Dakurah and Jessi Cisewski-Kehe
- Email: [sdakurah, jjkehe]@wisc.edu
================================================================================
