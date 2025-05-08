# Data and Code for "Evaluating mate encounter and walking dispersal dynamics of termites using posture tracking and behavioral simulation"
 
## Article information

This repository provides access to the data and source code used for the manuscript  
**Evaluating mate encounter and walking dispersal dynamics of termites using posture tracking and behavioral simulation**

**Author:**  
**Nobuaki Mizumoto**  
Department of Entomology & Plant Pathology, Auburn University, Auburn, AL, 36849, USA<br>
**Paper DOI:** [TBA](XXX)

This study reanalyzed the videos of termite mate searching behavior published in [Mizumoto and Dobata 2019](https://doi.org/10.1126/sciadv.aau6108), using the deep-learning posture tracking software, [SLEAP](https://sleap.ai). The analyis unwrapped the movement trajectories in a petri-dish arena to the open space, to compare the dispersal ability of termite mate searchers and nest searchers. Then, using the movement paramters, I developed a movement simulation to investigate the effect of density and light attraction on termite mating encounters.  
This repository includes raw tracking data and the Python, R, and C++ scripts.  
The models and labels are available at TBA.

## Repository Structure

- **`/analysis/codes/`**: Contains scripts for analysis.
  - **`data_prep.py`**: Python script to process and clean data from SLEAP tracking outputs.
  - **`trajectory_align.py`**: Python script to combine trajectories output to create Figure S1-3.
  - **`fmt_trajectory.R`**: R script for preprocessing all data to generate RDA files for statistical analysis and figure production.
  - **`output.R`**: R script for conducting statistical analysis and generating figures.
  - **`sim.R`**: R script for running simulations and output simulation results in figures.
  - **`sim.cpp`**: C++ script for simulations in a form of function.
- **`/data_raw/`**: Contains raw `.h5` data produced by SLEAP.
- **`/data_fmt/`**: Contains formatted datasets generated during analysis.
- **`/output/`**: Contains output files, including figures and analysis results.

## Setup & Dependencies
This project is written in R, Python, and C++, tested on Windows 11 (64-bit). Following is the environments.

### R Session Info

R version R version 4.4.1 (2024-06-14 ucrt)

#### Attached Packages
Rcpp, data.table, arrow, viridis, viridisLite, patchwork, ggplot2, circular, CircStats, boot, nlme, fitdistrplus, survival, MASS, car, carData, lme4, Matrix, dplyr

#### Reproducing R Environment
```r
remotes::install_version("Rcpp", version = "1.0.14")
remotes::install_version("data.table", version = "1.17.0")
remotes::install_version("arrow", version = "19.0.1")
remotes::install_version("viridis", version = "0.6.5")
remotes::install_version("viridisLite", version = "0.4.2")
remotes::install_version("patchwork", version = "1.3.0")
remotes::install_version("ggplot2", version = "3.5.1")
remotes::install_version("circular", version = "0.5-1")
remotes::install_version("CircStats", version = "0.2-6")
remotes::install_version("boot", version = "1.3-30")
remotes::install_version("nlme", version = "3.1-164")
remotes::install_version("fitdistrplus", version = "1.2-2")
remotes::install_version("survival", version = "3.6-4")
remotes::install_version("MASS", version = "7.3-60.2")
remotes::install_version("car", version = "3.1-3")
remotes::install_version("carData", version = "3.0-5")
remotes::install_version("lme4", version = "1.1-36")
remotes::install_version("Matrix", version = "1.7-0")
remotes::install_version("dplyr", version = "1.1.4")
```

### Python Environment
Python 3.11.4 | sys, os, glob, math, h5py 3.13.0, numpy 1.25.0, pandas 2.2.3, scipy 1.15.2, feather 0.4.1, PIL 11.2.0

## Citation
TBA
@article{mizumoto2025, title={Evaluating mate encounter and walking dispersal dynamics of termites using posture tracking and behavioral simulation}, author={Mizumoto, Nobuaki}, journal={TBA}, year={2025}, doi={DOI} }

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Nobuaki Mizumoto: nzm0095@auburn.edu
