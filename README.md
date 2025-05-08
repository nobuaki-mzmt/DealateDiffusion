# Data and Code for "Evaluating mate encounter and walking dispersal dynamics of termites using posture tracking and behavioral simulation"
 
## Article information

This repository provides access to the data and source code used for the manuscript  
**Evaluating mate encounter and walking dispersal dynamics of termites using posture tracking and behavioral simulation**

**Author:**  
**Nobuaki Mizumoto**  
Department of Entomology & Plant Pathology, Auburn University, Auburn, AL, 36849, USA<br>
**Paper DOI:** [TBA](XXX)

This study reanalyzed the videos
describes the tandem running behavior in the termite *Hodotermopsis sjostedti*, using high-resolution position data extracted from video tracking. The analysis focuses on identifying tandem runs, measuring their durations, detecting leader roles (male or female), and analyzing how behaviors vary with arena size.  
This repository includes raw tracking data and the Python and R scripts to analyze them.  
The models and labels are available at TBA.

## Repository Structure

- **`/analysis/codes/`**: Contains scripts for analysis.
  - **`analysis.R`**: R script for conducting statistical analysis and generating figures.
  - **`sleap_processing.py`**: Python script to process and clean data from SLEAP tracking outputs.
  
- **`/data_raw/`**: Contains raw `.h5` data produced by SLEAP.
- **`/data_fmt/`**: Contains formatted datasets generated during analysis.
- **`/output/`**: Contains output files, including figures and analysis results.

## Setup & Dependencies


This project is written in R. You’ll need the following packages:

```r
install.packages(c("stringr", "data.table", "arrow", "dplyr", "MASS", "ggplot2",
                   "patchwork", "knitr", "survival", "survminer", "zoo",
                   "cowplot", "coxme", "tidyr"))
```
This project also uses Python. You’ll need the following Python packages:

```ini
pandas==1.3.5
h5py==3.1.0
numpy==1.19.5
scipy==1.7.3
```

## Citation
TBA
@article{mizumoto2025, title={Observation of tandem running behavior in mating pairs of Asian dampwood termite, Hodotermopsis sjostedti}, author={Mizumoto, Nobuaki and Chambliss, William and Elijah, Carroll P and Nakazono, Tomohiro and Kanao, Taisuke}, journal={TBA}, year={2025}, doi={DOI} }

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
William Chambliss: wlc0018@auburn.edu  
Nobuaki Mizumoto: nzm0095@auburn.edu
