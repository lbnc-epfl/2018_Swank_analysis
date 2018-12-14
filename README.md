# 2018_Swank_analysis
This repository contains data and code to accompany the article 'Cell-free gene regulatory network engineering with synthetic transcription factors' by Zoe Swank, Nadanai Laohakunakorn, and Sebastian Josef Maerkl (https://doi.org/10.1101/407999). More information can be found in the paper or on our [lab website](http://lbnc.epfl.ch). Notebooks may be viewed directly in GitHub or using the Jupyter [nbviewer](https://nbviewer.jupyter.org/github/lbnc-epfl/2018_Swank_analysis/tree/master/), which tends to work more reliably. 

## Requirements

All code was tested on `python 3.6.4`. We use `numpy 1.13.3`, `scipy 1.0.0`, and `jupyter 1.0.0`. MCMC fitting was carried out using `emcee 2.2.1`. Please see the requirements file for more information. 

## Notebooks

`NB_bar_charts.ipynb` generates bar charts for Figures 4 and 5A-B.

`NB_models_thermo.ipynb` defines and fits thermodynamic models to cooperative dose response curves, which are used to generate plots for Figure 5 and Supplementary Figure S4.

`NB_models_Hill.ipynb` defines and fits Hill function models to cooperative dose response curves, and generates Supplementary Figure S6.

`NB_models_helix.ipynb` defines and fits the helical promoter model, and generates plots for Figure 5.

`NB_models_sens.ipynb` determines the dose response sensitivity for cooperative repressors and generates Supplementary Figure S5.

`NB_models_sitetuning.ipynb` plots the results of MCMC sampling for the site tuning data, used to generate Supplementary Figure S4.

`NB_TXTLdynamics.ipynb` analyzes the dynamics of TXTL on-chip and on the plate reader, and generates Supplementary Figure S3.

`helper.py` contains a few useful functions which are required by the notebooks.

## File structure

All notebooks and `helper.py` are contained in the top-level directory; additionally the notebooks require `/data/`, `/plots/`, and `/output/` subdirectories from which to read data, and into which plots and MCMC sampling output files are written. Files written into the `/output/` the output directory are ignored in this repository due to their large size but may be regenerated locally by running the notebooks.

## Data

The `/data/` directory contains the following files:

Cooperative dose response data:

	2site_coop_PDZ.csv
	2site_coop_GCN.csv
	2site_coop_AAGCN.csv

Helical library data:

	distance_chip_8A.csv

TXTL dynamics data from chip and plate reader: 

	dynamics_chip.csv
	dynamics_PR.csv

Repressor characterization data on-chip:

	repressor_1site_FR.csv
	repressor_1site_KD.csv
	repressor_1siteslide_KD.csv
	repressor_1siteslidedwn_FR.csv
	repressor_1siteslideup_FR.csv
	repressor_2site_FR.csv
	repressor_2site_KD.csv

Cooperative repressor characterization plate reader data:

	repressor_platereader.csv

## Output

All generated plots go to the `/plots/` subfolder. MCMC samples go into the `/output/` folder. These files may be very large.