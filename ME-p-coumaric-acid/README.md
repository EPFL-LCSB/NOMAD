# NOMAD application: improve p-coumaric acid production in S. cerevisiae
## About this subrepository
This subrepository contains the code required to reproduce the work presented
in the paper entitled "Kinetic-model-guided engineering of multiple S. cerevisiae strains 
improves p-coumaric acid production" (https://doi.org/10.1101/2024.12.15.628543).

It contains
1. The code used to develop engineering strains for the S. cerevisiae strain ST10284 built at DTU. (strain-design/scripts)
2. The 9 kinetic models representative of ST10284 (kinetic-modelling/parameters/kinetic_parameters_top_9_models.hdf5)
3. The code for plotting Figures 1, 2, 4, S2, and S3 (strain-design/scripts/plot)
4. The 5,000 steady states out of which we used the one closest to the mean for constructing the models (data/steady_states.csv)
5. The steady state profiles for the Vmax sensitivity analysis (data/steady_states_vmax_sensitivity.csv)
6. The output designs from each of the 9 kinetic models, along with the output of the robustness analysis (strain-design/output/designs)
7. The code used to generate kinetic models around the mean steady-state (kinetic-modelling). Note that this need not be run again and if you do run it, please name the output something else so you do not overwrite the existing parameters.
8. The codes used to compare our KMs with those in BRENDA (kinetic-modelling)

The output of the ODE simulations can be found in our Zenodo repository (https://doi.org/10.5281/zenodo.15432261). If you want to re-run the data postprocessing /plotting codes, please download them into the path strain-design/output/nonlinear-verification.
This will be necessary for running strain-design/scripts/IV_pca_yield_analysis.py and strain-design/scripts/plot/figure_2BC.py

The subrepository uses code from the mother repository - please follow the instructions in the NOMAD README to install NOMAD, pytfa and Skimpy.
We strongly recommend using the pre-built docker image instead of building your own.

For any queries, contact Bharath at bn312@cam.ac.uk



