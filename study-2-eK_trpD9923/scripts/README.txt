**Files
1. simulate_top_35_models.py
- This file contains the code to simulate the 35 models calibrated on the behavior of two strains - W3110 trpD9923 and W3110 trpD9923/pJLaroGfbr
- The output of this file is used to generate figure 6a
2. single_species.yaml
- This is the file with the medium data with which we start the bioreactor simulations

**Subfolders
1. eK_trpD9923 that contains all the scripts for generating data for the enhanced models of the reference strain, W3110 trpD9923
2. eK_trpD9923_d2 that contains the scripts for generating data for the enhanced models of the double mutant, W3110 trpD9923/pJLaroGfbrtktA
3. plot that contains all the scripts used to generate the plots for the following figures:
-- Main text figures --> Figure 6, Figure 7
-- Supplementary figures --> S.10, S.11, S.12, S.13, S.14

**Important points:
- The order in which the scripts are to be executed is given by the prefix attached to each file name in each folder.
- It is better to run all the scripts in eK_trpD9923 and eK_trpD9923_d2 before running the plots
- The plotting scripts themselves can be run in any order
- Some figures require the installation of seaborn/scikitlearn. This will be mentioned in the corresponding plot file
- Scripts starting with S. are for studies presented in the Supplementary notes.
- Further details about the various files are provided in the subfolders.
