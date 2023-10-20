This folder has the python scripts needed to generate the designs and simulations for the enhanced models, eK_trpD9923
The order in which the scripts are to be executed is given by the prefix attached to each file name.
Scripts starting with S. are for studies presented in the Supplementary notes.

** I_design_generation.py
- In this script we generate all the NOMAD designs for a set of given design constraints using the 13 eK_trpD9923 models
- The allowable fold change in concentration is set as a parameter inside the script
- Each run for an allowable fold change in concentrations takes typically around 1 hour

** II_design_robustness
- In this script we go through all the designs that were generated for each kinetic model (within 5% of the maximal objective
for that model)
- We identify the unique designs my membership
- We then apply each of the 34 unique designs to each of the kinetic models in an NRA setting (by setting the min. fold change
in activity for that enzyme to 1e-3)
- We then store the solution value that NRA gives us.
- In this manner we obtain a matrix of 34 x 13 NRA solutions
- The results of this will be used to get the top 5 designs (by average predicted inc. in anthranilate yield as per NRA)
- A run takes around 30 minutes.

** III_verification_of_top_designs_in_reactor
- We first select the top 5 designs based on the average predicted increase in anthranilate yield
- We will then simulate these designs in a bioreactor setting using the NRA suggested fold changes for each model
- The results of this script are used for plotting figure S.12 A and B
- A run takes around 40 minutes.

** IV_expression_sensitivity_of_top_designs
- Here we test the robustness of the designs to errors in experimental implementatoin
- we apply a 50% random perturbation 10 times to each combination of design and model
- Note that for the fold changes we use the mean of the NRA suggested enzyme fold changes for each enzyme in a design across the 13 kinetic models.
- We can apply to perturbation to either a single enzyme expression about its mean NRA predicted value while keeping the other 2 
at the mean NRA proposed values.
- To plot figure S.12C-F you need to run this script 4 times, for perturbing each enzyme individually and then all 3 together
- You can change which enzyme to perturb within the script.
- Each run consisting of 10x10x5 ODE simulations takes around 4 hours. 

** S.XI_old_designs_in_eK_trpD9923
- In this file, we implement the four designs generated using K_trpD9923 in the enhanced models, eK_trpD9923
- We use NRA to suggest the fold changes in enzyme activity for each of these designs in eK_trpD9923
- Note that we ues the same constraints that were used to generate the designs for eK_trpD9923 in I_design_generation.py
- We actually implement the top 5 designs from K_trpD9923 but we finally only plot 4 of them
- The results of this script are used to plot figure S.14
