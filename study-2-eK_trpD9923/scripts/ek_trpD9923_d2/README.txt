This folder has the python scripts needed to generate the designs and simulations for the enhanced models of the double mutant, eK_trpD9923_d2
The order in which the scripts are to be executed is given by the prefix attached to each file name.
Scripts starting with S. are for studies presented in the Supplementary notes.

** I_design_generation.py
- In this script we generate all the NOMAD designs for a set of given design constraints using the 13 eK_trpD9923_d2 models
- The allowable fold change in concentration is set as a parameter inside the script
- We first simulate the double mutant
- We then obtain the steady-state concentrations and fluxes at t=9 hours
- We calculate the control coefficients at this point and use them for NRA
- Each run for an allowable fold change in concentrations takes typically around 1 hour

** II_design_robustness
- In this script we go through all the designs that were generated for each kinetic model (within 5% of the maximal objective
for that model)
- We identify the unique designs my membership
- We then apply each of the 13 unique designs to each of the kinetic models in an NRA setting (by setting the min. fold change
in activity for that enzyme to 1e-3)
- We then store the solution value that NRA gives us.
- In this manner we obtain a matrix of 13 x 13 NRA solutions
- The results of this will be used to get the top 5 designs (by average predicted inc. in anthranilate yield as per NRA)
- A run takes around 20 minutes.

** III_verification_of_top_designs_in_reactor
- We first select the top 5 designs based on the average predicted increase in anthranilate yield
- We will then simulate these designs in a bioreactor setting using the NRA suggested fold changes for each model
- The results of this script are used for plotting figure S.13A
- A run takes around 40 minutes.

** IV_expression_sensitivity_of_top_designs
- Here we test the robustness of the designs to errors in experimental implementatoin of all three enzyme simultaneously
- we apply a 50% random perturbation 10 times to each combination of design and model
- Note that for the fold changes we use the mean of the NRA suggested enzyme fold changes for each enzyme in a design across the 13 kinetic models.
- We only run this once for all 3 enzymes together. The models all perform well in this case, obviating the need for testing enzyme individually
- The results are used for figure S.13B
- Each run consisting of 10x10x5 ODE simulations takes around 4 hours.
