This folder has the python scripts needed to generate the data that is presented in the paper.
The order in which the scripts are to be executed is given by the prefix attached to each file name.
Scripts starting with S. are for studies presented in the Supplementary notes.
The folder 'plot' contains the scripts to plot each figure.

** I_generate_all_designs.py
- In this script we generate all the NOMAD designs for a set of given design constraints using the 10 chosen kinetic models
- The results of this script are used in Figures 2, 3 and 5.
- For Figure 3, you need to run this script 4 times, with different allowable fold changes in concentration (2, 3, 10, and 20)
- The allowable fold change in concentration is set as a parameter inside the script
- For Figures 2 and 5, you only need the results for 3 fold change in concentrations.
- Each run for an allowable fold change in concentrations takes typically around 1 hour

** II_design_robustness
- In this script we go through all the designs that were generated for each kinetic model (within 5% of the maximal objective
for that model)
- We identify the unique designs my membership
- We then apply each of the unique designs to each of the kinetic models in an NRA setting (by setting the min. fold change 
in activity for that enzyme to 1e-3)
- We then store the solution value that NRA gives us.
- In this manner we obtain a matrix of 39 x 10 NRA solutions
- The results of this will be used to get the top 5 designs (by average predicted inc. in anthranilate yield as per NRA)
- The results will also be used to plot figures 4 and 6A
- A run takes around 30 minutes.

** III_get_fold_changes_top_designs
- This is a small script that just imposes the top 5 designs in all the models and gets the actual enzymatic fold changes 
suggested by NRA for each of those models
- These average proposed NRA fold changes are required to study the validity of the NRA predictions for the top 5 designs
in a nonlinear setting across different models and the sensitivity of the designs to enzymatic perturbations (scripts IV and V)
- A run takes around 30 minutes.

** IV_verification_of_top_designs_in_reactor
- We apply each of the top 5 designs to each kinetic model using the NRA suggested fold changes for that model and that design.
- We then simulate the response of the models to each design.
- This is used to plot Figures 7A and 7B.
- A run takes around 40 minutes.

** V_expression_sensitivity_of_top_designs
- Here we test the robustness of the designs to errors in experimental implementatoin
- we apply a 50% random perturbation 10 times to each combination of design and model
- Note that for the fold changes we use the mean of the NRA suggested enzyme fold changes for each enzyme in a design across the 10 kinetic models.
- We can apply to perturbation to either a single enzyme expression about its mean NRA predicted value while keeping the other 2 
at the mean NRA proposed values.
- To plot Figure 7C-F you need to run this script 4 times, for perturbing each enzyme individually and then all 3 together
- You can change which enzyme to perturb within the script.
- Each run consisting of 10x10x5 ODE simulations takes around 4 hours. 

** S.VI_phenotype_perturbation_analysis
- This is a sensitivity analysis test where we want to study how the NRA predicted anthranilate yield and the actual anthranilate production rate in nonlinear simulations
varies depending on how close / far we are from the reference phenotype.
- We vary the allowable fold changes in concentration from 2-fold to 20-fold
- We also vary the allowable fold changes in enzyme activity from 2-fold to 10-fold
- For each combination of enzyme activity and concentration fold changes, we generate the top design using NRA and simulate it in a bioreactor
- You can either run this script once or run it in batches on parallel CPUs
- For a given combination of concentration fold change and enzyme activity change, it takes around 20-30 minutes. Note the lower time, since we only simulate a single design (the top design).
- The results of this script are used for plotting figure S.9

** S.VII_double_mutant
- Two enzymes - DDPA and GLUDy - appear in all of our top 5 targets (see II_design_robustness and III_get_fold_changes_top_design)
- We wanted to know (prompted by a reviewer) what happens if we only target these two and how this double mutant performs compared to the triple mutant
- The results of the reactor simulations in this script as used to plot figure S.10

** S.VIII_i_aroGtktA_analysis_design_generation
- Prompted by reviews, we wanted to find out why our in-silico version of the experimentally implemented strain aroG^fbrtktA had such poor titers compared to the experimental implementation
- We use NRA to propose 4 additional modifications (apart from DDPA, TKT1 and TKT2) to increase anthranilate yield for each model
- We generated designs within 99% of the maximum predicted increase in anthranilate yield for each model
- The designs are then used by the next script for pruning

** S.VIII_ii_aroGtktA_analysis_design_collation
- After running the previous script, this script takes in all the generated designs, collates them and finds the most frequently occurring one.

** S.VIII_iii_aroGtktA_analysis_design_verification
- In this script, we verify the top design that occurs 6/10 times across the models.
- We do this in a reactor setting.
- The output of this script is used to generate figure S.11 in the supplementary material.
