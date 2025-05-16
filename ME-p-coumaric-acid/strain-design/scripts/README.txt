I_strain_design.py
- Execute this first.
- Goes through each of the 9 kinetic models, and uses NRA to generate designs under the following constaints
-- objective: pca yield wrt glucose
-- max 3 enzyme modifications
-- max 2 fold change in enzyme activity
-- max 5 fold change in intracellular concentations
-- designs should provide increases in yield that are within 5% of the maximum possible increase in pca yield for that model

II_design_robustness.py
- Takes the unique designs across all the models (unique by membership)
- Applies them to each model in NRA, and calculates the predicted yield of pca on glucose for that design and model combination
-- for example, design 1 could be among the best designs for model 1. but when we enforce it on model 2, it might not be among the top 5% but top 20% instead.
-- However, there could be a particular design that is only 'good' for one model. i.e. produces 95% of the maximal yield for one model but for other models, it is not able to produce
   even 10% of the maximal yield
-- This outputs a matrix of solutions N_designs x N_models

III_nonlinear_verification.py
- Takes the design-model matrix from II and choses the top 10 designs based on their mean predicted increase in pca yield across all 9 models.
- Runs reactor simulations for each design model combination for the top 10 designs (90 simulations in total + 9 wild type simulations)
- You can either run it one model at a time or for all models. Note that it can take quite a long time for a single simulation, with times varying greatly between
    10 minutes and 4 hours!
- The results of this are used to plot Figures 2B, 2C and Supplementary figure S3.

IV_pca_yield_analysis.py
- Takes in the results of the ODE simulations, and identifies the increase in pca yield at the t=5 hour mark.
- We do this to check if the NRA predictions (which should be valid in the exponential growth phase) about the increase in pca yield hold true in a nonlinear context
- The results of this are used to plot the heatmap (plots/Figure_2A)

single_species.yaml
- Reactor starting info

** Plot **
The scripts in the plot folder are self explanatory and they can be run without running any of the above scripts as all
the results from the runs are provided.