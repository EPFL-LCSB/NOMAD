**** NOTE: You can directly run scripts II and IV as all the samples are provided.
**** Rerunning I and III will result in a different set of samples

I_kinetic_parameter_sampling.py
- This script generates kinetic parameter sets that are consistent with the mean across all the generated steady state profiles (profile 3191)
- It first samples 15,000 parameter sets and only outputs that have dominant time constants that are ~3x faster than the doubling time of
ST10284 in fed-batch fermentation conditions.

II_kinetic_parameter_analysis.py
- This is to address some reviewer comments
- We study the variation in KMs and Vmaxs acorss the 9 kinetic models to see how much phenotypic coverage there is.
- We also compare the Km distribution across the kinetic models with the reported Km ranges in BRENDA.

III_kinetic_parameter_sampling_vmax_sensitivity.py
- For the vmax sensitivity study in Supplementary Note 1
- Generates 2,000 models around 100 steady states that are generated by 3x all the exofluxomics from the fed-batch data,
3x all the intracellular fluxes as well +-20%, and allowing the concentrations to vary within 20% of their original steady state built around the fed-batch info
- The objective of this study is to see how far 'off' we are when we simply 3x all the Vmaxs to simulate batch conditions from fed-batch models, knowing that
 the growth rate in batch conditions is roughly 3x that of the fed-batch conditions.

IV_vmax_sensitivity_analysis.py
- This script uses the results from scripts 1 and 2 to study how far the Vmaxs across the 2,000 perturbed kinetic models vary from 3x the Vmax of one of the original models

** The parameters folder contains all the models with good linearised dynamics, along with the top 9 models with their Vmaxs multiplied by a factor of 3


