Data used for building models and for plotting

steady_states.csv (available at https://doi.org/10.5281/zenodo.15432261)
- The set of 5,000 steady state profiles consistent with the strain ST10284 built using pyTFA
- The sample closest to the mean of all the profiles (index 3191) was used for kinetic model construction

steady_states_vmax_sensitivity.csv
- 100 sets of steady state profiles sampled for the vmax sensitivity study
- Details are in Supplementary Note I of the manuscript

ST10284_cons_relations_annotated.csv
- A matrix of conserved moieties. This is necessary so that we don't model dependent variables.
- For example, the total adenylate pool remains constant. Therefore if we model atp and adp, then amp is implicitly determined.

kms_in_our_model_manual_curation_.csv
- KMs pulled from BRENDA for S.cerevisiae.
- The KMS were then manually curated, with metabolite names changed to metabolite ids in our model.
- This data is used for comparing the KMs across our 9 kinetic models with those in BRENDA

reaction_ec_pairing.csv
- This is used in conjunction with kms_in_our_model_manual_curation_.csv
- We need it to pair the EC class outputted by BRENDA to our model reaction IDs.
- Note that a single EC can sometimes correspond to multiple reactions.

batch_experiments_biomass.csv / batch_experiments_pca.csv / batch_experiments_glucose.csv
- Experimental data from experiments conducted in triplicate
- This is for all the 10 design strains, and the reference strain - ST10284
- Contains the mean and standard deviation
- Used for plotting Figure 4 and Supplementary figures 4-6


