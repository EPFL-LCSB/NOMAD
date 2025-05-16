constraint_based_model_scaffold.json
- Used for all kinetic modelling and strain design codes
- All the omics constraints are already built in
- The thermo constraints were applied by using data from the equilibrator-api
- The omics intracellular omics constraints were applied using data from Park et al. (see references in main manuscript)
- Extracellular glucose, growth, and pca secretion rates were constrained using fed-batch data for ST10284

kinetic_model_scaffold.yaml
- Contains the scaffold model which will then be paired with sampled kinetic parameters
- Contains info on the reactions, metabolites, and reaction mechanisms in the model.