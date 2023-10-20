In this folder you will find the following files that are required for the NOMAD workflow
** tfa_samples 
- This is a csv file with around 5000 computationally generated steady state profiles
- These were generated using pyTFA and contain the following variables
-- concentrations
-- fluxes
-- Standard Gibb's free energy for all reactions
-- Thermodynamic displacement for all reactions 
-- Sample IDs 1618, 3717, 3970, 4011, 4012, 4015, 4024, 4026, 4027, 4028, 4029, 4057, 4693 are the steady-state profiles around which the
13 enhanced kinetic models were parametrized

** allosteric_regulations
- This file contains a list of compiled allosteric regulations that include mixed, simple, and competitive inhibition.
- For mixed inhibition, the same inhibitor, reaction_id pair appears twice, once with simple and once with competitive inhibition.

** enhanced_kinetic_models
- This is an hdf5 file containing all the kinetic parameters corresponding to each of the 13 enhanced kinetic models
- The parameters are the kms, kis (inhibition constants), maximal velocities (Vmaxs) and equilibrium constants.

** general_rxn_subsystem
- This is a file with the different reactions and the subsystems they belong to.
- It is used when generating figure 4.

There is also a subfolder, fermentation-experiments, that contains all the fermentation data reported by Hernandez et al. (2009).
The curves were obtained using webplotdigitizer.



