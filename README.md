# NOMAD
** About this repo **
This repository contains all the codes and softwares required to reproduce the work presented in the paper entitled
"Rational strain design with minimal phenotype perturbation".
The work is a property of the Laboratory for Computational Systems Biotechnology at the EPFL in Lausanne.
The repository is a stand alone that uses a docker based installation to ensure that no other softwares/dependencies need to be installed.
Note that it relies on dev versions of other softwares developed by the LCSB, namely

1. SKiMPy 
-- (Weilandt, Daniel R., et al. "Symbolic Kinetic Models in Python (SKiMpy): Intuitive modeling of large-scale biological kinetic models." bioRxiv (2022).)
-- https://github.com/EPFL-LCSB/skimpy
2. pyTFA
-- (Salvy, Pierre, et al. "pyTFA and matTFA: a Python package and a Matlab toolbox for Thermodynamics-based Flux Analysis." Bioinformatics 35.1 (2019): 167-169.)
-- https://github.com/EPFL-LCSB/pytfa/tree/master/pytfa

** Operating conditions **
- To ensure compatibility across different operating systems, a docker container is provided with all the necessary dependencies installed 
during the building of the environment.
- You will need to install Docker https://www.docker.com/products/docker-desktop/ to build and run the codes in a self-contained environment.
- The container runs with python version 3.6
- The exact versions of the various dependencies and packages can be found in the requirements.txt file and the Dockerfile.
- For the MILP solver, we have used IBM CPLEX Studio 12.8.
- We have tested the codes within the docker environment on 4 different computers, two with Windows 10 Pro 64-bit,
 and two operating Ubuntu 20.04.2 LTS 

** Instructions for installation and running ** 
1. Clone the repository on your local machine
2. Download and install docker https://www.docker.com/products/docker-desktop/
3. Download an appropriate solver such as CPLEX Studio 12.8
--> Copy the contents of the solver to NOMAD/docker/solvers or nomad/docker/solvers
--> Further instructions are in the solver folder
4. Build your NOMAD docker container
--> cd <base-directory>/nomad/docker
--> build.bat OR ./build
5. It takes around 15 minutes to build the docker on the above mentioned machines.  
6. Run your docker container
--> run.bat or ./run
7. Once you are in the container, go to the folder with all the scripts
--> cd ../../../NOMAD/anthranilate-study/scripts
--> RUN!!!   
8. Run times for each of the scripts is provided in the readme file in the scripts folder. 

** Data provided with this repo **
1. The parameters that characterize the 10 chosen kinetic models are provided in ./anthranilate-study/data/kinetic_params_top_10_models.hdf5
2. The final set of 41 unique designs is in the csv file ./anthranilate-study/data/all_unique_designs.csv
