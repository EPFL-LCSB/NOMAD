# NOMAD
## About this repo
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

## Operating conditions
- To ensure compatibility across different operating systems, the codes are run inside a docker environment.
- You will need to install Docker https://www.docker.com/products/docker-desktop/ to build and run the codes in a self-contained environment.
- The container runs with python version 3.6
- You can either use the ready made Docker image or build the image yourself  
- The exact versions of the various dependencies and packages can be found in the requirements.txt file and the Dockerfile.
- For the MILP solver, we have used IBM CPLEX Studio 12.8.
- We have tested the codes within the docker environment on 5 different computers, two with Windows 10 Pro 64-bit, two operating Ubuntu 20.04.2 LTS, and a MacOS Ventura 13.5.2 (m2 chip)
- For the MacOS Ventura, we used the prebuilt image. For all the others, we built the image using the Dockerfile

## Downloading and preparing the repository
1. Clone the repository on your local machine by going to the desired folder and typing '''git clone https://github.com/EPFL-LCSB/NOMAD.git'''
2. Download and install docker https://www.docker.com/products/docker-desktop/
3. Download an appropriate solver such as CPLEX Studio 12.8
--> Copy the contents of the solver to NOMAD/docker/solvers or nomad/docker/solvers
--> Further instructions are in the solver folder

## Installation
### Building your own docker image
1. '''cd <base-directory>/nomad/docker'''
2. '''build.bat''' OR '''./build'''
3. Common issues:
-- If you are building the dockerfile on windows, make sure all the files in the utils directory have unix EOL (end of line)
-- Makes sure that the root user is added to the Dockergroup 
4. It takes around 15 minutes to build the docker on the above mentioned machines.
### Using the prebuilt image
You can also use the prebuilt image by following the instructions below:
1. pull the image using the follwing command: 
'''
docker pull bharathn1990/nomad_docker:final
'''
2. Change the run file, either run.bat in Windows or run in Unix by replacing "nomad_docker" with "bharathn1990/nomad_docker:final"

## Running the codes
1. Once you have built the image or prepped the ready made one, create a container using code(run.bat) or code(./run) 
2. Once you are in the container, go to the folder with all the scripts
--> cd ../../../NOMAD/study-1-K_trpD9923/scripts for the first study
--> cd ../../../NOMAD/study-2-eK_trpD9923/scripts for the second study
--> RUN!!!   
Run times for each of the scripts is provided in the readme file in the scripts folder. 

## Data provided with this repo
1. The parameters that characterize the 10 kinetic models for the first study, K_trpD9923, are provided in ./study-1-K_trpD9923/data/kinetic_params_top_10_models.hdf5
2. The final set of 41 unique designs from K_trpD9923 is in the csv file ./study-1-K_trpD9923/data/all_unique_designs.csv
3. The parameters for the 13 enhanced kinetic models for the second study, eK_trpD9923, are provided in ./study-2-K_trpD9923/data/enhanced_kinetic_models.hdf5
4. The final set of 34 unique designs from eK_trpD9923 is in the csv file ./study-2-K_trpD9923/output/data/all_unique_designs_eK_trpD9923.csv
5. The final set of 13 unique designs from eK_trpD9923_d2 is in the csv file ./study-2-K_trpD9923/output/data/all_unique_designs_eK_trpD9923_d2.csv
