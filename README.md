# Remote Microgrid Optimization project repository

## About
This is a public repository containing the python code and working files in support of the following paper:
Zolan, A., Scioletti, M., Morton, D., and Newman, A., Decomposing Mixed-integer Programs for Optimal Microgrid Design (2019).  INFORMS Journal on Computing, accepted.

A link to the paper will be included as soon as a first look of the paper is available.

The code was written in the CPLEX/Python API, using Python 2.7.9.  Additional python libraries, like SciPy, are used here, and so a Python distribution like Anaconda is recommended.  Anaconda distributions for Python 2 are available here: 
https://docs.anaconda.com/anaconda/install/hashes/win-2-64/

CPLEX distributions may be obtained via the IBM website (free for academics):
https://www.ibm.com/analytics/cplex-optimizer  

Please note, the intended purpose of this repository is to share the code with researchers and allow users to replicate the results from the IJOC paper.  The code itself has room for improvement in efficiency and in-line documentation, and the latter will be a focus of future updates.

## Structure

There are three folders within this repository: 
 - 'inputs' -- contains all files read as input by the model
 - 'python' -- contains all python code required to run the analyses performed in the paper
 - 'sol_files' -- contains output files, some of which are read by the program in runtime
 
## Replicating analyses from the paper

### General use of the model

Due to the limits of the multiprocessing library in interactive Python as of the time of this publication, it is best to run the scripts in Python 2 from the command line.  The module to run to obtain meaningful results is RunPH.py, and after navigating to the directory **path_to_repository/RemoteMicrogridOpt/python**, an example command takes the form 

**python RunPH.py scen m n model**

in which: 

- **scen** is the string "ll" followed by an integer between 1 and 14 that denotes the location index;
- **m** is a positive integer that denotes the number of partitions on battery state of charge, with valid entries being positive integers; 
- **n** is a positive integer that denotes the number of partitions on battery current; and,
- **model** denotes the partitioning model to solve. keywords for model are as follows: 

========================
keyword  | Model
========================
univariate | Univariate implementation of Model (U)
bivariate | Bivariate implementation of Model (U)
gounaris | Univariate implemnetation of Model (G)
nagarajan | Univariate implementation of Model (N)
========================

Example: **python Run

### 