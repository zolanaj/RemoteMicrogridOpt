# Remote Microgrid Optimization project repository

## About
This is a public repository containing the python code and working files in support of the following paper:
Zolan, A., Scioletti, M., Morton, D., and Newman, A., Decomposing Mixed-integer Programs for Optimal Microgrid Design (2019).  INFORMS Journal on Computing, accepted.

The code was written in the CPLEX/Python API, using Python 2.7.9.  Additional python libraries, like SciPy, are used here, and so a Python distribution like Anaconda is recommended.  Anaconda distributions for Python 2 are available here: 
https://docs.anaconda.com/anaconda/install/hashes/win-2-64/

CPLEX distributions may be obtained via the IBM website (free for academics):
https://www.ibm.com/analytics/cplex-optimizer  

## Structure

There are three folders within this repository: 
 - 'inputs' -- contains all files read as input by the model
 - 'python' -- contains all python code required to run the analyses performed in the paper
 - 'sol_files' -- contains output files, some of which are read by the program in runtime
 


