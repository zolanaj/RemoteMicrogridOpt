# Remote Microgrid Optimization project repository

## About
This is a public repository containing the python code and working files in support of the following paper:
Zolan, A., Scioletti, M., Morton, D., and Newman, A., Decomposing Mixed-integer Programs for Optimal Microgrid Design (2019).  INFORMS Journal on Computing, accepted.

A link to the paper will be included as soon as a first look of the paper is available.  For now, if you want to use the paper and these models for reserach or analysis, please use the citation above.  If you want to use the data contained in the case studies, please use the following citation: 

M. Engels, P. A. Boyd, T. M. Koehler, S. Goel, D. R. Sisk, D. D. Hatley, V. V. Mendon, and J. C. Hail. Smart and Green Energy (SAGE) for base camps final report. Technical report, Pacific Northwest National Laboratory (PNNL), Richland, WA (US), 2014

A link to the above paper is located here: https://www.osti.gov/biblio/1160200

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

## General use of the model

Due to the limits of the multiprocessing library in interactive Python as of the time of this publication, it is best to run the scripts in Python 2 from the command line.  The module to run to obtain meaningful results is RunPH.py, and after navigating to the directory ```path_to_repository/RemoteMicrogridOpt/python``` an example command takes the form 

```python RunPH.py scen m n model minlp_flag```

in which: 

- **scen** is the string "ll" followed by an integer between 1 and 14 that denotes the location index;
- **m** is a positive integer that denotes the number of partitions on battery state of charge, with valid entries being positive integers; 
- **n** is a positive integer that denotes the number of partitions on battery current; 
- **model** denotes the partitioning model to solve. 
- **minlp_flag** solves the MINLP nonlinear model to within the optimality criterion if the argument is **minlp**, and solves the MILP approximation otherwise.
- **mingen_flag** does not use a mincap generator valid inequality if the argument is equal to **none** and uses it otherwise

Keywords for **model** are as follows (see the paper for model descriptions):


|keyword  | Model |
| -------- | ----- | 
| univariate | Univariate implementation of Model (U) |
| bivariate | Bivariate implementation of Model (U) |
| gounaris | Univariate implementation of Model (G) |
| nagarajan | Bivariate implementation of Model (N) |
| -------- | ----- | 

Example: ```python RunPH.py ll12 1 4 univariate minlp```

## Replicating analyses from the paper

Using the section on general use of the model above as a guide, the results and figures from the above may be obtained using the arguments below.  The outputs for each case are a results file with some basic outputs on the model, and an iteration file that records lower and upper bounds on the optimal objective value of a model instance as new bounds are obtained.  The former may be used to build all the tables, while the latter is useful in developing figures 6 and 7 from the paper.

### Table 1, Algorithm 1
Column For Algorithm 1: Let **m=1, n=1, model=univariate, minlp_flag=n, mingen_flag=y** for all **scen** from **ll1** to **ll14**
Column "Solve Directly": instead of the command above, use the command: ```python Partition_Univariate.py scen```, for all **scen** from **ll1** to **ll14**

### Table 2
Column (U): Let **m=1, n=4, model=univariate, minlp_flag=n, mingen_flag=y** for all **scen** from **ll1** to **ll14**
Column (G): Let **m=1, n=4, model=gounaris, minlp_flag=n, mingen_flag=y** for all **scen** from **ll1** to **ll14**
Column (N): Let **m=1, n=4, model=nagarajan, minlp_flag=n, mingen_flag=y** for all **scen** from **ll1** to **ll14**

### Table 3
Let **m=1, n=4, model=univariate, minlp_flag=minlp, mingen_flag=y** for all **scen** from **ll1** to **ll14**

### Figures 6 and 7
Let **m,n** correspond to each instance in the table, **model=univariate, minlp_flag=minlp, mingen_flag=y** if **n=1** and **model=bivariate** otherwise, for all **scen** from **ll1** to **ll14**

### Table 4
First Column: Let **m=1, n=4, model=univariate, minlp_flag=minlp, mingen_flag=none** for all **scen** from **ll1** to **ll14**
Second Column: Let **m=1, n=4, model=univariate, minlp_flag=minlp, mingen_flag=y** for all **scen** from **ll1** to **ll14**
