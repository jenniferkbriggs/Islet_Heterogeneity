# Islet_Heterogeneity
This repository contains sample code for studying islet heterogeneity. Main scripts include network/hub cell analysis, wave/phase analysis, and first responder analysis. Code is available in Matlab and Python. All code was written by Jennifer Briggs (2024) and Richard Benninger's lab at the University of Colorado Anschutz Medical Campus.

Note: The python code is executed under virtual environment for version control purposes. The best way to avoid bugs is to run it in the same virtual environment.
The first time you clone this code, run: *conda env create -f environment.yml* in the terminal.
Then, before every run type *conda activate ./envs* to activate the environment before you run the files.

## Examples 
Examples to run all scripts and an example calcium file are found under the folder "examples"

## Network Analysis in Matlab
To run general network analysis, you will use the function *RunNetworkAnalysis.m* 

If you want to determine the optimal threshold for getting a log-log plot, input calcium into *findoptRth*. Please note, if you want to have the threshold based on number of links, please email me and I will add this to the code!

## Network analysis in Python
### General Network Analysis: **Run_Network.py**.
The first section (lines 5-17) is title options for you to change. This is wheter you define whether you would like figures, where to save the file to. ]

Try typing python3 Run_Network.py in the terminal and you should have success running this code.


### Dependencies:
The code depends on the following python packages: csvkit, collection, networkx, numpy, matplotlib.pyplot, scipy, easygui, pandas, python-math, tk, pyvis. **If you get an error with that looks like: Exception has occurred: ModuleNotFoundError: no module named '<modulename>'. Then that module must be downloaded using either Conda or pip. To download the module, go to the terminal and type (for example) 'conda install -c conda-forge easygui'.




