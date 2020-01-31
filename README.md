# protective_decol_
## Code Repository for "Assessing the Potential Impact of a Long-acting Skin Disinfectant in the Prevention of Methicillin-resistant Staphylococcus aureus Transmission"


The Python scripts in this repository require the 'os', 'stochpy', 'numpy', 'stats', 'argparse', and 'random' libraries, all of which are available via pip and/or conda.
Note that a seed has not been specified in these files, so if the simulations are re-run, some changes in the results are to be expected due to stochastic variability

## Contents:

* 'protective_decol.psc': Configuration file for the stochastic simulation

* 'model_code.py': Python code for the running of the model itself and generating the data used in the analysis

* 'analysis.R' : A small R file to read the simulated data and fit a Poisson GLM to the outcomes of interest

* 'scatterplots.m' : MATLAB file to generate the scatter plots used in the analysis

* 'decol_simulations.csv' : The data used in the analysis. Note if this file is not renamed before running model_code.py the file will be overwritten