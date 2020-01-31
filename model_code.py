import os
import stochpy
import random
import numpy as numpy
from scipy import stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--suffix", help="suffix for files")
args = parser.parse_args()
suffix = args.suffix

workingdir = os.getcwd()

# Simulation parameters
start_time = 0.0
end_time = 8760
nruns=10000
cases = numpy.zeros([nruns, 5])

# Run is a single run of the model that returns the number of incident cases
# This outputs the parameter change notice every iteration
def Baserun(iterations):
    for k in range(0, iterations):
        pdict = {'epsilon1':random.uniform(0.015,0.417),
             'epsilon2':random.uniform(0.1,96)
            } 
        model = stochpy.SSA()
        model.Model(model_file='protective_decol.psc', dir=workingdir)
        model.ChangeParameter('epsilon1',pdict['epsilon1'])
        model.Endtime(end_time)
        model.ChangeParameter('epsilon2',1/pdict['epsilon2'])
        model.DoStochSim()
        model.GetRegularGrid(n_samples=end_time)
        outcomes = model.data_stochsim_grid.species
        cases[k,0] = outcomes[16][0][-1]
        cases[k,1] = outcomes[34][0][-1]
        cases[k,2] = outcomes[40][0][-1]
        cases[k,3] = pdict['epsilon1']
        cases[k,4] = pdict['epsilon2']
    return cases
  
#Fit delta
base_sweep = Baserun(nruns)
print("Complete")

print("Saving Files")
numpy.savetxt('decol_simulations.csv',base_sweep,delimiter=',',header="Colonization Cases,Bacteremia Cases, Bacteremia Related Deaths,CHG Replacement Effectiveness,Duration of Protection",comments='')