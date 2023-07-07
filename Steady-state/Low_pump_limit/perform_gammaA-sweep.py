#import
#from cmath import sqrt
#from qutip import *
import numpy as np
#import matplotlib.pyplot as plt
#import implementation of solving N-emitter
#from pump_sweep_fixed_NH import pump_sweep_fixed_NH
from gammaA_sweep_variable_NH import gammaA_sweep_variable_NH
from gammaA_sweep_spec import gammaA_sweep_spec
from getLinewidth_gammaA import getLinewidth_gammaA


###parameters
#system
N_em=1  #number of emitters
g=2 #coupling constant [THz]
gammaD=0
gamma2=gammaD/2 #decay rate [THz] for coupling from cavity to environment
kappa=0.1 #background decay [THz]
#g=2 #coupling constant [THz]
#kappa=0.1 #decay rate [THz] for coupling from cavity to environment
#gamma=0.01 #background decay [THz]
#gamma2=0 #pure dephasing [THz]
gammaA_logmin=-2  #minimum pump/g=10^(pump_logmin)
gammaA_logmax=3   #maximum pump/g=10^(pump_logmax)
pump=0.00001
#numerical
#N_Hilbert=10

###perform pump-sweep
#pump_sweep_fixed_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2,N_Hilbert)#get nP and g(0)^(2)
#for kappa in [0.1,0.5,1.0,1.1,1.5,2.0,4.0,8.0,16.0]:
#for gammaD in [0,1,5,10,30,45,60,120,200]:
#    gamma2=gammaD/np.sqrt(2)
#for g in [1,1.2,1.5,1.8,2.0,2.5,3.0,4.0,5.0]:
#for g in [0.2,0.5,1.0,1.5,2.0,3.0,5.0,10.0,50.0]:
#for gamma in [0.1,0.5,1.0,1.1,1.5,2.0,4.0,8.0,16.0,32.0,64.0,128.0]:
#for gamma in [256.0]:
for g in [g]:
    gammaA_sweep_variable_NH(N_em,g,gamma2,gammaA_logmin,gammaA_logmax,kappa,pump)
    gammaA_sweep_spec(N_em,g,gamma2,gammaA_logmin,gammaA_logmax,kappa,pump)#get spectrum
    getLinewidth_gammaA(N_em, g, gamma2, kappa, pump)