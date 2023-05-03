#import
#from cmath import sqrt
#from qutip import *
import numpy as np
#import matplotlib.pyplot as plt
#import implementation of solving N-emitter
#from pump_sweep_fixed_NH import pump_sweep_fixed_NH
from pump_sweep_variable_NH import pump_sweep_variable_NH
from pump_sweep_spec import pump_sweep_spec
from getLinewidth import getLinewidth


###parameters
#system
N_em=1  #number of emitters
g=2 #coupling constant [THz]
kappa=0.02 #decay rate [THz] for coupling from cavity to environment
gamma=0.012*0 #background decay [THz]
gamma2=0/np.sqrt(2) #pure dephasing [THz]
#g=0.3 #coupling constant [THz]
#kappa=0.1 #decay rate [THz] for coupling from cavity to environment
#gamma=0.01 #background decay [THz]
#gamma2=0 #pure dephasing [THz]
pump_logmin=-3  #minimum pump/g=10^(pump_logmin)
pump_logmax=0   #maximum pump/g=10^(pump_logmax)
#numerical
N_Hilbert=10

###perform pump-sweep
#pump_sweep_fixed_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2,N_Hilbert)#get nP and g(0)^(2)
#for gammaD in [0,0.01,0.1,0.5,1]:
#    gamma2=gammaD/np.sqrt(2)
#[0.1,0.5,1.0,1.1,1.5,2.0,4.0,8.0,16.0]:
for kappa in [2.0,4.0,8.0,16.0]:
    pump_sweep_variable_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2)
    pump_sweep_spec(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2)#get spectrum
    getLinewidth(N_em, g, kappa, gamma, gamma2)