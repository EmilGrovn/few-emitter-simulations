#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from pump_sweep_fixed_NH import pump_sweep_fixed_NH
from pump_sweep_variable_NH import pump_sweep_variable_NH
from pump_sweep_spec import pump_sweep_spec
###
#pretty plot
from matplotlib import rcParams
# rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.sans-serif": ["Computer Modern Roman"],
#     "font.size": 16})
# rcParams['axes.titlepad'] = 20

###parameters
#system
N_em=2
g=0.1 #coupling constant [THz]
#g=0.3*np.sqrt(2)
#g=6
kappa=0.02 #decay rate [THz] for coupling from cavity to environment
pump=0
gamma=0.012
gamma2=1
pump_logmin=-3
pump_logmax=3
#numerical
N_Hilbert=10

###perform pump-sweep
#pump_sweep_fixed_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2,N_Hilbert)#get nP and g(0)^(2)
pump_sweep_variable_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2)
pump_sweep_spec(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2,N_Hilbert)#get spectrum