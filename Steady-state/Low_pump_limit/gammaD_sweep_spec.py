#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyStateSpectrum import getSteadyStateSpectrum
#
def gammaD_sweep_spec(N_em,g,kappa,gammaD_logmin,gammaD_logmax,gamma,pump):
    #gammaD-list
    dloggD=0.1
    loggD=np.arange(gammaD_logmin,gammaD_logmax+dloggD,dloggD)
    gammaD_List=10**(loggD)
    #perform pump-sweep
    Spectrum_matrix=[[]]*len(gammaD_List)
    #get N_Hilbert from prior test
    out_ss=np.load('./data/{}-emitter_gammaD-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump))
    N_Hilb_list=out_ss['N_H_list']
    
    for i,gammaD in enumerate(gammaD_List):
        gamma2=gammaD/2        
        print(i+1,' out of ',len(gammaD_List))#output to show how far the simulation is
        N_Hilbert=N_Hilb_list[i]
        wlist, Spectrum = getSteadyStateSpectrum(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
        Spectrum_matrix[i]=Spectrum/np.max(Spectrum)
    #save
    np.savez('./data/{}-emitter_gammaD-sweep_spectrum_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump),w=wlist,spec=Spectrum_matrix,gammaD_List=gammaD_List)
