#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyStateSpectrum import getSteadyStateSpectrum
#
def gammaA_sweep_spec(N_em,g,gamma2,gammaA_logmin,gammaA_logmax,kappa,pump):
    #gammaD-list
    dloggA=0.1
    loggA=np.arange(gammaA_logmin,gammaA_logmax+dloggA,dloggA)
    gammaA_List=10**(loggA)
    #perform pump-sweep
    Spectrum_matrix=[[]]*len(gammaA_List)
    #get N_Hilbert from prior test
    out_ss=np.load('./data/{}-emitter_gammaA-sweep_ss_NH_intelligent_g={}_gam2={}_kap={}_pump={}.npz'.format(N_em,g,gamma2,kappa,pump))
    N_Hilb_list=out_ss['N_H_list']
    
    for i,gamma in enumerate(gammaA_List):      
        print(i+1,' out of ',len(gammaA_List))#output to show how far the simulation is
        N_Hilbert=N_Hilb_list[i]
        wlist, Spectrum = getSteadyStateSpectrum(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
        Spectrum_matrix[i]=Spectrum/np.max(Spectrum)
    #save
    np.savez('./data/{}-emitter_gammaA-sweep_spectrum_g={}_gam2={}_kap={}_pump={}.npz'.format(N_em,g,gamma2,kappa,pump),w=wlist,spec=Spectrum_matrix,gammaA_List=gammaA_List)
