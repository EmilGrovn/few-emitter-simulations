#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyStateSpectrum import getSteadyStateSpectrum
#
def kappa_sweep_spec(N_em,g,gamma2,kappa_logmin,kappa_logmax,gamma,pump):
    #gammaD-list
    dlogkap=0.1
    logkap=np.arange(kappa_logmin,kappa_logmax+dlogkap,dlogkap)
    kappa_List=10**(logkap)
    #perform pump-sweep
    Spectrum_matrix=[[]]*len(kappa_List)
    #get N_Hilbert from prior test
    out_ss=np.load('./data/{}-emitter_kappa-sweep_ss_NH_intelligent_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
    N_Hilb_list=out_ss['N_H_list']
    
    for i,kappa in enumerate(kappa_List):      
        print(i+1,' out of ',len(kappa_List))#output to show how far the simulation is
        N_Hilbert=N_Hilb_list[i]
        wlist, Spectrum = getSteadyStateSpectrum(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
        Spectrum_matrix[i]=Spectrum/np.max(Spectrum)
    #save
    np.savez('./data/{}-emitter_kappa-sweep_spectrum_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump),w=wlist,spec=Spectrum_matrix,kappa_List=kappa_List)
