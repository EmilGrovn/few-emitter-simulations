#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyStateSpectrum import getSteadyStateSpectrum
#
def pump_sweep_spec(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2):
    #pump-list
    dlogP=0.1
    logP=np.arange(pump_logmin,pump_logmax+dlogP,dlogP)
    pumpList=10**(logP)
    #perform pump-sweep
    Spectrum_matrix=[[]]*len(pumpList)
    #get N_Hilbert from prior test
    out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
    N_Hilb_list=out_ss['N_H_list']
    
    for i,pump in enumerate(pumpList):
        print(i+1,' out of ',len(pumpList))#output to show how far the simulation is
        N_Hilbert=N_Hilb_list[i]
        wlist, Spectrum = getSteadyStateSpectrum(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
        Spectrum_matrix[i]=Spectrum/np.max(Spectrum)
    #save
    np.savez('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2),w=wlist,spec=Spectrum_matrix,pumpList=pumpList)
