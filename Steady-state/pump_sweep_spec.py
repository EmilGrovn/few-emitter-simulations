#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyStateSpectrum import getSteadyStateSpectrum
#
def pump_sweep_spec(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2,N_Hilbert):
    #pump-list
    dlogP_over_g=0.1
    logP_over_g=np.arange(pump_logmin,pump_logmax+dlogP_over_g,dlogP_over_g)
    pump_over_g_List=10**(logP_over_g)
    pumpList=pump_over_g_List*g
    #perform pump-sweep
    Spectrum_matrix=[[]]*len(pumpList)
    #get N_Hilbert from prior test
    tol=10
    out_ss=np.load('.\data\{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_tol={}.npz'.format(N_em,g,kappa,tol))
    N_Hilb_list=out_ss['N_H_list']
    
    for i,pump in enumerate(pumpList):
        print(i)
        N_Hilbert=N_Hilb_list[i]
        wlist, Spectrum = getSteadyStateSpectrum(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
        Spectrum_matrix[i]=Spectrum/np.max(Spectrum)
    #save
    np.savez('.\data\{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(N_em,g,kappa,gamma,gamma2,tol),w=wlist,spec=Spectrum_matrix,pump_over_g_save=pump_over_g_List)
