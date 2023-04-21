#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyState import getSteadyState
#
def pump_sweep_fixed_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2,N_Hilbert):
    #pump-list
    dlogP_over_g=0.1
    logP_over_g=np.arange(pump_logmin,pump_logmax+dlogP_over_g,dlogP_over_g)
    pump_over_g_List=10**(logP_over_g)
    pumpList=pump_over_g_List*g
    #perform pump-sweep
    npList=np.zeros((len(pumpList)))
    g2List=np.zeros(len(pumpList))
    for i,pump in enumerate(pumpList):
        print(i)
        rho_ss, c_occ, c_occ_p2, emitter_occ_list = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma,gamma2=gamma2, N_Hilbert=N_Hilbert)
        nP=(c_occ*rho_ss).tr()
        nP_p2=(c_occ_p2*rho_ss).tr()
        g2=(nP_p2-nP)/nP**2
        npList[i]=nP
        g2List[i]=g2
    #save
    np.savez('.\data\{}-emitter_pump-sweep_ss_NH={}_g={}.npz'.format(N_em,N_Hilbert,g),g2=g2List,nP=npList,pump_over_g_save=pump_over_g_List)
