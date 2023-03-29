#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyState import getSteadyState
#
def pump_sweep_variable_NH(N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2):
    #pump-list
    dlogP_over_g=0.1
    logP_over_g=np.arange(pump_logmin,pump_logmax+dlogP_over_g,dlogP_over_g)
    pump_over_g_List=10**(logP_over_g)
    pumpList=pump_over_g_List*g
    #perform pump-sweep
    npList=np.zeros((len(pumpList)))
    g2List=np.zeros(len(pumpList))
    #initialize intelligent loop
    N_Hilbert=10
    inc=2
    N_Hilbert_list=[[]]*len(pumpList)
    cav_capacity_occ=[[]]*len(pumpList)
    g_state=basis(2,0)
    e_state=basis(2,1)

    for i,pump in enumerate(pumpList):
        print(i)
        rho_ss, c_occ, c_occ_p2 = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
        nP=(c_occ*rho_ss).tr()
        N_Hilbert_list[i]=N_Hilbert
        tol=10
        if (2*nP+tol)>N_Hilbert:
            while (2*nP+tol)>N_Hilbert:
                N_Hilbert+=inc
                rho_ss, c_occ, c_occ_p2 = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
                nP=(c_occ*rho_ss).tr()
                print('N_Hilbert={}'.format(N_Hilbert))
                N_Hilbert_list[i]=N_Hilbert
        elif max([10,2*nP+tol])<N_Hilbert-inc:
            while max([10,2*nP+tol])<N_Hilbert-inc:
                N_Hilbert+=-inc
                rho_ss, c_occ, c_occ_p2 = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
                nP=(c_occ*rho_ss).tr()
                print('N_Hilbert={}'.format(N_Hilbert))
                N_Hilbert_list[i]=N_Hilbert
        nP_p2=(c_occ_p2*rho_ss).tr()
        g2=(nP_p2-nP)/nP**2
        npList[i]=nP
        #g2=(nP_p2-nP)/nP**2
        g2List[i]=g2
        if N_em==1:
            a_max_state=basis(N_Hilbert,N_Hilbert-1)
            aMax=tensor(g_state,a_max_state)+tensor(e_state,a_max_state)
            print((aMax.dag()*rho_ss*aMax).shape)
            cav_capacity_occ[i]=(aMax.dag()*rho_ss*aMax).tr()
    #save
    np.savez('.\data\{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(N_em,g,kappa,gamma,gamma2,tol),g2=g2List,nP=npList,pump_over_g_save=pump_over_g_List,N_H_list=N_Hilbert_list,cav_conv=cav_capacity_occ)
    print(cav_capacity_occ)