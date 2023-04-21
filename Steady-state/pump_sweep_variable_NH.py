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
    neMatrix=np.zeros((len(pumpList),N_em))
    #initialize intelligent loop
    N_Hilbert=10
    inc=2
    N_Hilbert_list=[[]]*len(pumpList)
    cav_capacity_occ_list=[[]]*len(pumpList)
    g_state=basis(2,0)
    e_state=basis(2,1)

    for i,pump in enumerate(pumpList):
        print(i+1,' out of ',len(pumpList))#output to show how far the simulation is
        #obtain steady state density operator and operators
        rho_ss, c_occ, c_occ_p2, emitter_occ_list = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)

        #calculate cavity population of highest included state
        a_max_state_cav=basis(N_Hilbert,N_Hilbert-1)
        a_max_state_em=tensor(g_state+e_state)
        for i_em in range(1,N_em):
            a_max_state_em=tensor(a_max_state_em,tensor(g_state+e_state))
        aMax=tensor(a_max_state_em,a_max_state_cav)
        print((aMax.dag()*rho_ss*aMax).shape)
        cav_capacity_occ=(aMax.dag()*rho_ss*aMax).tr()

        #find smallest N_Hilbert satisfying our convergence criterion
        maxtol_cav_occ=0.0002
        cav_occ_Optimized=False
        if cav_capacity_occ>maxtol_cav_occ:
            while cav_capacity_occ>maxtol_cav_occ:
                N_Hilbert+=inc
                rho_ss, c_occ, c_occ_p2, emitter_occ_list = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
                a_max_state_cav=basis(N_Hilbert,N_Hilbert-1)
                aMax=tensor(a_max_state_em,a_max_state_cav)
                cav_capacity_occ=(aMax.dag()*rho_ss*aMax).tr()
                print('N_Hilbert={}'.format(N_Hilbert))
                N_Hilbert_list[i]=N_Hilbert
            cav_occ_Optimized=True
        elif 10<N_Hilbert-inc and cav_occ_Optimized==False:
            while 10<N_Hilbert and cav_capacity_occ<maxtol_cav_occ:
                N_Hilbert+=-inc
                rho_ss, c_occ, c_occ_p2, emitter_occ_list = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
                a_max_state_cav=basis(N_Hilbert,N_Hilbert-1)
                aMax=tensor(a_max_state_em,a_max_state_cav)        
                cav_capacity_occ=(aMax.dag()*rho_ss*aMax).tr()
                print('N_Hilbert={}'.format(N_Hilbert))
                N_Hilbert_list[i]=N_Hilbert
            #went one step too far
            N_Hilbert+=inc
            rho_ss, c_occ, c_occ_p2, emitter_occ_list = getSteadyState(N_em=N_em, g=g, kappa=kappa, pump=pump, gamma=gamma, gamma2=gamma2, N_Hilbert=N_Hilbert)
            a_max_state_cav=basis(N_Hilbert,N_Hilbert-1)
            aMax=tensor(a_max_state_em,a_max_state_cav) 
            cav_capacity_occ=(aMax.dag()*rho_ss*aMax).tr()
            print('N_Hilbert={}'.format(N_Hilbert))
            N_Hilbert_list[i]=N_Hilbert
        #optimal state reached
        #save relevant quantitites
        nP=(c_occ*rho_ss).tr()
        nP_p2=(c_occ_p2*rho_ss).tr()
        g2=(nP_p2-nP)/nP**2
        npList[i]=nP
        g2List[i]=g2
        for i_em in range(N_em):
            e_occ=emitter_occ_list[i_em]
            ne=(e_occ*rho_ss).tr()
            neMatrix[i,i_em]=ne
        cav_capacity_occ_list[i]=cav_capacity_occ
        N_Hilbert_list[i]=N_Hilbert
    #save
    np.savez('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2),g2=g2List,nP=npList,pump_over_g_save=pump_over_g_List,N_H_list=N_Hilbert_list,cav_conv=cav_capacity_occ,neMatrix=neMatrix)
    #tests
    #print(cav_capacity_occ_list)
    #print(max(cav_capacity_occ_list))
    #print(neMatrix[20:25,:])