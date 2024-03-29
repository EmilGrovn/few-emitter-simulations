#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from getSteadyState import getSteadyState
#
def gammaD_sweep_variable_NH(N_em,g,kappa,gammaD_logmin,gammaD_logmax,gamma,pump):
    #gammaD-list
    dloggD=0.1
    loggD=np.arange(gammaD_logmin,gammaD_logmax+dloggD,dloggD)
    gammaD_List=10**(loggD)
    #perform gammaD-sweep
    Npoints=len(gammaD_List)
    npList=np.zeros((Npoints))
    g2List=np.zeros(Npoints)
    neMatrix=np.zeros((Npoints,N_em))
    #initialize intelligent loop
    N0=2
    N_Hilbert=N0
    inc=1
    N_Hilbert_list=[[]]*Npoints
    cav_capacity_occ_list=[[]]*Npoints
    g_state=basis(2,0)
    e_state=basis(2,1)

    for i,gammaD in enumerate(gammaD_List):
        gamma2=gammaD/2
        print(i+1,' out of ',Npoints)#output to show how far the simulation is
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
        maxtol_cav_occ=0.0001   #10^(-4) chosen by judging from Matias' project and my project
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
        elif N0<N_Hilbert-inc and cav_occ_Optimized==False:
            while N0<N_Hilbert and cav_capacity_occ<maxtol_cav_occ:
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
    np.savez('./data/{}-emitter_gammaD-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump),g2=g2List,nP=npList,gammaD_List=gammaD_List,N_H_list=N_Hilbert_list,cav_conv=cav_capacity_occ,neMatrix=neMatrix)
    #tests
    #print(cav_capacity_occ_list)
    #print(max(cav_capacity_occ_list))
    #print(neMatrix[20:25,:])