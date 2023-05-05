import numpy as np
from qutip import *

def solveMasterEq_Nemitter(N_em, g, kappa, pump, gamma, gamma2, tlist,N_Hilbert):
    #operators (emitter1, ...,emitterN_em, cavity)-basis
    #a destroys a photon in the cavity
    a=tensor(identity(2))
    for i in range(1,N_em):
        a=tensor(a,identity(2))
    a=tensor(a, destroy(N_Hilbert))#destroys a photon in the cavity
    #build Hamiltonian and collapse operators
    #initialize
    C_list=[None]*(1+3*N_em)
    C_kappa=np.sqrt(kappa)*a    #cavity leak
    C_list[0]=C_kappa
    #build
    for i in range(0,N_em):
        if i==0:
            sm=destroy(2)
        else:
            sm=identity(2)
        for j in range(1,N_em):
            if i == j:
                sm=tensor(sm,destroy(2))
            else:
                sm=tensor(sm,identity(2))
        sm=tensor(sm,identity(N_Hilbert))
        #Hamiltonian
        if i==0:
            H=g*(a.dag()*sm+a*sm.dag())
        else:
            H+=g*(a.dag()*sm+a*sm.dag())
        #collapse operators
        C_pump=np.sqrt(pump)*sm.dag()   #pump
        C_gamma=np.sqrt(gamma)*sm       #background decay from excited state
        C_gamma2=np.sqrt(gamma2*2)*sm.dag()*sm  #pure dephasing
        C_list[3*i+1]=C_pump
        C_list[3*i+2]=C_gamma
        C_list[3*i+3]=C_gamma2
    
    #initial state operator
    g_state=basis(2,0)
    e_state=basis(2,1)
    c1=basis(N_Hilbert,0)
    c2=basis(N_Hilbert,1)
    #construct state where emitter 1 is occupied
    rho0=e_state
    for i in range(1,N_em):
        rho0=tensor(rho0,g_state)
    rho0=tensor(rho0,c1)

    #operators we want to calculate
    operator_list=[None]*(2+N_em)
    c_occ=a.dag()*a#cavity occupied
    c_occ_power2=c_occ*c_occ
    operator_list[0]=c_occ
    operator_list[1]=c_occ_power2
    for i in range(N_em):
        measure=(-sigmaz()+1)/2
        if i==0:
            e_occ=tensor(measure)
        else:
            e_occ=tensor(identity(2))
        for j in range(1,N_em):
            if j==i:
                e_occ=tensor(e_occ,measure)
            else:
                e_occ=tensor(e_occ,identity(2))
        e_occ=tensor(e_occ,identity(N_Hilbert))
        operator_list[i+2]=e_occ

    #calculate
    result = mesolve(H, rho0, tlist, C_list, operator_list)
    return result