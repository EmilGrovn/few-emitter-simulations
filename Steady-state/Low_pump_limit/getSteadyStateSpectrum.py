import numpy as np
from qutip import *

def getSteadyStateSpectrum(N_em, g, kappa, pump, gamma, gamma2, N_Hilbert):
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
    #list of frequencies
    #create tlist
    N=20000
    #Ww=2/0.3*g
    Ww=3*g
    #Ww=3*g*np.sqrt(N_em)
    dw=Ww/N
    #tlist = np.arange(0,Tw,dt)
    wlist=dw*np.arange(-N/2,N/2,1)
    #wlist=dw*np.arange(0,N,1)

    Spectrum=spectrum(H,wlist,c_ops=C_list,a_op=a.dag(),b_op=a)
    #Spectrum=spectrum(H,wlist,c_ops=[],a_op=a.dag(),b_op=a)
    return wlist, Spectrum