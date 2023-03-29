# compute spectrum for 
from cmath import sqrt
from qutip import *
import numpy as np

def compute2op2t_spectrum_trapz(N_em,g,kappa, pump, gamma,gamma2,N,Tw,N_Hilbert):
    #operators (emitter1, ...,emitterN_em, cavity)-basis
    #a destroys a photon in the cavity
    a=tensor(identity(2))
    for i in range(1,N_em):
        a=tensor(a,identity(2))
    a=tensor(a, destroy(N_Hilbert))#destroys a photon in the cavity
    #build Hamiltonian and collapse operators
    #initialize
    C_list=[None]*(1+3*N_em)
    C_kappa=np.sqrt(kappa)*a
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
        C_gamma2=np.sqrt(gamma2*2)*sm.dag()*sm
        C_list[3*i+1]=C_pump
        C_list[3*i+2]=C_gamma
        C_list[3*i+3]=C_gamma2

    #initial state operator
    g_state=basis(2,0)
    e_state=basis(2,1)
    c1=basis(N_Hilbert,0)
    #c2=basis(N_Hilbert,1)
    #construct state where emitter 1 is occupied
    rho0=e_state
    for i in range(1,N_em):
        rho0=tensor(rho0,g_state)
    rho0=tensor(rho0,c1)

    #create tlist
    dt=Tw/N
    #dw=2*np.pi/Tw
    #Ww=2*np.pi/dt
    tlist = np.arange(0,Tw,dt)
    #wlist=dw*np.arange(-N/2,N/2,1)
    Ww=0.8
    Nw=800
    Ww=0.6
    #Nw=600
    dw=Ww/Nw
    wlist=np.arange(-Ww/2,Ww/2+dw,dw)
    

    #calculate 2-time correlation function using QuTiP
    A=a.dag()
    B=a
    C2=correlation_2op_2t(H,rho0,tlist,tlist,C_list,A,B,reverse=False,solver='me')

    Int1=np.trapz(C2,dx=dt,axis=0)
    phase=np.outer(wlist,tlist)
    Integr2=Int1*np.exp(-1j*phase)
    Int2=np.trapz(Integr2,dx=dt,axis=1)
    S=kappa*Int2
    print(S[:5])
    return wlist, S