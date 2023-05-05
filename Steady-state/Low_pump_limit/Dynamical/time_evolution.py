from cmath import sqrt
from qutip import *
import numpy as np
#import implementation of solving single-emitter
from solveMasterEq_Nemitter import solveMasterEq_Nemitter

def time_evolution(N_em,g,kappa,pump,gamma,gamma2,Tw,Nt,N_Hilbert):
    #create time-list
    tlist=np.linspace(0.0, Tw, Nt)#time list [ps]

    result=solveMasterEq_Nemitter(N_em=N_em, g=g, kappa=kappa,pump=pump, gamma=gamma,gamma2=gamma2, tlist=tlist, N_Hilbert=N_Hilbert)
    c_occ=result.expect[0]
    e_occ_Matrix=[[]]*N_em
    for i in range(N_em):
        e_occ_Matrix[i]=result.expect[i+2]
    #save
    np.savez('.\data\{}-emitter_time_evolution_pump={}_Tw={}_Nt={}_NH={}.npz'.format(N_em,pump,Tw,Nt,N_Hilbert),c_occ=c_occ,e_occM=e_occ_Matrix,tlist=tlist)

