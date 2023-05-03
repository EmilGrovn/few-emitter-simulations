#calculate linewidth

import numpy as np
from qutip import *

def getLinewidth(N_em, g, kappa, gamma, gamma2):
    ###IMPORT DATA
    out_ss=np.load('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
    wlist=out_ss['w']
    spectrum_matrix=out_ss['spec']
    pg_list_ss=out_ss['pump_over_g_save']
    ## calculate linewidth ##
    linewidth=np.zeros(len(pg_list_ss))
    for i in np.arange(0,len(pg_list_ss)):
        spec_i=spectrum_matrix[i,:]     #choose spectrum for given pump
        wlistAboveHM=wlist[spec_i>=0.5] #extract frequencies where spectrum is above half max
        specAboveHM=spec_i[spec_i>=0.5] #extract corresponding spectrum values

        linewidth[i]=np.max(wlistAboveHM)-np.min(wlistAboveHM)  #simple case with one pulse at w=0

        if len(specAboveHM[wlistAboveHM==0])==0:
            print('Oh')
            #then there are two seperate peaks for negative and positive omega
            linewidth[i]=np.max(wlistAboveHM[wlistAboveHM>=0])-np.min(wlistAboveHM[wlistAboveHM>=0])
            #I am not sure this is entirely correct because there might still be a contribution from the other peak
        elif specAboveHM[wlistAboveHM==0]<1:
            print('Oh',i,pg_list_ss[i])
            #then there are two partly seperated peaks for negative and positive omega
            wlistAboveHMPosOmega=wlistAboveHM[wlistAboveHM>=0] #Positive frequencies
            specAboveHMPosOmega=specAboveHM[wlistAboveHM>=0]   #corresponding spectrum values
            linewidth[i]=2*(np.max(wlistAboveHMPosOmega)-wlistAboveHMPosOmega[np.max(specAboveHMPosOmega)==specAboveHMPosOmega])
    np.savez('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2),wlist=wlist,linewidth=linewidth)