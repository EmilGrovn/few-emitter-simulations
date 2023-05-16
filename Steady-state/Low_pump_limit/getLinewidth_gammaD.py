#calculate linewidth

import numpy as np
from qutip import *

def getLinewidth_gammaD(N_em, g, kappa, gamma, pump):
    ###IMPORT DATA
    out_ss=np.load('./data/{}-emitter_gammaD-sweep_spectrum_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump))
    wlist=out_ss['w']
    spectrum_matrix=out_ss['spec']
    pump_list_ss=out_ss['gammaD_List']
    ## calculate linewidth ##
    linewidth=np.zeros(len(pump_list_ss))
    for i in np.arange(0,len(pump_list_ss)):
        spec_i=spectrum_matrix[i,:]     #choose spectrum for given pump
        #use symmetry
        wlistPos=wlist[wlist>=0]
        specPos=spec_i[wlist>=0]
        #identify maximum
        specMax=max(specPos)
        indexVec=np.arange(len(specPos))
        maxIndex=indexVec[specMax==specPos]
        wMax=wlistPos[maxIndex]
        #initialize while-loop
        #check for increasing omega
        Dw=wlist[1]-wlist[0]
        wInc=wMax
        specInc=specMax
        currentIndex=maxIndex
        wInc2=wlistPos[currentIndex+1]
        specInc2=specPos[currentIndex+1]
        boundary=max(wlistPos)
        while specInc2 >=0.5 and specInc2<specInc and wInc2<boundary:
            currentIndex+=1
            wInc=wInc2
            wInc2=wlistPos[currentIndex+1]
            specInc=specInc2
            specInc2=specPos[currentIndex+1]
        #check for increasing omega
        wDec=wMax
        specDec=specMax
        currentIndex=maxIndex
        wDec2=wDec-Dw
        if wDec2>=0:
            specDec2=specPos[currentIndex-1]
            while specDec2 >=0.5 and specDec2<specDec and wDec2>0:
                currentIndex-=1
                wDec=wDec2
                wDec2=wlistPos[currentIndex-1]
                specDec=specDec2
                specDec2=specPos[currentIndex-1]
        linewidth[i]=2*max([abs(wInc-wMax),abs(wMax-wDec)])
        print(linewidth[i])


        # spec_i=spectrum_matrix[i,:]     #choose spectrum for given pump
        # wlistAboveHM=wlist[spec_i>=0.5] #extract frequencies where spectrum is above half max
        # specAboveHM=spec_i[spec_i>=0.5] #extract corresponding spectrum values

        # linewidth[i]=np.max(wlistAboveHM)-np.min(wlistAboveHM)  #simple case with one pulse at w=0

        # if len(specAboveHM[wlistAboveHM==0])==0:
        #     print('Oh')
        #     then there are two seperate peaks for negative and positive omega
        #     linewidth[i]=np.max(wlistAboveHM[wlistAboveHM>=0])-np.min(wlistAboveHM[wlistAboveHM>=0])
        #     I am not sure this is entirely correct because there might still be a contribution from the other peak
        # elif specAboveHM[wlistAboveHM==0]<1:
        #     print('Oh',i,pg_list_ss[i])
        #     #then there are two partly seperated peaks for negative and positive omega
        #     wlistAboveHMPosOmega=wlistAboveHM[wlistAboveHM>=0] #Positive frequencies
        #     specAboveHMPosOmega=specAboveHM[wlistAboveHM>=0]   #corresponding spectrum values
        #     linewidth[i]=2*(np.max(wlistAboveHMPosOmega)-wlistAboveHMPosOmega[np.max(specAboveHMPosOmega)==specAboveHMPosOmega])
    np.savez('./data/{}-emitter_gammaD-sweep_linewidth_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump),wlist=wlist,linewidth=linewidth)