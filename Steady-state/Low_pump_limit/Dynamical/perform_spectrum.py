#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from compute2op2t_spectrum_trapz import compute2op2t_spectrum_trapz
#
from matplotlib import rcParams
rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 20})
rcParams['axes.titlepad'] = 20
    #"font.size": 16})
    #system
N_em=1
pump=0
g=0.4 #coupling constant [THz]
kappa=0.1 #decay rate [THz] for coupling from cavity to environment
gammaA=0
gammaD=1000
gamma=gammaA
gamma2=gammaD/np.sqrt(2)
#numerical
N_Hilbert=N_em+1
Nt=2000
Tw=600 #[ps]
#
wlist, Spectrum = compute2op2t_spectrum_trapz(N_em,g,kappa, pump, gamma,gamma2,Nt,Tw,N_Hilbert)
S= Spectrum.real
S=S/max(S)
#get linewidth
#use symmetry
wlistPos=wlist[wlist>=0]
specPos=S[wlist>=0]
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
linewidth=2*max([abs(wInc-wMax),abs(wMax-wDec)])

np.savez('.\data\{}-emitter_spectrum_g={}_kap={}_gamA={}_gamD={}_Tw={}_Nt={}_NH={}.npz'.format(N_em,g,kappa,gammaA,gammaD,Tw,Nt,N_Hilbert),S=S,wlist=wlist,linewidth=linewidth)

# #
# fig = plt.figure(1,figsize=(6,3))
# print(max(S))
# #plt.semilogy(wlist/g, S)
# for i in range(len(Q_factor_list)):
#     plt.plot(wlist/g, S_matrix[i],styles[i])
# #plt.xlim(0,1)
# plt.xlabel(r'$\left(\omega-\omega_{\rm eg}\right)/g_c$') 
# plt.ylabel('Spectrum [a.u.]')
# plt.text(-0.3/g, 0.5, '$\hbar\gamma_2=${:.0f}$\mu$eV'.format(gamma2))
# figureLabel='(b)'
# if gamma2==0:
#     figureLabel='(a)'
#     plt.legend(Legend,prop={'size': 16})
#     #plt.legend(Legend)
# plt.text(-4,0.9,figureLabel, weight='bold')
# #plt.title('g={}'.format(g))
# plt.xticks(np.arange(-4, 4+1, 1.0))
# plt.yticks(np.arange(0, 1+0.5, 0.5))
# #plt.legend()
# plt.grid()
# plt.tight_layout()
# plt.savefig('./plots/spectrum_normalized_long_gamma2={}.pdf'.format(gamma2))
# #plt.show()