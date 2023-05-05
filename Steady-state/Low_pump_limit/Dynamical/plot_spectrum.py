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
g=0.4 #coupling constant [THz]
kappa=0.1 #decay rate [THz] for coupling from cavity to environment
pump=0
gammaA=0
gammaD=0
gamma=gammaA
gamma2=gammaD/np.sqrt(2)
#numerical
N_Hilbert=N_em+1
Nt=2000
Tw=500 #[ps]
#

out=np.load('.\data\{}-emitter_spectrum_g={}_kap={}_gamA={}_gamD={}_Tw={}_Nt={}_NH={}.npz'.format(N_em,g,kappa,gammaA,gammaD,Tw,Nt,N_Hilbert))
wlist=out['wlist']
Spectrum=out['S']
linewidth=out['linewidth']
#
fig = plt.figure(1,figsize=(6,3))
#plt.semilogy(wlist/g, S)
plt.plot(wlist, Spectrum,'b*')
#plt.xlim(0,1)
plt.xlabel(r'$\left(\omega-\omega_{\rm eg}\right)$') 
plt.ylabel('Spectrum [a.u.]')
plt.xticks(np.arange(-1, 1+0.5, 0.5))
plt.yticks(np.arange(0, 1+0.25, 0.25))
print(linewidth)
plt.title('linewidth={:.5f} THz'.format(linewidth[0]/(2*np.pi)))
#plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('./plots/spectrum_normalized.pdf')
#plt.show()
print(linewidth/(2*np.pi))
dw=wlist[1]-wlist[0]
print(linewidth/dw)
print(dw)
print(dw/(2*np.pi))