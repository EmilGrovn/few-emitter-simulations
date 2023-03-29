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
#Jespers system
hbar=6.582*10**2 #micro eV picoseconds
N_em=1
g=50/hbar
gamma=1/hbar
gamma2=10/hbar*0
pump=0
Q_factor_list=[3000,8000,16000,32000]
omega_resonance=1880.487454 #estimeret vha plot i Jespers note
#numerical
N_Hilbert=2
scale=2
Nt=600*scale
Tw=200*scale

S_matrix=[[]]*len(Q_factor_list)
Legend=[[]]*len(Q_factor_list)
styles=['b--','r--','g--','k--']
styles=['b','r','g','k']
#
for i,Q_factor in enumerate(Q_factor_list):
    kappa=omega_resonance/Q_factor
    wlist, Spectrum = compute2op2t_spectrum_trapz(N_em,g,kappa, pump, gamma,gamma2,Nt,Tw,N_Hilbert)
    S= Spectrum.real
    S=S/max(S)
    S_matrix[i]=S
    #Legend[i]='g/$\kappa$={:.2f}'.format(g/kappa)
    Legend[i]='Q={:.0f}'.format(Q_factor)
#
fig = plt.figure(1,figsize=(6,3))
print(max(S))
#plt.semilogy(wlist/g, S)
for i in range(len(Q_factor_list)):
    plt.plot(wlist/g, S_matrix[i],styles[i])
#plt.xlim(0,1)
plt.xlabel(r'$\left(\omega-\omega_{\rm eg}\right)/g_c$') 
plt.ylabel('Spectrum [a.u.]')
plt.text(-0.3/g, 0.5, '$\hbar\gamma_2=${:.0f}$\mu$eV'.format(hbar*gamma2))
figureLabel='(b)'
if gamma2==0:
    figureLabel='(a)'
    plt.legend(Legend,prop={'size': 16})
    #plt.legend(Legend)
plt.text(-4,0.9,figureLabel, weight='bold')
#plt.title('g={}'.format(g))
plt.xticks(np.arange(-4, 4+1, 1.0))
plt.yticks(np.arange(0, 1+0.5, 0.5))
#plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('./plots_rep/spectrum_normalized_long_gamma2={}.pdf'.format(gamma2*hbar))
#plt.show()