#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
###
#pretty plot
from matplotlib import rcParams
rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 16})
rcParams['axes.titlepad'] = 20

###PARAMETERS
#system
N_em=1
g=2 #coupling constant [THz]
gamma2=0/np.sqrt(2) #decay rate [THz] for coupling from cavity to environment
gammaD=np.sqrt(2)*gamma2
gamma=0
pump=0.0001


###
fig, ax1 = plt.subplots(1,figsize=(11,7))

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_kappa-sweep_spectrum_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
wlist=out_ss['w']
spectrum_matrix=out_ss['spec']
out_ss=np.load('./data/{}-emitter_kappa-sweep_ss_NH_intelligent_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
kappa_list=out_ss['kappa_List']


##ax3
ax1.set(xlabel=r'$\mathrm{log}(\kappa)\, [\mathrm{ps^{-1}}]$') 
ax1.set(ylabel=r'$(\omega-\omega_{eg})\,\mathrm{[ps^{-1}]}$')
#ax3.set(ylim=[-1.5*g,1.5*g])
#ax3.set_xticks(ticks=np.arange(-3,4,1),labels=10^(np.arange(-3,4,1)))
#g=np.sqrt(N_em)*g
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(2)-1)),(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(3)-np.sqrt(2))),(spectrum_matrix).transpose(), 20, cmap='RdGy')
ax1.contourf(np.log10(kappa_list), wlist,(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss)[np.log10(pg_list_ss)<0], wlist,(spectrum_matrix[np.log10(pg_list_ss)<0]).transpose(), 20, cmap='RdGy')
plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=ax1)
#ax1.axes.xaxis.set_ticklabels([])

#augmentation
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,gammaD,gamma,pump)
ax1.set(title=TITLE)
#ax1.set(xlim=[min(pump_list),max(pump_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
#ax1.set(ylim=[min(linewidth)/5,max([max(linewidth)*1.3])])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#save plot
plt.savefig('./plots/spectrumplot_{}-emitter_kappa-sweep_g={}_gamD={}_gam={}_pump={}.pdf'.format(N_em,g,gammaD,gamma,pump))
plt.show()