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
kappa=0.1 #decay rate [THz] for coupling from cavity to environment
gamma=0
gamma2=0/np.sqrt(2)
gammaD=np.sqrt(2)*gamma2

#gammaDList=[0,1,5,10,30,45,60,120,200]
#gammaDList=[0,1,5,10,30,60,120,1000,5000]
gammaDList=[0,1,5,10,30,120,1000,5000,10000]

###
fig, axs = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(11,7))
for i in range(3):
    for j in range(3):
        gammaD=gammaDList[i*3+j]
        gamma2=gammaD/np.sqrt(2)
        print(gammaD)
        ###IMPORT DATA
        out_ss=np.load('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
        wlist=out_ss['w']
        spectrum_matrix=out_ss['spec']
        out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
        nP_list_ss=out_ss['nP']
        neMatrix=out_ss['neMatrix']
        pump_list=out_ss['pumpList']

        #axs[i,j].set(xlabel=r'$P [THz]$') 
        #axs[i,j].set(ylabel=r'$(\omega-\omega_{eg})$')
        axs[i,j].contourf(np.log10(pump_list), wlist,(spectrum_matrix).transpose(), 20, cmap='RdGy')
        #plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=axs[i,j])
        if i==2:
            axs[i,j].axes.set_xticks([-4,-3,-2,-1], labels=[r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$'])
        axs[i,j].text(min(np.log10(pump_list)), max(wlist)*1.1, '$\gamma_D$={} ps$^-$$^1$'.format(gammaD))
        

#augmentation
TITLE='{}-emitter. $\kappa$={:.2f} ps$^-$$^1$. $\gamma_A$={:.2f} ps$^-$$^1$. $g$={} ps$^-$$^1$'.format(N_em,kappa,gamma,g)
fig.suptitle(TITLE)
fig.supxlabel(r'$P$ [ps$^-$$^1$]')
fig.supylabel(r'$(\omega-\omega_{eg})\,\mathrm{[ps^{-1}]}$')
#ax1.set(xlim=[min(pump_list),max(pump_list)])

#save plot
plt.savefig('./plots/spectrum_gammaD_multiplot_{}-emitter_pump-sweep_g={}_kap={}_gam={}.pdf'.format(N_em,g,kappa,gamma))