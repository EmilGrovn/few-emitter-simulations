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
kappa=0.02 #decay rate [THz] for coupling from cavity to environment
gamma=0*0.012
gamma2=0/np.sqrt(2)
gammaD=np.sqrt(2)*gamma2

kappaList=[0.1,0.5,1.0,1.1,1.5,2.0,4.0,8.0,16.0]


###
fig, axs = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(11,7))
for i in range(3):
    for j in range(3):
        kappa=kappaList[i*3+j]
        ###IMPORT DATA
        out_ss=np.load('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
        wlist=out_ss['w']
        spectrum_matrix=out_ss['spec']
        out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
        nP_list_ss=out_ss['nP']
        neMatrix=out_ss['neMatrix']
        pg_list_ss=out_ss['pump_over_g_save']
        pump_list=g*pg_list_ss

        #axs[i,j].set(xlabel=r'$P [THz]$') 
        #axs[i,j].set(ylabel=r'$(\omega-\omega_{eg})$')
        axs[i,j].contourf(np.log10(pump_list), wlist,(spectrum_matrix).transpose(), 20, cmap='RdGy')
        plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=axs[i,j])
        #axs[i,j].axes.xaxis.set_ticklabels([])
        axs[i,j].text(0, 0, '$\kappa$={} THz'.format(kappa))
        

#augmentation
TITLE='{}-emitter. $\gamma_D$={:.2f} THz. $g$={:.2f} THz'.format(N_em,gammaD,g)
fig.suptitle(TITLE)
#fig.supxlabel(r'$P [THz]$')
#fig.supylabel(r'$(\omega-\omega_{eg})$')
#ax1.set(xlim=[min(pump_list),max(pump_list)])

#save plot
plt.savefig('./plots/spectrum_kappa_multiplot_{}-emitter_pump-sweep_g={}_gam={}_gam2={}.pdf'.format(N_em,g,gamma,gamma2))