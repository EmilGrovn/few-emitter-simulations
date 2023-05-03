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
g=0.1 #coupling constant [THz]
kappa=0.02 #decay rate [THz] for coupling from cavity to environment
gamma=0.012
gammaD=1
gamma2=gammaD/np.sqrt(2)

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
pg_list_ss=out_ss['pump_over_g_save']
pump_list=g*pg_list_ss

###
fig, ax1 = plt.subplots(1,figsize=(11,7))
#analytiske linjebredder
gammaR=4*g**2/(pump_list+gamma+gammaD+kappa)
Reg1_linewidth=(3*pump_list+gamma+gammaD+kappa)/2
Reg3_linewidth=kappa/(nP_list_ss+1)
neList=neMatrix[:,0]
ST_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em


##plot linewidths
ST_unmodified=ST_linewidth*2
ax1.loglog(pump_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
ax1.loglog(pump_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
ax1.loglog(pump_list, ST_linewidth/(2*np.pi),'g:')    #plot ST (regime 2) linewidth
ax1.loglog(pump_list, kappa*np.ones((len(pg_list_ss)))/(2*np.pi),'b')  #simple approx
ax1.loglog(pump_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax1.loglog(pump_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth

#augmentation
ax1.legend(['Simulation','ST unmodified','ST modified','$\kappa$','Regime 1','Regime 3'])
ax1.grid()
TITLE='{}-emitter. $\gamma_D$={:.2f} THz'.format(N_em,gammaD)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$P\,[\mathrm{THz}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(pump_list),max(pump_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_pump-sweep_g={}_kap={}_gam={}_gam2={}.pdf'.format(N_em,g,kappa,gamma,gamma2))
plt.show()