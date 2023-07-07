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
gammaA=gamma
pump=0.0001

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_gammaD-sweep_linewidth_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('./data/{}-emitter_gammaD-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump))
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
gammaD_list=out_ss['gammaD_List']
N_Hilbert=N_em+1

###
fig, ax1 = plt.subplots(1,figsize=(11,7))
#analytiske linjebredder
gammaR=4*g**2/(pump+gamma+gammaD_list+kappa)
Reg1_linewidth=(3*pump+gamma+gammaD_list+kappa)/2
Reg3_linewidth=kappa/(nP_list_ss+1)
neList=neMatrix[:,0]
ST_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em
#analytical linewidth
gammaD_list=gammaD_list
FWHM_strong=(kappa+(gammaD_list+gammaA))/2
FWHM_linewidth=(gammaD_list+gammaA)/2*(1-np.sqrt(1-16*g**2/(kappa-(gammaD_list+gammaA))**2))+kappa/2*(1+np.sqrt(1-16*g**2/(kappa-(gammaA+gammaD_list))**2))
#FWHM_linewidth=gamma+gammaD_list+4*g**2/(kappa-(gammaD_list+gamma))
#FWHM_linewidth=kappa-4*g**2/(kappa-(gammaD_list+gamma))
##plot linewidths
ST_unmodified=ST_linewidth*2
ax1.loglog(gammaD_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
ax1.loglog(gammaD_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
ax1.loglog(gammaD_list, ST_linewidth/(2*np.pi),'g:')    #plot ST (regime 2) linewidth
ax1.loglog(gammaD_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax1.loglog(gammaD_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth
ax1.loglog(gammaD_list, kappa*np.ones((len(gammaD_list)))/(2*np.pi),'b--')  #simple approx
ax1.loglog(gammaD_list, FWHM_strong/(2*np.pi),'m--')  #simple approx
ax1.loglog(gammaD_list, FWHM_linewidth/(2*np.pi),'m--')  #simple approx

#augmentation
ax1.legend(['Simulation','ST unmodified','ST modified','Rabi peak','Quenching','$\kappa$','P=0 analytical'])
ax1.grid()
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\kappa$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,kappa,gamma,pump)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\gamma_D\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(gammaD_list),max(gammaD_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_gammaD-sweep_g={}_kap={}_gam={}_pump={}.pdf'.format(N_em,g,kappa,gamma,pump))
plt.show()
print(linewidth)
dw=wlist[1]-wlist[0]
print(linewidth/dw)

g2_list_ss=out_ss['g2']
plt.plot(gammaD_list,g2_list_ss)
#plt.show()
print(g2_list_ss)
print(nP_list_ss*10**4)
print(neMatrix*10**4)