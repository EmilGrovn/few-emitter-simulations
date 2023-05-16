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
kappa=0.1
pump=0.001

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_gammaA-sweep_linewidth_g={}_gam2={}_kap={}_pump={}.npz'.format(N_em,g,gamma2,kappa,pump))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('./data/{}-emitter_gammaA-sweep_ss_NH_intelligent_g={}_gam2={}_kap={}_pump={}.npz'.format(N_em,g,gamma2,kappa,pump))
nP_list_ss=out_ss['nP']
g2_list_ss=out_ss['g2']
neMatrix=out_ss['neMatrix']
gammaA_list=out_ss['gammaA_List']
N_Hilbert=N_em+1

###
fig, ax1 = plt.subplots(1,figsize=(11,7))
#analytiske linjebredder
gammaR=4*g**2/(pump+gammaA_list+gammaD+kappa)
Reg1_linewidth=(3*pump+gammaA_list+gammaD+kappa)/2
Reg3_linewidth=gammaA_list/(nP_list_ss+1)
neList=neMatrix[:,0]
ST_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em


##plot linewidths
ST_unmodified=ST_linewidth*2
ax1.loglog(gammaA_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
ax1.loglog(gammaA_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
ax1.loglog(gammaA_list, ST_linewidth/(2*np.pi),'g:')    #plot ST (regime 2) linewidth
ax1.loglog(gammaA_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax1.loglog(gammaA_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth
kappa_list=kappa*np.ones(len(gammaA_list))
ax1.loglog(gammaA_list, kappa*np.ones(len(gammaA_list))/(2*np.pi),'b--')  #simple approx
#ax1.loglog(gammaA_list, linewidth[0]/(2*np.pi)+(linewidth[27]-linewidth[0])/(kappa_list[27]-kappa_list[0])*(kappa_list-kappa_list[0])/(2*np.pi),'y--')
kap1=kappa_list**(-1)
#ax1.loglog(kappa_list, linewidth[-20]/(2*np.pi)+(linewidth[-20]-linewidth[-1])/(kap1[-20]-kap1[-1])*(kap1-kap1[-20])/(2*np.pi),'y--')

#ax1.loglog(kappa_list, linewidth[min(kappa_list-1000)==(kappa_list-1000)]/(2*np.pi)-(kappa_list-kappa_list[min(kappa_list-1000)==(kappa_list-1000)])/(2*np.pi),'m--')

#augmentation
ax1.legend(['Simulation','ST unmodified','ST modified','Rabi peak','Quenching','$\kappa$','$\propto \kappa$','$\propto 1/\kappa$'])
ax1.grid()
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$. $\kappa$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,gammaD,kappa,pump)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\gamma_A\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(gammaA_list),max(gammaA_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_gammaA-sweep_g={}_gamD={}_kap={}_pump={}.pdf'.format(N_em,g,gammaD,kappa,pump))
#plt.show()
print(linewidth)
dw=wlist[1]-wlist[0]
print(linewidth/dw)

#plt.plot(kappa_list,g2_list_ss)
#plt.xscale('log')
plt.show()
print(g2_list_ss)
print(nP_list_ss*10**4)
print(neMatrix)