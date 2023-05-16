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
g=0.4 #coupling constant [THz]
kappa=0.1 #decay rate [THz] for coupling from cavity to environment
gamma=0
gammaD=0
gamma2=gammaD/np.sqrt(2)

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
pump_list=out_ss['pumpList']
N_Hilbert=N_em+1
includeZeroPump=True
try:
    Nt=2000
    Tw=500 #[ps]
    out=np.load('.\Dynamical\data\{}-emitter_spectrum_g={}_kap={}_gamA={}_gamD={}_Tw={}_Nt={}_NH={}.npz'.format(N_em,g,kappa,gamma,gammaD,Tw,Nt,N_Hilbert))
    linewidth_dyn=out['linewidth']
except FileNotFoundError:
    print('FileNotFound', '.\Dynamical\data\{}-emitter_spectrum_g={}_kap={}_gamA={}_gamD={}_Tw={}_Nt={}_NH={}.npz'.format(N_em,g,kappa,int(gamma),gammaD,Tw,Nt,N_Hilbert))
    includeZeroPump=False

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
ax1.loglog(pump_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax1.loglog(pump_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth
ax1.loglog(pump_list, kappa*np.ones((len(pump_list)))/(2*np.pi),'b')  #simple approx
if includeZeroPump==True:
    ax1.loglog(pump_list, [linewidth_dyn/(2*np.pi)]*len(pump_list),'m')

#augmentation
ax1.legend(['Simulation','ST unmodified','ST modified','Rabi peak','Quenching','$\kappa$','P=0'])
ax1.grid()
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\kappa$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$'.format(N_em,g,kappa,gamma,gammaD)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$P\, \mathrm{[ps^{-1}]}$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(pump_list),max(pump_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_pump-sweep_g={}_kap={}_gam={}_gam2={}.pdf'.format(N_em,g,kappa,gamma,gamma2))
#plt.show()
print(linewidth)
dw=wlist[1]-wlist[0]
print(linewidth/dw)