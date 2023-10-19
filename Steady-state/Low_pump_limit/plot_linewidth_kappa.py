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
N_em=2
g=2 #coupling constant [THz]
#gamma2=0/np.sqrt(2) #decay rate [THz] for coupling from cavity to environment
g=2 #coupling constant [THz]
gammaD=2
gamma2=gammaD/2 #decay rate [THz] for coupling from cavity to environment
gamma=1 #background decay [THz]
gammaA=gamma
# gammaD=0
# gamma2=gammaD/2 #decay rate [THz] for coupling from cavity to environment
# gamma=11
# gammaA=gamma
# gammaD=0
# gamma2=gammaD/2 #decay rate [THz] for coupling from cavity to environment
# gamma=11
# gammaA=gamma
pump=0.0001

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_kappa-sweep_linewidth_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('./data/{}-emitter_kappa-sweep_ss_NH_intelligent_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
nP_list_ss=out_ss['nP']
g2_list_ss=out_ss['g2']
neMatrix=out_ss['neMatrix']
kappa_list=out_ss['kappa_List']
N_Hilbert=N_em+1
try:
    out_ss=np.load('./data/spectrum_formula_linewidth_g={}_gamD={}_gamA={}.npz'.format(g,gammaD,gammaA))
    linewidth_list=out_ss['linewidth_list']
except:
    print('Hey:)')

###
fig, ax1 = plt.subplots(1,figsize=(11,7))
#analytiske linjebredder
gammaR=4*g**2/(pump+gamma+gammaD+kappa_list)
Reg1_linewidth=(3*pump+gamma+gammaD+kappa_list)/2
Reg3_linewidth=kappa_list/(nP_list_ss+1)
neList=neMatrix[:,0]
ST_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em


##plot linewidths
ST_unmodified=ST_linewidth*2
print(kappa_list[min(abs(kappa_list-10))==abs(kappa_list-10)])
print('linewidth at kappa=1000',linewidth[min(abs(kappa_list-1000))==abs(kappa_list-1000)]/(2*np.pi))
ax1.loglog(kappa_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
ax1.loglog(kappa_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
ax1.loglog(kappa_list, ST_linewidth/(2*np.pi),'g:')    #plot ST (regime 2) linewidth
ax1.loglog(kappa_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax1.loglog(kappa_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth
ax1.loglog(kappa_list, kappa_list/(2*np.pi),'b--')  #simple approx
#ax1.loglog(kappa_list, linewidth[0]/(2*np.pi)+(linewidth[27]-linewidth[0])/(kappa_list[27]-kappa_list[0])*(kappa_list-kappa_list[0])/(2*np.pi),'y--')
kap1=kappa_list**(-1)
#ax1.loglog(kappa_list, linewidth[-20]/(2*np.pi)+(linewidth[-20]-linewidth[-1])/(kap1[-20]-kap1[-1])*(kap1-kap1[-20])/(2*np.pi),'y--')
###ax1.loglog(kappa_list[kappa_list<4*g],kappa_list[kappa_list<4*g]/2/(2*np.pi),'y:')
#construct analytic FWHM prediction
a0=-kappa_list/4
a1=-(gammaD+gammaA)/4
a2=a0-a1
b0=np.sqrt(a2**2-g**2)
FWHM_linewidth=np.zeros((len(kappa_list)))
FWHM_linewidth_approx=np.zeros((len(kappa_list)))
FWHM_extended=np.zeros((len(kappa_list)))
for i in range(len(kappa_list)):
    a2i=a2[i]
    if abs(a2i)<g:
        FWHM_linewidth[i]=(kappa_list[i]+gammaD+gammaA)/2
        FWHM_linewidth_approx[i]=(kappa_list[i]+gammaD+gammaA)/2
        FWHM_extended[i]=(kappa_list[i]+gammaD+gammaA)/2
    elif a2i>g:
        #gammaA and gammaD dominates
        FWHM_linewidth[i]=abs(2*a0[i]*(1+np.sqrt(1-g**2/a2i**2))+2*a1*(1-np.sqrt(1-g**2/a2i**2)))
        FWHM_linewidth_approx[i]=kappa_list[i]+4*g**2/(gammaD+gammaD-kappa_list[i])
        FWHM_extended[i]=np.sqrt(-b0[i]**2+ np.sqrt(4*b0[i]**4+4*((a0[i]+a1)**2-b0[i]**2)**2))
    elif -a2i>g:
        #gammaA and gammaD dominates
        FWHM_linewidth[i]=abs(2*a0[i]*(1-np.sqrt(1-g**2/a2i**2))+2*a1*(1+np.sqrt(1-g**2/a2i**2)))
        FWHM_linewidth_approx[i]=gammaD+gammaA+4*g**2/(kappa_list[i]-gammaD-gammaA)
        FWHM_extended[i]=np.sqrt(-b0[i]**2+ np.sqrt(4*b0[i]**4+4*((a0[i]+a1)**2-b0[i]**2)**2))   

#abs af spektrum
#FWHM_linewidth=np.sqrt(3)*FWHM_linewidth
#a2=a0-a1
#sign_a2=a2/np.abs(a2)
#FWHM_linewidth=2*a0+2*a1+2*np.sqrt(a2**2-g**2)
#The
#FWHM_linewidth=(gammaD+gammaA)/2*(1-np.sqrt(1-16*g**2/(kappa_list-(gammaD+gammaA))**2))+kappa_list/2*(1+np.sqrt(1-16*g**2/(kappa_list-(gammaA+gammaD))**2))
#FWHM_linewidth=kappa_list/2*(1+np.sqrt(1-16*g**2/(kappa_list-(gammaA+gammaD))**2))
#FWHM_linewidth=2*a0*(1-sign_a2*np.sqrt(1-g**2/a2**2))+2*a1*(1+sign_a2*np.sqrt(1-g**2/a2**2))
#ax1.loglog(kappa_list[kappa_list>=4*g],kappa_list[kappa_list>=4*g]/2*(1-np.sqrt(1-(4*g/kappa_list[kappa_list>=4*g])**2))/(2*np.pi),'y--')
###ax1.loglog(kappa_list[kappa_list>=4*g],kappa_list[kappa_list>=4*g]/2*(1-np.sqrt(1-(4*g/kappa_list[kappa_list>=4*g])**2))/(2*np.pi),'y--')
ax1.loglog(kappa_list,FWHM_linewidth/(2*np.pi),'m.')
ax1.loglog(kappa_list,FWHM_linewidth_approx/(2*np.pi),'y*')
try:
    ax1.loglog(kappa_list,linewidth_list/(2*np.pi),'k.')
except NameError:
    print('Hola')
###ax1.loglog(kappa_list,FWHM_extended/(2*np.pi),'b.')
#ax1.loglog(kappa_list, linewidth[min(kappa_list-1000)==(kappa_list-1000)]/(2*np.pi)-(kappa_list-kappa_list[min(kappa_list-1000)==(kappa_list-1000)])/(2*np.pi),'m--')

#augmentation
#ax1.legend(['Simulation','ST unmodified','ST modified','Rabi peak','Quenching','$\kappa$','$\propto \kappa$','$\propto 1/\kappa$'])
###ax1.legend(['Simulation','ST unmodified','ST modified','Rabi peak','Quenching','$\kappa$','P=0 analytical strong','P=0 analytical weak','P=0 analytical large $\kappa$: $4g^2/\kappa$'])
ax1.legend(['Simulation','ST unmodified','ST modified','Rabi peak','Quenching','$\kappa$','P=0 analytical','P=0 analytical approx','full spectrum formula'])
ax1.grid()
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,gammaD,gamma,pump)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(kappa_list),max(kappa_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
#ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_kappa-sweep_g={}_gamD={}_gam={}_pump={}.pdf'.format(N_em,g,gammaD,gamma,pump))
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