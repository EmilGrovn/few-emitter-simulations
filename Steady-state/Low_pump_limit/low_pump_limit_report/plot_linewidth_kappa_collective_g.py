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
gCol=g*np.sqrt(N_em)
gammaD=8
gamma2=gammaD/2
gamma=5
gammaA=gamma
pump=0.0001

###IMPORT DATA
out_ss=np.load('../data/{}-emitter_kappa-sweep_linewidth_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('../data/{}-emitter_kappa-sweep_ss_NH_intelligent_g={}_gam2={}_gam={}_pump={}.npz'.format(N_em,g,gamma2,gamma,pump))
nP_list_ss=out_ss['nP']
g2_list_ss=out_ss['g2']
neMatrix=out_ss['neMatrix']
kappa_list=out_ss['kappa_List']
N_Hilbert=N_em+1
try:
    out_ss=np.load('./data/spectrum_formula_linewidth_g={}_gamD={}_gamA={}.npz'.format(gCol,gammaD,gammaA))
    linewidth_list=out_ss['linewidth_list']
except NameError:
    print('No analytical spectrum prediction:(')

###
fig, ax1 = plt.subplots(1,figsize=(11,7))
#analytical prediction
gammaR=4*g**2/(pump+gamma+gammaD+kappa_list)
Rabi_broadening_linewidth=(3*pump+gamma+gammaD+kappa_list)/2   #rabi_broadening
Quenching_linewidth=kappa_list/(nP_list_ss+1)
neList=neMatrix[:,0]
ST_modified_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em
ST_unmodified_linewidth=ST_modified_linewidth*2

##plot linewidths
print(kappa_list[min(abs(kappa_list-10))==abs(kappa_list-10)])
print('linewidth at kappa=1000',linewidth[min(abs(kappa_list-1000))==abs(kappa_list-1000)]/(2*np.pi))
#simulation
ax1.loglog(kappa_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
#analytical linewidths from literature
ax1.loglog(kappa_list, ST_unmodified_linewidth/(2*np.pi),'g-')    #plot ST (regime 2) linewidth
ax1.loglog(kappa_list, Rabi_broadening_linewidth/(2*np.pi),'r-')  #regime 1 linewidth
ax1.loglog(kappa_list, Quenching_linewidth/(2*np.pi),'k-')  #regime 3 linewidth
#ax1.loglog(kappa_list, kappa_list/(2*np.pi),'b-')  #simple approx
ax1.loglog(kappa_list, (gammaA+gammaD)/(2*np.pi)*np.ones(len(kappa_list)),'b--')  #simple approx
#construct analytic FWHM prediction
a0=-kappa_list/4
a1=-(gammaD+gammaA)/4
a2=a0-a1
b0=np.sqrt(a2**2-gCol**2)
FWHM_linewidth=np.zeros((len(kappa_list)))
FWHM_linewidth_approx=np.zeros((len(kappa_list)))
for i in range(len(kappa_list)):
    a2i=a2[i]
    l1=a0[i]+a1+b0[i]
    l2=a0[i]+a1-b0[i]
    a3i=a0[i]+a1
    b0i=b0[i]
    if abs(a2i)<gCol:
        FWHM_linewidth[i]=(kappa_list[i]+gammaD+gammaA)/2
        #FWHM_linewidth_approx[i]=(kappa_list[i]+gammaD+gammaA)/2
    elif a2i>gCol:
        #gammaA and gammaD dominates
        #FWHM_linewidth[i]=abs(2*a0[i]*(1+np.sqrt(1-g**2/a2i**2))+2*a1*(1-np.sqrt(1-g**2/a2i**2)))
        #FWHM_linewidth[i]=np.sqrt(2*(l1**2+l2**2)*np.sqrt(-1+np.sqrt(1+4*l1**2*l2**2/(l1**2+l2**2)**2)))
        #FWHM_linewidth[i]=np.sqrt(-2*l1**2-2*l2**2+2*np.sqrt(l1**4+l2**4+6*l1**2*l2**2))
        FWHM_linewidth[i]=2*np.sqrt(-a3i**2-b0i**2+np.sqrt(2*a3i**4+2*b0i**4))
        FWHM_linewidth_approx[i]=kappa_list[i]+4*gCol**2/(gammaA+gammaD-kappa_list[i])
    elif -a2i>gCol:
        #kappa dominates
        #FWHM_linewidth[i]=abs(2*a0[i]*(1-np.sqrt(1-g**2/a2i**2))+2*a1*(1+np.sqrt(1-g**2/a2i**2)))
        #FWHM_linewidth[i]=np.sqrt((l1**2+l2**2)/2*np.sqrt(-1+np.sqrt(1+4*l1**2*l2**2/(l1**2+l2**2)**2)))
        #FWHM_linewidth[i]=np.sqrt(-2*l1**2-2*l2**2+2*np.sqrt(l1**4+l2**4+6*l1**2*l2**2))
        FWHM_linewidth[i]=2*np.sqrt(-a3i**2-b0i**2+np.sqrt(2*a3i**4+2*b0i**4))
        FWHM_linewidth_approx[i]=gammaD+gammaA+4*gCol**2/(kappa_list[i]-gammaD-gammaA)   
#predictions of our model
ax1.loglog(kappa_list,linewidth_list/(2*np.pi),'k.')
ax1.loglog(kappa_list,FWHM_linewidth/(2*np.pi),'c.')
ax1.loglog(kappa_list,FWHM_linewidth_approx/(2*np.pi),'m.')

#augmentation
ax1.legend(['Simulation','ST unmodified','Rabi peak','Quenching','$\gamma_A+\gamma_D$','Full spectrum prediction','P=0 analytical','P=0 analytical large $|a_2|$'])
ax1.grid()
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,gammaD,gamma,pump)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(kappa_list),max(kappa_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
#ax1.axes.xaxis.set_ticklabels([])

#mark different regimes
Q12=-4*gCol+(gammaA+gammaD)
Q23=4*gCol+(gammaA+gammaD)
Q2len=(Q23-Q12)
print(Q12)

Min=ax1.get_ylim()[0]
Max=ax1.get_ylim()[1]
plt.semilogx([Q12,Q12],[Min,Max],'k--')
plt.text(Q12/10, Max/1.3, '$a_2>g\sqrt{N_{em}}$')
plt.semilogx([Q23,Q23],[Min,Max],'k--')
plt.text(Q12+Q2len/10, Max/1.3, '$|a_2|<g\sqrt{N_{em}}$')
#plt.text(Q23/10, Max/1.3, '$|a_2|<g$')
plt.text(Q23*10, Max/1.3, '$-a_2>g\sqrt{N_{em}}$')

# plt.loglog([Q12,Q12],[ax1.get_ylim()[0],ax1.get_ylim()[1]],'k--')
# plt.text(Q12/10, ax1.get_ylim()[1]/1.5, '$a_2>g$')
# plt.loglog([Q23,Q23],[ax1.get_ylim()[0],ax1.get_ylim()[1]],'k--')
# plt.text(Q23/10, ax1.get_ylim()[1]/1.5, '$|a_2|<g$')
# plt.text(Q23*10, ax1.get_ylim()[1]/1.5, '$-a_2>g$')

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_kappa-sweep_g={}_gamD={}_gam={}_pump={}.pdf'.format(N_em,g,gammaD,gamma,pump))
plt.show()
# print(linewidth)
# dw=wlist[1]-wlist[0]
# print(linewidth/dw)

# #plt.plot(kappa_list,g2_list_ss)
# #plt.xscale('log')
# plt.show()
# print(g2_list_ss)
# print(nP_list_ss*10**4)
# print(neMatrix)
print(gammaA+gammaD-linewidth_list)