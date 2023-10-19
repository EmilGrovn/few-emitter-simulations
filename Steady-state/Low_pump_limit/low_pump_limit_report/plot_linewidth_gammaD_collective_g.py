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
N_em=4
g=2 #coupling constant [THz]
gCol=g*np.sqrt(N_em)
kappa=8.5 #decay rate [THz] for coupling from cavity to environment
gamma=0.1
gammaA=gamma
pump=0.0001

###IMPORT DATA
out_ss=np.load('../data/{}-emitter_gammaD-sweep_linewidth_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('../data/{}-emitter_gammaD-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_pump={}.npz'.format(N_em,g,kappa,gamma,pump))
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
gammaD_list=out_ss['gammaD_List']
N_Hilbert=N_em+1
try:
    out_ss=np.load('./data/spectrum_formula_linewidth_g={}_kap={}_gamA={}.npz'.format(gCol,kappa,gammaA))
    linewidth_list=out_ss['linewidth_list']
except NameError:
    print('No analytical spectrum prediction:(')

###
fig, ax1 = plt.subplots(1,figsize=(11,7))
#analytical prediction
gammaR=4*g**2/(pump+gamma+gammaD_list+kappa)
Rabi_broadening_linewidth=(3*pump+gamma+gammaD_list+kappa)/2   #rabi_broadening
Quenching_linewidth=kappa/(nP_list_ss+1)
neList=neMatrix[:,0]
ST_modified_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em
ST_unmodified_linewidth=ST_modified_linewidth*2

##plot linewidths
#simulation
ax1.loglog(gammaD_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
#analytical linewidths from literature
ax1.loglog(gammaD_list, ST_unmodified_linewidth/(2*np.pi),'g-')    #plot ST (regime 2) linewidth
ax1.loglog(gammaD_list, Rabi_broadening_linewidth/(2*np.pi),'r-')  #regime 1 linewidth
ax1.loglog(gammaD_list, Quenching_linewidth/(2*np.pi),'k-')  #regime 3 linewidth
#ax1.loglog(kappa_list, kappa_list/(2*np.pi),'b-')  #simple approx
ax1.loglog(gammaD_list, kappa/(2*np.pi)*np.ones(len(gammaD_list)),'b--')  #simple approx
#ax1.loglog(gammaD_list, (gammaA+gammaD_list)/(2*np.pi)*np.ones(len(gammaD_list)),'m--')

#construct analytic FWHM prediction
a0=-kappa/4
a1=-(gammaD_list+gammaA)/4
a2=a0-a1
b0=np.sqrt(a2**2-gCol**2)
FWHM_linewidth=np.zeros((len(gammaD_list)))
FWHM_linewidth_approx=np.zeros((len(gammaD_list)))
for i in range(len(gammaD_list)):
    a2i=a2[i]
    a3i=a0+a1[i]
    b0i=b0[i]
    if abs(a2i)<gCol:
        FWHM_linewidth[i]=(kappa+gammaD_list[i]+gammaA)/2
    elif a2i>gCol:
        #gammaA and gammaD dominates
        FWHM_linewidth[i]=2*np.sqrt(-a3i**2-b0i**2+np.sqrt(2*a3i**4+2*b0i**4))
        FWHM_linewidth_approx[i]=kappa+4*gCol**2/(gammaD_list[i]+gammaA-kappa)
    elif -a2i>gCol:
        #kappa dominates
        FWHM_linewidth[i]=2*np.sqrt(-a3i**2-b0i**2+np.sqrt(2*a3i**4+2*b0i**4))
        FWHM_linewidth_approx[i]=gammaD_list[i]+gammaA+4*gCol**2/(kappa-gammaD_list[i]-gammaA)   
#predictions of our model
ax1.loglog(gammaD_list,linewidth_list/(2*np.pi),'k.')
ax1.loglog(gammaD_list,FWHM_linewidth/(2*np.pi),'c.')
ax1.loglog(gammaD_list,FWHM_linewidth_approx/(2*np.pi),'m.')

#augmentation
ax1.legend(['Simulation','ST unmodified','Rabi peak','Quenching','$\kappa$','Full spectrum prediction','P=0 analytical','P=0 analytical large $|a_2|$'])
#ax1.legend(['Simulation','ST unmodified','Rabi peak','Quenching','$\kappa$','$\gamma_A+\gamma_D$','Full spectrum prediction','P=0 analytical','P=0 analytical large $a_2$'])
ax1.grid()
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\kappa$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,kappa,gamma,pump)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\gamma_D\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(xlim=[min(gammaD_list),max(gammaD_list)])
ax1.set(ylim=[min(linewidth)/(2*np.pi)/2,max(linewidth)/(2*np.pi)*2])

#mark different regimes
Q12=-4*gCol+(-gammaA+kappa)
Q23=4*gCol+(-gammaA+kappa)
Q2len=(Q23-Q12)
print(Q12)

Min=ax1.get_ylim()[0]
Max=ax1.get_ylim()[1]
plt.semilogx([Q12,Q12],[Min,Max],'k--')
plt.text(Q12/10, Max/1.5, '$-a_2>g\sqrt{N_{em}}$')
plt.semilogx([Q23,Q23],[Min,Max],'k--')
plt.text(Q12+Q2len/10, Max/1.5, '$|a_2|<g\sqrt{N_{em}}$')
plt.text(Q23*10, Max/1.5, '$a_2>g\sqrt{N_{em}}$')

#save plot
plt.savefig('./plots/linewidthplot_{}-emitter_gammaD-sweep_g={}_kap={}_gam={}_pump={}.pdf'.format(N_em,g,kappa,gamma,pump))
plt.show()