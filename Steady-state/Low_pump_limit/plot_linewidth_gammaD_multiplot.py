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

gammaDList=[0,1,5,10,30,45,60,120,200]
gammaDList=[0,1,5,10,30,60,120,1000,5000]
gammaDList=[0,1,5,10,30,120,1000,5000,10000]

###
fig, axs = plt.subplots(3, 3, sharex=True, sharey=False, figsize=(11,7))
for i in range(3):
    for j in range(3):
        gammaD=gammaDList[i*3+j]
        gamma2=gammaD/np.sqrt(2)
        print(gammaD)
        ###IMPORT DATA
        out_ss=np.load('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
        wlist=out_ss['wlist']
        linewidth=out_ss['linewidth']
        out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
        nP_list_ss=out_ss['nP']
        neMatrix=out_ss['neMatrix']
        pump_list=out_ss['pumpList']

        #analytiske linjebredder
        gammaR=4*g**2/(pump_list+gamma+gammaD+kappa)
        Reg1_linewidth=(3*pump_list+gamma+gammaD+kappa)/2
        Reg3_linewidth=kappa/(nP_list_ss+1)
        neList=neMatrix[:,0]
        ST_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em

        ##plot linewidths
        ST_unmodified=ST_linewidth*2
        axs[i,j].loglog(pump_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
        axs[i,j].loglog(pump_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
        axs[i,j].loglog(pump_list, ST_linewidth/(2*np.pi),'g:')    #plot ST (regime 2) linewidth
        axs[i,j].loglog(pump_list, kappa*np.ones((len(pump_list)))/(2*np.pi),'b')  #simple approx
        if min(Reg1_linewidth/(2*np.pi))/1000<max(linewidth/(2*np.pi)):
            axs[i,j].loglog(pump_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
        axs[i,j].loglog(pump_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth

        #axs[i,j].set(xlabel=r'$P [THz]$') 
        #axs[i,j].set(ylabel=r'$(\omega-\omega_{eg})$')
        #axs[i,j].contourf(np.log10(pump_list), wlist,(spectrum_matrix).transpose(), 20, cmap='RdGy')
        #plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=axs[i,j])
        if i==2:
            axs[i,j].axes.set_xticks([0.0001,0.001,0.01,0.1], labels=[r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$'])
        if j==0:
            pass
        Min=min(linewidth)/(2*np.pi)
        Max=max(linewidth)*1.1/(2*np.pi)
        #axs[i,j].axes.set_ylim(Min,Max)
        axs[i,j].text(0.0001, axs[i,j].axes.get_ylim()[1]+(axs[i,j].axes.get_ylim()[1]-axs[i,j].axes.get_ylim()[0])*0.01, '$\gamma_D$={} ps$^-$$^1$'.format(gammaD))
        

#augmentation
TITLE='{}-emitter. $\kappa$={:.2f} ps$^-$$^1$. $\gamma_A$={:.2f} ps$^-$$^1$. $g$={} ps$^-$$^1$'.format(N_em,kappa,gamma,g)
print(TITLE)
fig.suptitle(TITLE)
fig.supxlabel(r'$P$ [ps$^-$$^1$]')
fig.supylabel(r'FWHM [THz]')
#ax1.set(xlim=[min(pump_list),max(pump_list)])

#save plot
plt.savefig('./plots/linewidth_gammaD_multiplot_{}-emitter_pump-sweep_g={}_kap={}_gam={}.pdf'.format(N_em,g,kappa,gamma))