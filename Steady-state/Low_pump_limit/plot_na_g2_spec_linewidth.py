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
gamma2=0.1/np.sqrt(2)
gammaD=np.sqrt(2)*gamma2


###IMPORT DATA
out_ss=np.load('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
wlist=out_ss['wlist']
linewidth=out_ss['linewidth']
out_ss=np.load('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
spectrum_matrix=out_ss['spec']
pump_list=out_ss['pumpList']
out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
g2_list_ss=out_ss['g2']
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
cav_conv_list=out_ss['cav_conv']
N_H_list=out_ss['N_H_list']

## identifying different regimes using g2 ##
delta=0.05
try:
    Pth1=min(pump_list[g2_list_ss>1-delta])
except ValueError:
    Pth1=max(pump_list)
try:
    Pth2=max(pump_list[g2_list_ss<1+delta])
except ValueError:
    Pth2=max(pump_list)
try:
    Pth_Qn=min(pump_list[g2_list_ss>2-delta])
except ValueError:
    Pth_Qn=max(pump_list)

##
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,figsize=(11,7))
#plot thresholds
#ax1.loglog([Pth1]*2, [min(nP_list_ss),max(nP_list_ss)],'k:')
#ax1.loglog([Pth2]*2, [min(nP_list_ss),max(nP_list_ss)],'k:')
ax2.loglog([Pth1]*2, [0.5,3],'k:')
ax2.loglog([Pth2]*2, [0.5,3],'k:')
#ax2.loglog([Pth_Qn]*2, [0.5,3],'k:')  #plot regime boundaries only on middle plot
#regimer
ax2.text(Pth1-Pth2*0.025, 1.4, 'I')
ax2.text(Pth1+Pth2*0.1, 1.4, 'II')
ax2.text(Pth_Qn, 1.6, 'III')
#ax3.loglog([Pth1]*2, [min(wlist/g),max(wlist/g)],'k:')
#ax3.loglog([Pth2]*2, [min(wlist/g),max(wlist/g)],'k:')
ax1.grid()
ax2.grid()  #we want a grid
ax4.grid()
##ax1
ax1.loglog(pump_list, nP_list_ss,'bo',mfc='none')   #cavity population
#ax1.xlabel('Pump [THz]')
#ax1.ylabel('Occupation probability of cavity') 
ax1.set(ylabel=r'$n_a$') 
TITLE='{}-emitter. $g$={:.3f} THz. $\kappa$={:.3f} THz. $\gamma_A$={:.3f} THz. $\gamma_D$={:.3f} THz'.format(N_em,g,kappa,gamma,gammaD)
ax1.set(title=TITLE)
ax1.set(xlim=[min(pump_list),max(pump_list)])
ax1.axes.xaxis.set_ticklabels([])

##ax2
ax2.loglog(pump_list, [[1]]*len(pump_list),'k--')          #1-line
ax2.loglog(pump_list, [[2]]*len(pump_list),'k--')          #2-line
ax2.loglog(pump_list, g2_list_ss,'bo',mfc='none')           #plot g2
ax2.set(ylabel=r'$g^{(2)}(0)$')
ax2.set(xlim=[min(pump_list),max(pump_list)])
ax2.set_yticks(np.arange(1,10))
ax2.set_yticklabels(np.arange(1,10)) #virker ikke lige nu
#ax2.axes.yaxis.set_ticklabels([1,1.5,2])
ax2.set(ylim=[min([-0.05+min(g2_list_ss),0.95]),max([max(g2_list_ss),2])+0.3])
print('g2 max',max(g2_list_ss))
##ax3
ax4.set(xlabel=r'$P$ [THz]') 
ax3.set(ylabel=r'$\omega-\omega_{eg}$')
#ax3.set(ylim=[-1.5*g,1.5*g])
#ax3.set_xticks(ticks=np.arange(-3,4,1),labels=10^(np.arange(-3,4,1)))
#g=np.sqrt(N_em)*g
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(2)-1)),(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(3)-np.sqrt(2))),(spectrum_matrix).transpose(), 20, cmap='RdGy')
ax3.contourf(np.log10(pump_list), wlist,(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss)[np.log10(pg_list_ss)<0], wlist,(spectrum_matrix[np.log10(pg_list_ss)<0]).transpose(), 20, cmap='RdGy')
#plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=ax3)
ax3.axes.xaxis.set_ticklabels([])

##ax4
#analytiske linjebredder
gammaR=4*g**2/(pump_list+gamma+gammaD+kappa)
Reg1_linewidth=(3*pump_list+gamma+gammaD+kappa)/2
Reg3_linewidth=kappa/(nP_list_ss+1)
neList=neMatrix[:,0]
print(np.shape(neMatrix))
ST_linewidth=gammaR/2*np.divide(neList,nP_list_ss)*N_em
print(neMatrix[0,:])

##plot linewidths
ST_unmodified=ST_linewidth*2
ax4.loglog(pump_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
ax4.loglog(pump_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
ax4.loglog(pump_list, ST_linewidth/(2*np.pi),'g:')    #plot ST (regime 2) linewidth
ax4.loglog(pump_list, kappa*np.ones((len(pump_list)))/(2*np.pi),'b')  #simple approx
ax4.loglog(pump_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax4.loglog(pump_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth

##plot thresholds
ax4.loglog([Pth1]*2, [ax4.axes.get_ylim()[0],ax4.axes.get_ylim()[1]],'k:')
ax4.loglog([Pth2]*2, [ax4.axes.get_ylim()[0],ax4.axes.get_ylim()[1]],'k:')
#ax4.loglog([Pth_Qn]*2, [ax4.axes.get_ylim()[0],ax4.axes.get_ylim()[1]],'k:')
#ax4.loglog([Pth_cl]*2/g, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'g:')
#ax4.loglog([Pth_qn]*2/g, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'b:')

##augmentation of plot
ax4.set(ylabel=r'$ \mathrm{FWHM}$')
ax4.set(xlim=[min(pump_list),max(pump_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax4.set(ylim=[min([min(linewidth)/2/(2*np.pi)]),max([max(linewidth)*1.3/(2*np.pi)])])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
ax4.axes.xaxis.set_ticklabels([])
#ax4.axes.yaxis.set_ticks([0.1,0.2,0.3,0.4])
#ax4.axes.yaxis.set_ticklabels([0.1,0.2,None,0.4])

plt.subplots_adjust(left=0.13,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
plt.savefig('./plots/na_g2_spetrum_linewidth_{}-emitter_pump-sweep_g={}_kap={}_gam={}_gam2={}.pdf'.format(N_em,g,kappa,gamma,gamma2))
#np.savez('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2),linewidth=linewidth)
#plt.show()
Dw=wlist[1]-wlist[0]
print()
print(max(wlist))
print(len(linewidth))
print(linewidth/Dw)