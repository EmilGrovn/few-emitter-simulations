#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from pump_sweep_fixed_NH import pump_sweep_fixed_NH
from pump_sweep_variable_NH import pump_sweep_variable_NH 
from pump_sweep_spec import pump_sweep_spec
###
#pretty plot
from matplotlib import rcParams
rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 16})
rcParams['axes.titlepad'] = 20

###PARAMETERS - DETERMINE SYSTEM
#PHYSICAL
N_em=1
g=0.1 #coupling constant [THz]
#g=2
#g=6
kappa=0.02 #decay rate [THz] for coupling from cavity to environment
#pump=0
gamma=0.012
gamma2=0.1/np.sqrt(2)
#gamma2=1/sqrt(2)
#pump_logmin=-3
#pump_logmax=3
#NUMERICAL
N_Hilbert=10
tol=10

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(N_em,g,kappa,gamma,gamma2,tol))
wlist=out_ss['w']
spectrum_matrix=out_ss['spec']
pg_list_ss=out_ss['pump_over_g_save']
out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(N_em,g,kappa,gamma,gamma2,tol))
g2_list_ss=out_ss['g2']
nP_list_ss=out_ss['nP']
ne_list=out_ss['neList']
cav_conv_list=out_ss['cav_conv']
N_H_list=out_ss['N_H_list']
#print(nP_list_ss[max(cav_conv_list)==cav_conv_list])
#print(N_H_list[max(cav_conv_list)==cav_conv_list])
#print(max(nP_list_ss))

## identifying different regimes using g2
delta=0.05
Pth1=min(pg_list_ss[g2_list_ss>1-delta])
Pth2=max(pg_list_ss[g2_list_ss<1+delta])
Pth_Qn=min(pg_list_ss[g2_list_ss>2-delta])
## calculate linewidth
print(np.shape(spectrum_matrix))
print(np.max(spectrum_matrix,axis=1)) #every spectrum has been normalized to have a unit peak
print(wlist[spectrum_matrix[0,:]>0.5]/g)
linewidth=np.zeros(len(pg_list_ss))
for i in np.arange(0,len(pg_list_ss)):
    spec_i=spectrum_matrix[i,:]
    wlistAboveHM=wlist[spec_i>0.5]
    specAboveHM=spec_i[spec_i>0.5]
    linewidth[i]=np.max(wlistAboveHM)-np.min(wlistAboveHM)

    if len(specAboveHM[wlistAboveHM==0])==0:
        print('Oh')
        #then there are two seperate peaks for negative and positive omega
        linewidth[i]=np.max(wlistAboveHM[wlistAboveHM>0])-np.min(wlistAboveHM[wlistAboveHM>0])
        print(linewidth[i])
    elif specAboveHM[wlistAboveHM==0]<1:
        print('Oh',i,pg_list_ss[i])
        #then there are two seperate peaks for negative and positive omega
        wlistAboveHMPosOmega=wlistAboveHM[wlistAboveHM>0]
        specAboveHMPosOmega=specAboveHM[wlistAboveHM>0]

        linewidth[i]=2*(np.max(wlistAboveHMPosOmega)-wlistAboveHMPosOmega[np.max(specAboveHMPosOmega)==specAboveHMPosOmega])
        print(linewidth[i])
print(linewidth)
## FIGURES
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4)
#plot g2=1-thresholds
#ax1.loglog([Pth1]*2, [min(nP_list_ss),max(nP_list_ss)],'k:')
#ax1.loglog([Pth2]*2, [min(nP_list_ss),max(nP_list_ss)],'k:')
ax2.loglog([Pth1]*2, [0.5,3],'k:')
ax2.loglog([Pth2]*2, [0.5,3],'k:')  #plot regime boundaries only on middle plot
#ax3.loglog([Pth1]*2, [min(wlist/g),max(wlist/g)],'k:')
#ax3.loglog([Pth2]*2, [min(wlist/g),max(wlist/g)],'k:')
ax4.loglog([Pth1]*2, [min(wlist/g),max(wlist/g)],'k:')
ax4.loglog([Pth2]*2, [min(wlist/g),max(wlist/g)],'k:')
ax1.grid()
ax2.grid()  #we want a grid
ax4.grid()
#AX1
ax1.loglog(pg_list_ss, nP_list_ss[:len(pg_list_ss)],'b*')    #cavity population
#ax1.xlabel('Pump/g [THz]')
#ax1.ylabel('Occupation probability of cavity')
#augmentation
ax1.set(ylabel=r'$n_a$') 
TITLE='{}-emitter. $\gamma_D$={:.2f} THz'.format(N_em,gamma2*np.sqrt(2))
ax1.set(title=TITLE)
ax1.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
ax1.axes.xaxis.set_ticklabels([])
#AX2
ax2.loglog(pg_list_ss, [[1]]*len(pg_list_ss),'k--')          #1-line
ax2.loglog(pg_list_ss, [[2]]*len(pg_list_ss),'k--')          #2-line
ax2.loglog(pg_list_ss, g2_list_ss[:len(pg_list_ss)],'b*')    #plot g2
#augmentation
ax2.set(ylabel=r'$g^{(2)}(0)$')
ax2.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
ax2.set_yticks(np.arange(1,10))
ax2.set_yticklabels(np.arange(1,10)) #virker ikke lige nu
#ax2.axes.yaxis.set_ticklabels([1,1.5,2])
ax2.set(ylim=[-0.05+min(g2_list_ss),max([max(g2_list_ss),2])+0.3])
ax2.text(10**(np.log10(Pth2)*0.2+np.log10(Pth1)*0.8),1.5, '$g^{(2)}=1\pm$'+'{}'.format(delta))
#AX3
#augmentation
ax4.set(xlabel=r'$P/g$',ylabel=r'$\omega/g$') 
ax3.set(ylim=[-1.5,1.5])
#ax3.set_xticks(ticks=np.arange(-3,4,1),labels=10^(np.arange(-3,4,1)))
#gCol=np.sqrt(N_em)*g   #collective g
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(2)-1)),(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(3)-np.sqrt(2))),(spectrum_matrix).transpose(), 20, cmap='RdGy')
ax3.contourf(np.log10(pg_list_ss), wlist/g,(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss)[np.log10(pg_list_ss)<0], wlist,(spectrum_matrix[np.log10(pg_list_ss)<0]).transpose(), 20, cmap='RdGy')
#plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=ax3)
ax3.axes.xaxis.set_ticklabels([])
#AX4
#analytiske linjebredder
gammaR=4*g**2/(pg_list_ss*g+gamma+np.sqrt(2)*gamma2+kappa)
Reg1_linewidth=(3*pg_list_ss*g+gamma+np.sqrt(2)*gamma2+kappa)/2
Reg3_linewidth=kappa/(nP_list_ss+1)
ST_linewidth=gammaR/2*np.divide(ne_list,nP_list_ss)
print(ne_list)
print(nP_list_ss)
print(gammaR)
#thresholds - forstår ikke helt fysikken nok til at finde ud af, 
# hvordan jeg får de rigtige thresholds
# får samme Pcl som Matias, så det tyder jo på, at jeg har gættet rigtig
# men alligevel er Ppr helt ude at skide
ne_th=ne_list[Pth1==pg_list_ss][0]
print('ne_th',ne_th)
gammaR_th=gammaR[Pth1==pg_list_ss][0]
gammaC=gammaR_th*(2*ne_th-1)
print(gammaC)
xi=gammaR_th/(2*gammaC)
beta=gammaR_th/(gammaR_th+gamma)
print(beta)
#beta_eff=beta/(1+2*xi*(1-beta))
factor=(2/(1-1/(2*xi)))/2*gammaC/beta
print(factor)
print((1+2*xi))
print((1+2*xi*(1-beta)))
#print(beta_eff)
Pth_cl=factor*(1+2*xi)
Pth_pr=factor*(1+2*xi*(1-beta))
#Pth_pr=0.03
Pth_qn=(4*g**2-kappa**2-2*gamma*kappa-kappa*np.sqrt(2)*gamma2+np.sqrt((4*g**2)**2-24*g**2*kappa**2-32*g**2*kappa*gamma-8*g**2*kappa*np.sqrt(2)*gamma2+kappa**4+2*kappa*3*np.sqrt(2)*gamma2+kappa*2*(np.sqrt(2)*gamma2)**2))/(2*kappa)
print('Pth1',Pth1)
print('g',g)
print(Pth_cl)
print(Pth_pr)
print(Pth_qn)
#plot linewidths
ax4.loglog(pg_list_ss, linewidth[:len(pg_list_ss)]/g,'b*')    #plot linewidth
ax4.loglog(pg_list_ss, ST_linewidth[:len(pg_list_ss)]/g,'g')    #plot ST (regime 2) linewidth
ax4.loglog(pg_list_ss, Reg1_linewidth[:len(pg_list_ss)]/g,'r')  #regime 1 linewidth
ax4.loglog(pg_list_ss, Reg3_linewidth[:len(pg_list_ss)]/g,'k')  #regime 3 linewidth
#plot thresholds
ax4.text(10**(np.log10(Pth2)*0.2+np.log10(Pth1)*0.8),1.5, '$g^{(2)}=1\pm$'+'{}'.format(delta))
ax4.loglog([Pth1]*2, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'r:')
ax4.loglog([Pth2]*2, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'r:')
ax4.loglog([Pth_Qn]*2, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'k:')
#ax4.loglog([Pth_cl]*2/g, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'g:')
#ax4.loglog([Pth_qn/g]*2, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'k:')
#ax4.text((Pth1)*0.02,min(Reg1_linewidth), 'RegI-FWHM')
#ax4.text(10**(np.log10(Pth2)*0.2+np.log10(Pth1)*0.8),min(Reg1_linewidth)/3, 'ST-FWHM')
#ax4.text(Pth_qn*1.1/g,min(Reg1_linewidth)/3, 'Qn-FWHM')
#augmentation of plot
ax4.set(ylabel=r'$ \mathrm{FWHM}/g$')
ax4.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax4.set(ylim=[min([min(linewidth/g/3)]),max([max(linewidth/g)+0.05])])
#ax4.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10))
ax4.axes.xaxis.set_ticklabels([])
#ax4.axes.yaxis.set_ticks([0.1,0.2,0.3,0.4])
#ax4.axes.yaxis.set_ticklabels([0.1,0.2,None,0.4])

plt.subplots_adjust(left=0.12,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
plt.savefig('./plots/superplot_{}-emitter_pump-sweep_g={}_kap={}_gam={}_gam2={}_tol={}.pdf'.format(N_em,g,kappa,gamma,gamma2,tol))
#plt.show()