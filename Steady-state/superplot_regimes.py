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
gamma2=0.1/np.sqrt(2)
gammaD=np.sqrt(2)*gamma2
correctLW=True
#pump_logmin=-3
#pump_logmax=3
#numerical
#N_Hilbert=10

###IMPORT DATA
out_ss=np.load('./data/{}-emitter_pump-sweep_spectrum_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
wlist=out_ss['w']
spectrum_matrix=out_ss['spec']
pg_list_ss=out_ss['pump_over_g_save']
out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2))
g2_list_ss=out_ss['g2']
nP_list_ss=out_ss['nP']
neMatrix=out_ss['neMatrix']
cav_conv_list=out_ss['cav_conv']
N_H_list=out_ss['N_H_list']
#deprecated print-statements
#print(nP_list_ss[max(cav_conv_list)==cav_conv_list])
#print(N_H_list[max(cav_conv_list)==cav_conv_list])
#print('ne-list', neMatrix)
#print(max(nP_list_ss))
#print(nP_list_ss)
#print(N_H_list)
#print(pg_list_ss)

## identifying different regimes using g2 ##
delta=0.05
try:
    Pth1=min(pg_list_ss[g2_list_ss>1-delta])
except ValueError:
    Pth1=max(pg_list_ss)
try:
    Pth2=min(pg_list_ss[g2_list_ss<1+delta])
except ValueError:
    Pth2=max(pg_list_ss)
try:
    Pth_Qn=min(pg_list_ss[g2_list_ss>2-delta])
except ValueError:
    Pth_Qn=max(pg_list_ss)


## calculate linewidth ##
linewidth=np.zeros(len(pg_list_ss))
if correctLW==False:
    for i in np.arange(0,len(pg_list_ss)):
        spec_i=spectrum_matrix[i,:]     #choose spectrum for given pump
        wlistAboveHM=wlist[spec_i>=0.5] #extract frequencies where spectrum is above half max
        specAboveHM=spec_i[spec_i>=0.5] #extract corresponding spectrum values

        linewidth[i]=np.max(wlistAboveHM)-np.min(wlistAboveHM)  #simple case with one pulse at w=0
        #check for multiple peaks
        DeltaW=wlist[1]-wlist[0]    #frequency spacing
        #if linewidth[i]/DeltaW>len(wlistAboveHM):
        Bool=wlistAboveHM==0    #check if w=0 is part of above half max spectrum
        print(i)
        print(linewidth[i])
        if all(Bool==False):
            print('Oh')
            #then there are two seperate peaks for negative and positive omega
            linewidth[i]=np.max(wlistAboveHM[wlistAboveHM>0])-np.min(wlistAboveHM[wlistAboveHM>0])
            print(linewidth[i])
        #does not account for the case with partly seperated peaks!
if correctLW==True:
    for i in np.arange(0,len(pg_list_ss)):
        spec_i=spectrum_matrix[i,:]     #choose spectrum for given pump
        wlistAboveHM=wlist[spec_i>=0.5] #extract frequencies where spectrum is above half max
        specAboveHM=spec_i[spec_i>=0.5] #extract corresponding spectrum values

        linewidth[i]=np.max(wlistAboveHM)-np.min(wlistAboveHM)  #simple case with one pulse at w=0

        if len(specAboveHM[wlistAboveHM==0])==0:
            print('Oh')
            #then there are two seperate peaks for negative and positive omega
            linewidth[i]=np.max(wlistAboveHM[wlistAboveHM>=0])-np.min(wlistAboveHM[wlistAboveHM>=0])
            #I am not sure this is entirely correct because there might still be a contribution from the other peak
        elif specAboveHM[wlistAboveHM==0]<1:
            print('Oh',i,pg_list_ss[i])
            #then there are two partly seperated peaks for negative and positive omega
            wlistAboveHMPosOmega=wlistAboveHM[wlistAboveHM>=0] #Positive frequencies
            specAboveHMPosOmega=specAboveHM[wlistAboveHM>=0]   #corresponding spectrum values
            linewidth[i]=2*(np.max(wlistAboveHMPosOmega)-wlistAboveHMPosOmega[np.max(specAboveHMPosOmega)==specAboveHMPosOmega])
##
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,figsize=(11,7))
#plot thresholds
#ax1.loglog([Pth1]*2, [min(nP_list_ss),max(nP_list_ss)],'k:')
#ax1.loglog([Pth2]*2, [min(nP_list_ss),max(nP_list_ss)],'k:')
ax2.loglog([Pth1]*2, [0.5,3],'k:')
ax2.loglog([Pth2]*2, [0.5,3],'k:')
ax2.loglog([Pth_Qn]*2, [0.5,3],'k:')  #plot regime boundaries only on middle plot
#ax3.loglog([Pth1]*2, [min(wlist/g),max(wlist/g)],'k:')
#ax3.loglog([Pth2]*2, [min(wlist/g),max(wlist/g)],'k:')
ax1.grid()
ax2.grid()  #we want a grid
ax4.grid()
pump_list=pg_list_ss*g
##ax1
ax1.loglog(pump_list, nP_list_ss,'bo',mfc='none')   #cavity population
#ax1.xlabel('Pump [THz]')
#ax1.ylabel('Occupation probability of cavity') 
ax1.set(ylabel=r'$n_a$') 
TITLE='{}-emitter. $\gamma_D$={:.2f} THz'.format(N_em,gammaD)
ax1.set(title=TITLE)
ax1.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
ax1.axes.xaxis.set_ticklabels([])

##ax2
ax2.loglog(pump_list, [[1]]*len(pg_list_ss),'k--')          #1-line
ax2.loglog(pump_list, [[2]]*len(pg_list_ss),'k--')          #2-line
ax2.loglog(pump_list, g2_list_ss,'bo',mfc='none')           #plot g2
ax2.set(ylabel=r'$g^{(2)}(0)$')
ax2.set(xlim=[min(pump_list),max(pump_list)])
ax2.set_yticks(np.arange(1,10))
ax2.set_yticklabels(np.arange(1,10)) #virker ikke lige nu
#ax2.axes.yaxis.set_ticklabels([1,1.5,2])
ax2.set(ylim=[min([-0.05+min(g2_list_ss),0.95]),max([max(g2_list_ss),2])+0.3])

##ax3
ax4.set(xlabel=r'$P [THz]$') 
ax3.set(ylabel=r'$(\omega-\omega_{eg})$')
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
#thresholds - forstår ikke helt fysikken nok til at finde ud af, 
# hvordan jeg får de rigtige thresholds
# får samme Pcl som Matias, så det tyder jo på, at jeg har gættet rigtig
# men alligevel er Ppr helt ude at skide
#neMatrix=neMatrix[0,:]
ne_th=neMatrix[Pth1==pg_list_ss][0][0]
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
print('classic',Pth_cl)
print(Pth_pr)
print(Pth_qn)

##plot linewidths
ST_unmodified=ST_linewidth*2
ax4.loglog(pump_list, linewidth/(2*np.pi),'bo',mfc='none')    #plot linewidth
ax4.loglog(pump_list, ST_unmodified/(2*np.pi),'g')    #plot ST (regime 2) linewidth
ax4.loglog(pump_list, kappa*np.ones((len(pg_list_ss))),'b')  #simple approx
ax4.loglog(pump_list, Reg1_linewidth/(2*np.pi),'r')  #regime 1 linewidth
ax4.loglog(pump_list, Reg3_linewidth/(2*np.pi),'k')  #regime 3 linewidth

##plot thresholds
ax4.loglog([Pth1]*2, [min([min(linewidth),0.001]),max([max(linewidth),10])],'k:')
ax4.loglog([Pth2]*2, [min([min(linewidth),0.001]),max([max(linewidth),10])],'k:')
ax4.loglog([Pth_Qn]*2, [min([min(linewidth),0.001]),max([max(linewidth),10])],'k:')
#ax4.loglog([Pth_cl]*2/g, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'g:')
#ax4.loglog([Pth_qn]*2/g, [min([min(linewidth/g),0.001]),max([max(linewidth/g),10])],'b:')

##augmentation of plot
ax4.set(ylabel=r'$ \mathrm{FWHM}$')
ax4.set(xlim=[min(pump_list),max(pump_list)])
#ax2.set_yticks(ticks=np.arange(1,10),labels=np.arange(1,10)) #virker ikke lige nu
ax4.set(ylim=[min([min(linewidth)/2]),max([max(linewidth)*1.3])])
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
plt.savefig('./plots/superplot_{}-emitter_pump-sweep_g={}_kap={}_gam={}_gam2={}.pdf'.format(N_em,g,kappa,gamma,gamma2))
#np.savez('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}.npz'.format(N_em,g,kappa,gamma,gamma2),linewidth=linewidth)
#plt.show()
Dw=wlist[1]-wlist[0]
print()
print(max(wlist))
print(len(linewidth))
print(linewidth/Dw)