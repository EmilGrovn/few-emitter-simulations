#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from pump_sweep_variable_NH import pump_sweep_variable_NH
###
#pretty plot
from matplotlib import rcParams
rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 12})
rcParams['axes.titlepad'] = 20

###parameters
#system
g=0.1 #coupling constant [THz]
kappa=0.02 #decay rate [THz] for coupling from cavity to environment
pump=0
gamma=0.012
gamma2=2/np.sqrt(2)
pump_logmin=-3
pump_logmax=3
#numerical
N_em_list=[2]
g2_matrix=[[]]*len(N_em_list)*2
nP_matrix=[[]]*len(N_em_list)*2
linewidth_matrix=[[]]*len(N_em_list)*2
tol=10
legend=[None]*len(N_em_list)*2
###plot

for i,N_em in enumerate(N_em_list):
    out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(N_em,g,kappa,gamma,gamma2,tol))
    g2_list_ss=out_ss['g2']
    nP_list_ss=out_ss['nP']
    pg_list_ss=out_ss['pump_over_g_save']
    g2_matrix[i]=g2_list_ss
    nP_matrix[i]=nP_list_ss
    out_ss=np.load('./data/{}-emitter_pump-sweep_ss_NH_intelligent_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(1,g,kappa/N_em,gamma,gamma2,tol))
    g2_list_ss=out_ss['g2']
    nP_list_ss=out_ss['nP']
    g2_matrix[i+len(N_em_list)]=g2_list_ss
    nP_matrix[i+len(N_em_list)]=nP_list_ss
    #linewidth
    out_ss=np.load('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(N_em,g,kappa,gamma,gamma2,tol))
    linewidth_matrix[i]=out_ss['linewidth']
    out_ss=np.load('./data/{}-emitter_pump-sweep_linewidth_g={}_kap={}_gam={}_gam2={}_tol={}.npz'.format(1,g,kappa/N_em,gamma,gamma2,tol))
    linewidth_matrix[i+len(N_em_list)]=out_ss['linewidth']
    legend[i]='{}-emitter'.format(N_em)
    legend[i+len(N_em_list)]=r'{}-emitter with $\kappa/N$ for $N=${}'.format(1,N_em)
style=['b*','g*','r*']
style1=['b-','g-','g--','g--','k:','k:']
style2=['b--','g--','r--']
style3=['b:','g:','r:']
fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(10,7))
TITLE='$\gamma_D$={:.2f} THz'.format(gamma2*np.sqrt(2))
ax1.set(title=TITLE)
#actual plot
for i in range(len(N_em_list)*2):
    ax1.loglog(pg_list_ss, nP_matrix[i],style1[i])
for i in range(len(N_em_list)*2):
    ax2.loglog(pg_list_ss, g2_matrix[i],style[i])
for i in range(len(N_em_list)*2):
    ax3.loglog(pg_list_ss, linewidth_matrix[i]*(i+1),style[i])
#specifications
ax1.grid()
ax2.grid()
#ax1.xlabel('Pump/g [THz]')
#ax1.ylabel('Occupation probability of cavity') 
ax1.set(ylabel=r'$n_a$') 
TITLE='{}-emitter'.format(N_em)
#ax1.set(title=TITLE)
ax1.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
ax1.axes.xaxis.set_ticklabels([])
#ax2
ax2.loglog(pg_list_ss, [[1]]*len(pg_list_ss),'k--')
ax2.loglog(pg_list_ss, [[2]]*len(pg_list_ss),'k--')
ax2.set(xlabel=r'$P/g$',ylabel=r'$g^{(2)}(0)$')
ax2.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
#ax2.set_yticks(ticks=[1,2,6,9],labels={1,2,6,9})
#ax2.set(ylim=[-0.05+min(g2_matrix[2]),max([max(g2_matrix[0]),2])+0.3])
ax3.set(xlabel=r'$P/g$',ylabel='FWHM')
ax3.set(xlim=[min(pg_list_ss),max(pg_list_ss)])
ax1.legend(legend,loc='best')
plt.tight_layout()
plt.savefig('./plots/N-emitter_plot_eff_single_em_g={}_kap={}_gam={}_gam2={}_tol={}.pdf'.format(g,kappa,gamma,gamma2,tol))
plt.show()