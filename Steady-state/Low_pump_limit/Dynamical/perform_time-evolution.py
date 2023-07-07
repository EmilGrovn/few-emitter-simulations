#import
from cmath import sqrt
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
#import implementation of solving single-emitter
from time_evolution import time_evolution
#from time_evolution_IC import time_evolution_IC
###
#pretty plot
from matplotlib import rcParams
rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 16})
rcParams['axes.titlepad'] = 20

###parameters
#system
N_em=1
g=2 #coupling constant [THz]
kappa=100 #decay rate [THz] for coupling from cavity to environment
pump=0
gammaA=0
gammaD=0
gamma=gammaA
gamma2=gammaD/np.sqrt(2)
#numerical
N_Hilbert=N_em+1
Nt=2000
Tw=60 #[ps]

print(g,kappa,gammaA,gammaD)
#numerical
#N_Hilbert=N_em+1
#Nt=600
#Tw=g*10**4 #[ps] scales with g
#Tw=100

###perform pump-sweep
time_evolution(N_em,g,kappa,pump,gamma,gamma2,Tw,Nt,N_Hilbert)

out=np.load('.\data\{}-emitter_time_evolution_pump={}_Tw={}_Nt={}_NH={}.npz'.format(N_em,pump,Tw,Nt,N_Hilbert))
c_occ_list=out['c_occ']
e_occ_list=out['e_occM']
tlist=out['tlist']
print('Last emitter occupation',e_occ_list[:,-1])
#
fig = plt.figure(1,figsize=(4,4))
styles=['b','b--','b:']
for i in range(N_em):
    plt.plot(tlist, e_occ_list[i].flatten(),styles[i],label="Emitter "+'{}'.format(i+1))
plt.plot(tlist, c_occ_list,'r',label="Cavity")
#plt.semilogy(tlist, e_occ_list.flatten(),'b',label="Emitter") 
#plt.semilogy(tlist, c_occ_list,'r',label="Cavity") 
plt.xlabel('Time [ps]') 
plt.ylabel('Occupation probability') 
#plt.legend(bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=1)
plt.legend()
plt.xlim(0,Tw)
#plt.xlim(0,0.25)
plt.ylim(0,1)
Title=r''
if g<abs(kappa-gamma)/4:
    Title+='Weak regime'
elif g==abs(kappa-gamma)/4:
    Title+='Weak-int boarder'
elif g<np.sqrt((kappa**2+gamma**2)/8):
    Title+='Intermediate regime'
elif g==np.sqrt((kappa**2+gamma**2)/8):
    Title+='Int-strong boarder'
else:
    Title+='Strong regime'
if gamma==0 and kappa > 0:
    Title+='. g/$\kappa$={:.2f}'.format(g/kappa)
else:
    Title+='. g={:.0f} GHz'.format(g*10**3)
plt.title(Title)
plt.tight_layout()
plt.grid()
#plt.savefig('./rep_plots/population_dynamics_g={}_kappa={}.pdf'.format(g,kappa))
#plt.savefig('./plots/population_dynamics_figaGod_g={}_kappa={}.pdf'.format(g,kappa))
plt.savefig('./plots/population_dynamics_2_Nem={}_g={}_kap={}_gamA={}_gamD={}.pdf'.format(N_em,g,kappa,gammaA,gammaD))
#plt.savefig('./plots/log_population_dynamics_g={}_kappa={}.pdf'.format(g,kappa))
#plt.show()