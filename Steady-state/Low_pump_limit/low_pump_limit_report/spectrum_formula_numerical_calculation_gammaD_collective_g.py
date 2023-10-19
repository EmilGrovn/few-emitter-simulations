#import packages
import numpy as np
import matplotlib.pyplot as plt
#consider an unpumped 1-emitter with parameters
N_em=2
g=2*np.sqrt(N_em) #coupling constant [ps^(-1)]
gammaA=0.1
kappa=8.5
#kappa=0.1
gammaD_list=10**np.arange(-2,3+0.1,0.1)
#gammaD_list=10**np.arange(-5,3+0.1,0.1)
#gammaD_list=10**np.arange(-2,10+0.1,0.1)
#omega-list
dw=0.0001*g
omega_list=np.arange(-2*g,2*g+dw,dw)
omega_list=np.arange(-3*g,3*g+dw,dw)
#omega_list=np.arange(-4*g,4*g+dw,dw)
#quantities which we will save
linewidth_list=np.zeros(len(gammaD_list))    #linewidth
spectrum_list=[[]]*len(gammaD_list)          #spectrum
K1_list=np.zeros(len(gammaD_list))           #K1
K2_list=np.zeros(len(gammaD_list))           #K2
peakloc_list=np.zeros(len(gammaD_list))      #peak location
B21A12_list=np.zeros(len(gammaD_list))
B22A22_list=np.zeros(len(gammaD_list))
B21A11_list=np.zeros(len(gammaD_list))
B22A21_list=np.zeros(len(gammaD_list))
Koef_lamb1=np.zeros(len(gammaD_list))        #coeffiecient of 1st spectrum term
Koef_lamb2=np.zeros(len(gammaD_list))        #coeffiecient of 2nd spectrum term
lambda1_list=np.zeros(len(gammaD_list))
lambda2_list=np.zeros(len(gammaD_list))
for ele,gammaD in enumerate(gammaD_list):
    print('kappa',kappa)
    #update matrix elements
    gammaCE=(kappa+gammaA+gammaD)/2
    gammaGE=(gammaA+gammaD)/2
    #the Liouvillian matrix in a basis of {EE,CC,EC,CE,GG,GE,GC,EG,CG}
    L=[[-gammaA ,      0,    1j*g,   -1j*g,      0,       0,       0,       0,       0],
    [0       , -kappa,   -1j*g,    1j*g,      0,       0,       0,       0,       0],
    [1j*g    ,  -1j*g,-gammaCE,       0,      0,       0,       0,       0,       0],
    [-1j*g   ,   1j*g,       0,-gammaCE,      0,       0,       0,       0,       0],
    [gammaA  ,  kappa,       0,       0,      0,       0,       0,       0,       0],
    [0       ,      0,       0,       0,      0,-gammaGE,    1j*g,       0,       0],
    [0       ,      0,       0,       0,      0,    1j*g,-kappa/2,       0,       0],
    [0       ,      0,       0,       0,      0,       0,       0,-gammaGE,   -1j*g],
    [0       ,      0,       0,       0,      0,       0,       0,   -1j*g,-kappa/2]]
    #calculate eigenvalues and -vectors
    Lambda,V=np.linalg.eig(L)
    #lambdaMatrix=np.diag(Lambda)    #diagonal matrix of eigenvalues
    B=V                             #eigenvector matrix
    A=np.linalg.inv(B)              #inverse eigenvector matrix
    
    #we are ready to calculate the spectrum using our formula
    #calculate the 4 infamous coefficients
    B21A11=B[6,5]*A[5,5]
    B22A21=B[6,6]*A[6,5]
    B21A11_list[ele]=B21A11
    B22A21_list[ele]=B22A21
    #check with our analytical calculations
    a0=-kappa/4
    a1=-(gammaA+gammaD)/4
    a2=a0-a1
    b0=np.sqrt(a2**2-g**2)
    print('coefficients', B21A11,B22A21,1j*g/(2*b0)) #check -B21A11=B22A21=pm ig/2b0 (fortegn afhænger af, hvilket fortegn vi vælger på egenvektoren)
    B21A12=B[6,5]*A[5,6]
    B22A22=B[6,6]*A[6,6]
    B21A12_list[ele]=B21A12
    B22A22_list[ele]=B22A22
    # #check
    # print('coefficient',B21A12,a2/(2*b0)+1/2) #check
    # print('coefficient',B22A22,-a2/(2*b0)+1/2) #check
    # Ban=np.array([[(-a2+b0)/np.sqrt((-a2+b0)**2+g**2),(-a2-b0)/np.sqrt((-a2-b0)**2+g**2)],
    #                             [1j*g/np.sqrt((-a2+b0)**2+g**2),1j*g/np.sqrt((-a2-b0)**2+g**2)]])
    # print('numerical eigenvectorMatrix',B[5:7,5:7])
    # print('analytical eigenvectorMatrix',Ban) #same apart from the second vector being scaled by i
    #eigenvalues
    lambda1=Lambda[5]
    lambda2=Lambda[6]
    lambda1_list[ele]=lambda1
    lambda2_list[ele]=lambda2

    #calculate K1 and K2
    #print(np.round(B,2))
    Alam=np.multiply(np.transpose(A[1:5,0]),1/Lambda[1:5])  #sum i=2 to 5 (i.e., 1 to 4 in python): 
    K1=np.sum(np.multiply(B[3,1:5],Alam))                   #row 4 of B, i=2 to 5
    K2=np.sum(np.multiply(B[1,1:5],Alam))                   #row 2 of B, i=2 to 5
    K1_list[ele]=K1
    K2_list[ele]=K2
    # print('K1',K1)
    # print('K2',K2)

    #S=kappa*((K2*B21A11+K1*B21A12)/(lambda1-1j*omega_list)+(K2*B22A21+K1*B22A22)/(lambda2-1j*omega_list))
    S=kappa*((K1*B21A11+K2*B21A12)/(lambda1-1j*omega_list)+(K1*B22A21+K2*B22A22)/(lambda2-1j*omega_list))
    # print('hey',K1*B21A11,K2*B21A12,K1*B22A21,K2*B22A22)
    # print('hey',K2*B21A11,K1*B21A12,K2*B22A21,K1*B22A22)
    # print('hey_brug',K1*B21A11+K2*B21A12,K1*B22A21+K2*B22A22)
    Koef_lamb1[ele]=K1*B21A11+K2*B21A12
    Koef_lamb2[ele]=K1*B22A21+K2*B22A22
    #print(S)#er komplex
    Spec=np.real(S)     #It turns out that we need the real part to agree with QuTiP. I do not know the physical argument
    Spec_norm=np.absolute(Spec)/np.max(np.absolute(Spec))   #normalize spectrum
    spectrum_list[ele]=Spec_norm                    #save normalized spectrum
    
    ### calculate linewidths
    spec_i=Spec_norm
    #use symmetry
    wlistPos=omega_list[omega_list>=0]
    specPos=spec_i[omega_list>=0]
    #identify maximum
    specMax=max(specPos)
    indexVec=np.arange(len(specPos))
    maxIndex=indexVec[specMax==specPos]
    wMax=wlistPos[maxIndex]
    #initialize while-loop
    #check for increasing omega
    Dw=omega_list[1]-omega_list[0]
    wInc=wMax
    specInc=specMax
    currentIndex=maxIndex
    wInc2=wlistPos[currentIndex+1]
    specInc2=specPos[currentIndex+1]
    boundary=max(wlistPos)
    while specInc2 >=0.5 and specInc2<specInc and wInc2<boundary:
        currentIndex+=1
        wInc=wInc2
        wInc2=wlistPos[currentIndex+1]
        specInc=specInc2
        specInc2=specPos[currentIndex+1]
    #check for increasing omega
    wDec=wMax
    specDec=specMax
    currentIndex=maxIndex
    wDec2=wDec-Dw
    if wDec2>=0:
        specDec2=specPos[currentIndex-1]
        while specDec2 >=0.5 and specDec2<specDec and wDec2>0:
            currentIndex-=1
            wDec=wDec2
            wDec2=wlistPos[currentIndex-1]
            specDec=specDec2
            specDec2=specPos[currentIndex-1]
    linewidth_list[ele]=2*max([abs(wInc-wMax),abs(wMax-wDec)])


#save linewidth and spectrum
spectrum_list=np.array(spectrum_list)
np.savez('./data/spectrum_formula_linewidth_g={}_kap={}_gamA={}.npz'.format(g,kappa,gammaA),gammaD_list=gammaD_list,linewidth_list=linewidth_list)
print(linewidth_list-kappa)
#prediction
a0=-kappa/4
a1=-(gammaD_list+gammaA)/4
a2=a0-a1
b0=np.sqrt(a2**2-g**2)
FWHM_linewidth=np.zeros((len(gammaD_list)))
FWHM_linewidth_approx=np.zeros((len(gammaD_list)))
for i in range(len(gammaD_list)):
    a2i=a2[i]
    a3i=a0+a1[i]
    b0i=b0[i]
    if abs(a2i)<g:
        FWHM_linewidth[i]=(kappa+gammaD_list[i]+gammaA)/2
    elif a2i>g:
        #gammaA and gammaD dominates
        FWHM_linewidth[i]=2*np.sqrt(-a3i**2-b0i**2+np.sqrt(2*a3i**4+2*b0i**4))
        FWHM_linewidth_approx[i]=kappa+4*g**2/(gammaD_list[i]+gammaA-kappa)
    elif -a2i>g:
        #kappa dominates
        FWHM_linewidth[i]=2*np.sqrt(-a3i**2-b0i**2+np.sqrt(2*a3i**4+2*b0i**4))
        FWHM_linewidth_approx[i]=gammaD_list[i]+gammaA+4*g**2/(kappa-gammaD_list[i]-gammaA)   

#plot quantities of interest
#make plots pretty
from matplotlib import rcParams
rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 16})
rcParams['axes.titlepad'] = 20

##linewidth
fig, ax1 = plt.subplots(1,figsize=(11,7))
ax1.loglog(gammaD_list,linewidth_list/(2*np.pi),'bo')
ax1.loglog(gammaD_list,(kappa)*np.ones(len(gammaD_list))/(2*np.pi),'b--')
ax1.loglog(gammaD_list,(gammaA+gammaD_list)*np.ones(len(gammaD_list))/(2*np.pi),'m--')
ax1.loglog(gammaD_list,2*np.absolute(np.real(lambda1_list))/(2*np.pi),'g')
ax1.loglog(gammaD_list,2*np.absolute(np.real(lambda2_list))/(2*np.pi),'r')
ax1.loglog(gammaD_list,FWHM_linewidth/(2*np.pi),'m.')
ax1.grid()
ax1.legend(['FWHM','$\kappa$','$\gamma_A+\gamma_D$','$2|\Re(\lambda_1)|$','$2|\Re(\lambda_2)|$','Analytical'])
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\kappa_D$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$.'.format(1,g,kappa,gammaA)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\gamma_D\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
##spectrum
fig, ax1 = plt.subplots(1,figsize=(11,7))
ax1.set(xlabel=r'$\mathrm{log}(\gamma_D)\, [\mathrm{ps^{-1}}]$') 
ax1.set(ylabel=r'$(\omega-\omega_{eg})\,\mathrm{[ps^{-1}]}$')
#ax3.set(ylim=[-1.5*g,1.5*g])
#ax3.set_xticks(ticks=np.arange(-3,4,1),labels=10^(np.arange(-3,4,1)))
ax1.contourf(np.log10(gammaD_list)[:50], omega_list,(spectrum_list).transpose()[:,:50], 20, cmap='RdGy')
plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=ax1)
plt.show()
# 
# fig, ax1 = plt.subplots(1,figsize=(11,7))
# Cvec1=(np.absolute(Koef_lamb2)-np.absolute(Koef_lamb1))/np.absolute(Koef_lamb1)
# Cvec2=(np.absolute(Koef_lamb1)-np.absolute(Koef_lamb2))/np.absolute(Koef_lamb2)
# bum=-3
# ax1.semilogx(gammaD_list[gammaD_list>bum],Cvec1,'m.')
# ax1.semilogx(gammaD_list[gammaD_list>bum],Cvec2,'b.')
# ax1.semilogx(gammaD_list[gammaD_list>bum],Koef_lamb1/np.absolute(Koef_lamb1),'m--')
# ax1.semilogx(gammaD_list[gammaD_list>bum],Koef_lamb2/np.absolute(Koef_lamb2),'b--')
# ax1.semilogx(gammaD_list[gammaD_list>bum],np.zeros(len(gammaD_list[gammaD_list>bum])),'k--')
# TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$.'.format(1,g,gammaD,gammaA)
# ax1.set(title=TITLE)
# ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$')#,ylabel=r'$(|C2|-|C1|)/|C1|$')
# ax1.legend(['$(|C2|-|C1|)/|C1|$','$(|C1|-|C2|)/|C2|$','sign(C1)','sign(C2)','0'])
# ax1.set(ylim=[-5,5])
# #plot regimes
# Q12=-4*g+kappa-gammaA
# Q23=4*g+kappa-gammaA
# Q2len=(Q23-Q12)
# print(Q12)
# Min=ax1.get_ylim()[0]
# Max=ax1.get_ylim()[1]
# plt.semilogx([Q12,Q12],[Min,Max],'k--')
# plt.text(Q12/10, Max/1.5, '$a_2>g$')
# plt.semilogx([Q23,Q23],[Min,Max],'k--')
# plt.text(Q12+Q2len/10, Max/1.5, '$|a_2|<g$')
# plt.text(Q23*10, Max/1.5, '$-a_2>g$')

##sign and relative size of coefficients
fig, ax1 = plt.subplots(1,figsize=(11,7))
Ko1=np.multiply(Koef_lamb1,np.real(lambda1_list))
Ko2=np.multiply(Koef_lamb2,np.real(lambda2_list))
Cvec1=(np.absolute(Ko2)-np.absolute(Ko1))/np.absolute(Ko1)
Cvec2=(np.absolute(Ko1)-np.absolute(Ko2))/np.absolute(Ko2)
print(Ko1)
print(Cvec1)
bum=-4
ax1.semilogx(gammaD_list[gammaD_list>bum],Cvec1,'m.')
ax1.semilogx(gammaD_list[gammaD_list>bum],Cvec2,'b.')
ax1.semilogx(gammaD_list[gammaD_list>bum],Ko1/np.absolute(Ko1),'m--')
ax1.semilogx(gammaD_list[gammaD_list>bum],Ko2/np.absolute(Ko2),'b--')
ax1.semilogx(gammaD_list[gammaD_list>bum],np.zeros(len(gammaD_list[gammaD_list>bum])),'k--')
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\kappa$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$.'.format(1,g,kappa,gammaA)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$')#,ylabel=r'$(|C2|-|C1|)/|C1|$')
ax1.legend(['$(|C2\Re{\lambda_2}|-|C1\Re{\lambda_1}|)/|C1\Re{\lambda_1}|$','$(|C1\Re{\lambda_1}|-|C2\Re{\lambda_2}|)/|C2\Re{\lambda_2}|$','sign(C1$\Re{\lambda_1}$)','sign(C2$\Re{\lambda_2}$)','0'])
ax1.set(ylim=[-2,2])
#plot regimes
Q12=-4*g+(-gammaA+kappa)
Q23=4*g+(-gammaA+kappa)
Q2len=(Q23-Q12)
print(Q12)
Min=ax1.get_ylim()[0]
Max=ax1.get_ylim()[1]
plt.semilogx([Q12,Q12],[Min,Max],'k--')
plt.text(Q12/10, Max/1.5, '$-a_2>g$')
plt.semilogx([Q23,Q23],[Min,Max],'k--')
plt.text(Q12+Q2len/10, Max/1.5, '$|a_2|<g$')
plt.text(Q23*10, Max/1.5, '$a_2>g$')
plt.show()
#
fig, ax1 = plt.subplots(1,figsize=(11,7))
print(Ko1)
print(Cvec1)
bum=-3
ax1.semilogx(gammaD_list[gammaD_list>bum],Ko1,'m.')
ax1.semilogx(gammaD_list[gammaD_list>bum],Ko2,'b.')
ax1.semilogx(gammaD_list[gammaD_list>bum],np.absolute(Ko2)-np.absolute(Ko1),'r.')
ax1.semilogx(gammaD_list[gammaD_list>bum],np.zeros(len(Ko1)),'k--')
TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\kappa$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$.'.format(1,g,kappa,gammaA)
ax1.set(title=TITLE)
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$')#,ylabel=r'$(|C2|-|C1|)/|C1|$')
ax1.legend(['$C2\Re{\lambda_2}$','$C1\Re{\lambda_1}$','$|C2\Re{\lambda_2}|-|C1\Re{\lambda_1}|$','0'])
#ax1.set(ylim=[-5,5])
#plot regimes
Q12=-4*g+(-gammaA+kappa)
Q23=4*g+(-gammaA+kappa)
Q2len=(Q23-Q12)
print(Q12)
Min=ax1.get_ylim()[0]
Max=ax1.get_ylim()[1]
plt.semilogx([Q12,Q12],[Min,Max],'k--')
plt.text(Q12/10, Max/1.5, '$-a_2>g$')
plt.semilogx([Q23,Q23],[Min,Max],'k--')
plt.text(Q12+Q2len/10, Max/1.5, '$|a_2|<g$')
plt.text(Q23*10, Max/1.5, '$a_2>g$')
plt.show()