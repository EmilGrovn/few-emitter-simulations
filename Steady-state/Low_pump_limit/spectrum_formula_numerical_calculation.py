import numpy as np
import matplotlib.pyplot as plt
#consider an unpumped 1-emitter with parameters
g=2 #coupling constant [ps^(-1)]
kappa=100 #decay rate [ps^(-1)] for coupling from cavity to environment
kappa_list=10**np.arange(-2,3+0.1,0.1)
linewidth_list=np.zeros(len(kappa_list))
dw=0.0001*g
omega_list=np.arange(-2*g,2*g+dw,dw)
omega_list=np.arange(-3*g,3*g+dw,dw)
spectrum_list=[[]]*len(kappa_list)
K1_list=np.zeros(len(kappa_list))
K2_list=np.zeros(len(kappa_list))
peakloc_list=np.zeros(len(kappa_list))
B21A12_list=np.zeros(len(kappa_list))
B22A22_list=np.zeros(len(kappa_list))
B21A11_list=np.zeros(len(kappa_list))
B22A21_list=np.zeros(len(kappa_list))
Koef_lamb1=np.zeros(len(kappa_list))
Koef_lamb2=np.zeros(len(kappa_list))
for ele,kappa in enumerate(kappa_list):
    print('kappa',kappa)
    gammaA=0
    gammaD=0
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
    Lambda,V=np.linalg.eig(L)
    lambdaMatrix=np.diag(Lambda)
    omega=0
    omegaImatrix=np.diag(np.ones(9)/1j*omega)
    B=V
    A=np.linalg.inv(B)
    #we are ready to calculate the spectrum using our formula

    #calculate the 4 infamous coefficients
    B21A11=B[6,5]*A[5,5]#+1
    B22A21=B[6,6]*A[6,5]
    B21A11_list[ele]=B21A11
    B22A21_list[ele]=B22A21
    a0=-kappa/4
    a1=-(gammaA+gammaD)/4
    a2=a0-a1
    b0=np.sqrt(a2**2-g**2)
    print('coefficients', B21A11,B22A21,1j*g/(2*b0)) #check -B21A11=B22A21=pm ig/2b0 (fortegn afhænger af, hvilket fortegn vi vælger på egenvektoren)
    B21A12=B[6,5]*A[5,6]
    B22A22=B[6,6]*A[6,6]
    print('coefficient',B21A12,a2/(2*b0)+1/2) #check
    print('coefficient',B22A22,-a2/(2*b0)+1/2) #check
    Ban=np.array([[(-a2+b0)/np.sqrt((-a2+b0)**2+g**2),(-a2-b0)/np.sqrt((-a2-b0)**2+g**2)],
                                [1j*g/np.sqrt((-a2+b0)**2+g**2),1j*g/np.sqrt((-a2-b0)**2+g**2)]])
    print('numerical eigenvectorMatrix',B[5:7,5:7])
    print('analytical eigenvectorMatrix',Ban) #same apart from the second vector being scaled by i
    #eigenvalues
    lambda1=Lambda[5]
    lambda2=Lambda[6]
    B21A12_list[ele]=B21A12
    B22A22_list[ele]=B22A22

    #calculate K1 and K2
    print(np.round(B,2))
    Alam=np.multiply(np.transpose(A[1:5,0]),1/Lambda[1:5])
    K1=np.sum(np.multiply(B[3,1:5],Alam))
    K2=np.sum(np.multiply(B[1,1:5],Alam))
    print('K1',K1)
    print('K2',K2)
    K1_list[ele]=K1
    K2_list[ele]=K2

    #S=kappa*((K2*B21A11+K1*B21A12)/(lambda1-1j*omega_list)+(K2*B22A21+K1*B22A22)/(lambda2-1j*omega_list))
    S=kappa*((K1*B21A11+K2*B21A12)/(lambda1-1j*omega_list)+(K1*B22A21+K2*B22A22)/(lambda2-1j*omega_list))
    print('hey',K1*B21A11,K2*B21A12,K1*B22A21,K2*B22A22)
    print('hey',K2*B21A11,K1*B21A12,K2*B22A21,K1*B22A22)
    print('hey_brug',K1*B21A11+K2*B21A12,K1*B22A21+K2*B22A22)
    Koef_lamb1[ele]=K1*B21A11+K2*B21A12
    Koef_lamb2[ele]=K1*B22A21+K2*B22A22
    #print(S)#er komplex
    #print(abs(S))
    Spec=np.real(S)
    #Spec=np.abs(S)**2
    #Spec=np.sqrt(np.real(S)**2+np.imag(S)**2) #dette er hvad abs gør
    Spec_norm=np.absolute(Spec)/np.max(np.absolute(Spec))
    spectrum_list[ele]=Spec_norm


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
    ###
    # fig, ax1 = plt.subplots(1,figsize=(11,7))
    # ax1.plot(omega_list,Spec_norm,'.')
    # ax1.grid()
    # plt.show()
spectrum_list=np.array(spectrum_list)
#
np.savez('./data/spectrum_formula_linewidth_g={}_gamD={}_gamA={}.npz'.format(g,gammaD,gammaA),kappa_list=kappa_list,linewidth_list=linewidth_list)
#
fig, ax1 = plt.subplots(1,figsize=(11,7))
ax1.loglog(kappa_list,linewidth_list/(2*np.pi))
ax1.grid()
#TITLE='{}-emitter. $g$={} ps$^-$$^1$. $\gamma_D$={} ps$^-$$^1$. $\gamma_A$={} ps$^-$$^1$. $P$={} ps$^-$$^1$'.format(N_em,g,gammaD,gamma,pump)
#ax1.set(title=TITLE)
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
#plt.show()
###
fig, ax1 = plt.subplots(1,figsize=(11,7))
##ax3
ax1.set(xlabel=r'$\mathrm{log}(\kappa)\, [\mathrm{ps^{-1}}]$') 
ax1.set(ylabel=r'$(\omega-\omega_{eg})\,\mathrm{[ps^{-1}]}$')
#ax3.set(ylim=[-1.5*g,1.5*g])
#ax3.set_xticks(ticks=np.arange(-3,4,1),labels=10^(np.arange(-3,4,1)))
#g=np.sqrt(N_em)*g
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(2)-1)),(spectrum_matrix).transpose(), 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss), wlist/(g*(np.sqrt(3)-np.sqrt(2))),(spectrum_matrix).transpose(), 20, cmap='RdGy')
ax1.contourf(np.log10(kappa_list)[:50], omega_list,(spectrum_list).transpose()[:,:50], 20, cmap='RdGy')
#plt.contourf(np.log10(pg_list_ss)[np.log10(pg_list_ss)<0], wlist,(spectrum_matrix[np.log10(pg_list_ss)<0]).transpose(), 20, cmap='RdGy')
plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap='RdGy'),ax=ax1)

fig, ax1 = plt.subplots(1,figsize=(11,7))
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(ylabel=r'$K1$')
ax1.semilogx(kappa_list,np.real(K1_list))
ax1.semilogx(kappa_list,np.imag(K1_list))
ax1.semilogx(kappa_list,np.absolute(K1_list))
ax1.legend(['Real','Imag','Abs'])
#ax3.set(ylim=[-1.5*g,1.5*g])

fig, ax1 = plt.subplots(1,figsize=(11,7))
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
ax1.set(ylabel=r'$K2$')
ax1.semilogx(kappa_list,np.real(K2_list))
ax1.semilogx(kappa_list,np.imag(K2_list))
ax1.semilogx(kappa_list,np.absolute(K2_list))
ax1.legend(['Real','Imag','Abs'])
print(K2_list[:])
print(K1_list[:])
fig, ax1 = plt.subplots(1,figsize=(11,7))
ax1.set(xlabel=r'$\kappa\,[\mathrm{ps^{-1}}]$',ylabel=r'$ \mathrm{FWHM}\,[\mathrm{THz}]$')
bum=4*g
bum=1
ax1.semilogx(kappa_list[kappa_list>bum],np.absolute(np.multiply(K2_list,B21A12_list)[kappa_list>bum]+np.multiply(K1_list,B21A11_list)[kappa_list>bum]))
ax1.semilogx(kappa_list[kappa_list>bum],np.absolute(np.multiply(K2_list,B22A22_list)[kappa_list>bum]+np.multiply(K1_list,B22A21_list)[kappa_list>bum]))
ax1.semilogx(kappa_list[kappa_list>bum],np.absolute(Koef_lamb1[kappa_list>bum]),'m.')
ax1.semilogx(kappa_list[kappa_list>bum],np.absolute(Koef_lamb2[kappa_list>bum]),'k.')
#ax1.semilogx(kappa_list[kappa_list>bum],)
#ax1.semilogx(kappa_list[kappa_list>bum],)
ax1.legend(['K1*B21A11+K2*B21A12','K1*B22A21+K2*B22A22','koeff1','koeff2'])
#ax1.legend(['B21A12','B22A22','B21A11','B22A21'])
ax1.set(xlim=[2,1000])
#print('koefficienter',np.multiply(K1_list,B21A12_list)[kappa_list>bum],np.multiply(K1_list,B22A22_list)[kappa_list>bum],np.multiply(K2_list,B21A11_list)[kappa_list>bum]*10**10,np.multiply(K2_list,B22A21_list)[kappa_list>bum]*10**10)
#ax3.set(ylim=[-1.5*g,1.5*g])
plt.show()