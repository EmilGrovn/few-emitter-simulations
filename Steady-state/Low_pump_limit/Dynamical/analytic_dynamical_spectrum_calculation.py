import numpy as np
#consider an unpumped 1-emitter with parameters
g=1 #coupling constant [ps^(-1)]
kappa=1 #decay rate [ps^(-1)] for coupling from cavity to environment
gammaA=0
gammaD=0
gammaCE=(kappa+gammaA+gammaD)/2
gammaGE=(gammaA+gammaD)/2
#the Liouvillian matrix in a basis of {EE,CC,EC,CE,GG,GE,EG,GC,CG}
L=[[-gammaA ,      0,    1j*g,   -1j*g,      0,       0,       0,       0,       0],
   [0       , -kappa,   -1j*g,    1j*g,      0,       0,       0,       0,       0],
   [1j*g    ,  -1j*g,-gammaCE,       0,      0,       0,       0,       0,       0],
   [-1j*g   ,   1j*g,       0,-gammaCE,      0,       0,       0,       0,       0],
   [gammaA  ,  kappa,       0,       0,      0,       0,       0,       0,       0],
   [0       ,      0,       0,       0,      0,-gammaGE,       0,    1j*g,       0],
   [0       ,      0,       0,       0,      0,       0,-gammaGE,       0,   -1j*g],
   [0       ,      0,       0,       0,      0,    1j*g,       0,-kappa/2,       0],
   [0       ,      0,       0,       0,      0,       0,   -1j*g,       0,-kappa/2]]
Lambda,V=np.linalg.eig(L)
lambdaMatrix=np.diag(Lambda)
omega=0
omegaImatrix=np.diag(np.ones(9)/1j*omega)
Ainv=V
A=np.linalg.inv(V)
aDag=np.zeros((9,9))
aDag[8,4]=1
aDag[3,5]=1
aDag[1,7]=1
a=np.transpose(aDag)
print(a)#yes, it looks correct. Good check
print(np.matmul(aDag,a))#det er også korrekt
print(Lambda)
#we are ready to calculate the spectrum now
#let us do it in steps
rho0=np.zeros((9,9))
rho0[0,0]=1
#print(np.matmul(Ainv,np.matmul(np.linalg.inv(lambdaMatrix),np.matmul(A,rho0))))
# print(A)
# print(np.shape(A))
# print(A[5:,:5])
# print(A[:5,5:])
# print(A[:5,:5])
# print(A[5:,5:])
# print(Lambda)
#print(np.matmul(A,rho0))
#cannot inverse the eigenvalue matrix, since it has a zero-eigenvalue.
# we could perhaps add an extremely small
#let us start with calculating the 2-time correlation function then
#ah here we have to puild up a matrix with t and tau as well
#we are capable of doing that (easiest (in terms of programming) way is to do two for loops)
#then we can calculate the spectrum by numerical integration
lambdaMatrixReg=lambdaMatrix
lambdaMatrixReg[0,0]=1 #get rid of singularity
#my analytical considerations show that this eigenvalue does not matter anyway
#so we should be able to calculate the spectrum correctly numerically now
wlist=np.arange(-2*g,2*g+g/100,g/100)
Slist=np.zeros(len(wlist))+1j*np.zeros(len(wlist))
wlist=[g]
for i,omega in enumerate(wlist):
   Cbuild=np.matmul(A,rho0)
   #print('Step1',Cbuild[:,0])
   Cbuild=np.matmul(np.linalg.inv(lambdaMatrixReg),Cbuild)
   #print(np.linalg.inv(lambdaMatrixReg))
   #print('Step2',Cbuild[:5,:5])
   Cbuild=np.matmul(Ainv,Cbuild)
   Cbuild=np.matmul(a,Cbuild)
   C=np.matmul(A,Cbuild)
   D=np.matmul((np.linalg.inv(lambdaMatrixReg)-1j*omega),C)
   E=np.matmul(Ainv,D)
   F=np.matmul(aDag,E)
   Trace=F[0,0]+F[1,1]+F[4,4]
   Slist[i]=kappa*Trace

#print(Slist)
#I have calculated a formula (using that there are a lot of zeros in our matrices)
omega=0
S=0
for j in [5,6,7,8]:
   Cj2=A[1,0]/Lambda[1]*(A[j,5]*Ainv[3,1]+A[j,7]*Ainv[1,1])
   Dj2=Cj2/(Lambda[j]-1j*omega)
   S+=Ainv[7,j]*Dj2
   #test
   Cbuild=np.matmul(A,rho0)
   Cbuild=np.matmul(np.linalg.inv(lambdaMatrixReg),Cbuild)
   Cbuild=np.matmul(Ainv,Cbuild)
   Cbuild=np.matmul(a,Cbuild)
   Cbuild=np.matmul(A,Cbuild)
   #print("Cj2: analytical",Cj2,"\n numerical",Cbuild[j,1])
   #giver ikke det samme lige nu

S=kappa*S
#print(S)