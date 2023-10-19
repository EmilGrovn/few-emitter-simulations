from qutip import *

a=basis(3,2)
b=basis(3,1)
c=basis(3,0)
print(tensor(a,c)+tensor(b,c))
print(tensor(a+b,c))
print(tensor(tensor(a+b),c))
print(tensor(tensor(a+b),tensor(c)))