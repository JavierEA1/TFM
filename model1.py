from gurobipy import *
import numpy as np
m=model('model1')
n=50
w=6
p=5
mu=20
mu_coeff = list(range(0,n+5,p))
long_l=len(mu_coeff)
so_coeff=[i**2 for i in mu_coeff]
Const_coeff=np.zeros((long_l-1,long_l+1))
for i in range(0,long_l-1):
    for j in range(0,(long_l-i)):
        if(j!=0):
            Const_coeff[i,long_l-j]=n-j*p + 5-5*i          
        else:
            Const_coeff[i,long_l-j]=mu-w
rest_coeff=np.ones((long_l,2))
for i in range(0,long_l):
    rest_coeff[i,1]=i*p

min(t)



