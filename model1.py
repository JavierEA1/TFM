from gurobipy import *
import numpy as np
m=model('model1')
n=50
w=6
p=5
sm=600
mu=20
pary=mu-w
distr=np.repeat(0.1,11)
mu_coeff = list(range(0,n+5,p))
long_l=len(mu_coeff)
binvar=list(range(1,long_l+1))

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
y=m.addVars(binvar,name="y",vtype=GRB.BINARY)


m.addConstrs(quicksum(distr[i]*mu_coeff[i] for i in mu_coeff)=mu)
m.addConstrs(quicksum(distr[i]*so_coeff[i] for i in so_coeff)=sm)
m.addConstrs(quicksum(Const_coeff[i][j]*distr[j] for j in distr )+pary*y[i] <=mu for i in y )
m.addConstrs(quicksum(y)=1)
m.addConstrs(obj-5*(i)*y[i]>0 for i in y)
m.setObjective(obj, GRB.MINIMIZE)
m.optimize()





