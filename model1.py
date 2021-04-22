import gurobipy as gp
from gurobipy import GRB

import numpy as np
m=gp.Model('model1')
n=50
w=6
p=5
sm=600
distr=np.zeros(11)
ct=np.zeros(11)
for i in range(0,len(distr)):
    distr[i]=0.1-0.002*i
    ct[i]=5*distr[i]*i
mu=sum(ct)
pary=mu-w

mu_coeff = list(range(0,n+5,p))
long_l=len(mu_coeff)
binvar=list(range(0,long_l))
so_coeff=[i**2 for i in mu_coeff]
Const_coeff=np.zeros((long_l,long_l+1))
for i in range(0,long_l):
    for j in range(0,(long_l-i)):
        if(j!=0):
            Const_coeff[i,long_l-j]=n-j*p + 5-5*i          
        else:
            Const_coeff[i,long_l-j]=mu-w
rest_coeff=np.ones((long_l,2))
for i in range(0,long_l):
    rest_coeff[i,1]=i*p
t=100000
y=m.addVars(binvar,name="y",obj=t,vtype=GRB.BINARY)


m.setObjective(t, GRB.MINIMIZE)

# m.addConstrs(gp.quicksum(distr[i]*mu_coeff[i] for i in range(len(mu_coeff)))==mu)
# m.addConstrs(gp.quicksum(distr[i]*so_coeff[i] for i in range(len(so_coeff)))==sm)
m.addConstrs(gp.quicksum(Const_coeff[i][j]*distr[j] for j in range(len(distr)) )+pary*y[i] <=mu for i in range(len(y)))
m.addConstr(gp.quicksum(y[i] for i in range(len(y)))==1)
m.addConstrs(gp.quicksum(t-5*i*y[i]>=0 for i in range(len(y))))
m.optimize()

for v in m.getVars():
    print("%s %s %8.2f %s %8.2f %s %8.2f %s %8.2f" % 
              (v.Varname, "=", v.X, ", reduced cost = ", abs(v.RC), ", from coeff = ", v.SAObjLow, "to coeff = ", v.SAObjUp))
    print(" ")



