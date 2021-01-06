import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from sympy.solvers import solve
from sympy import Symbol

#useful dict for later use, containing arguments and their values
ar={'k_mLb': 0.0082, 'k_mTb': 0.0149, 'k_mLa': 1, 'k_mTa': 0.3865, 'theta_L': 600, 'theta_T': 500, 'n_L': 4, 'n_T': 4, 'gamma_mL': 0.04, 'gamma_mT': 0.04, 'k_pL': 0.1, 'k_pT': 0.1, 'gamma_pL': 0.002, 'gamma_pT': 0.002}
ary={'mL':0,'mT':0,'pL':0,'pT':0}
#function that returns the differential equations at hand
def toggle_derivative (y,t,args):
    
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pT, gamma_pL, gamma_pT = args
    
    mL, mT, pL, pT = y
    
    dmLdt = k_mLb+k_mLa*(1/(1+((pT**n_T)/(theta_T**n_T))))-gamma_mL*mL
    
    dmTdt = k_mTb+k_mTa*(1/(1+((pL**n_L)/(theta_L**n_L))))-gamma_mT*mT
    
    dpLdt = k_pL*mL - gamma_pL*pL
    
    dpTdt = k_pT*mT - gamma_pT*pT
    
    return [dmLdt, dmTdt, dpLdt, dpTdt]
#define arguments in tuples for odeint
y=(0,0,0,0)
args=(0.0082,0.0149,1,0.3865,600,500,4,4,0.04,0.04,0.1,0.1,0.002,0.002)
#generate values for t axis
t=np.linspace(0,3000,3001)
#integrate our diff. eq. for our parameters
model=odeint(toggle_derivative,y,t,args=(args,))
#plot the integration results, both for mRNA and protein, on two different subplots
fig, g=plt.subplots(1,2)
g[0].plot(t,model[:,0],label="LacI")
g[0].plot(t,model[:,1],label="TetR")
g[1].plot(t,model[:,2],label="LacI")
g[1].plot(t,model[:,3],label="TetR")
#labeling
plt.suptitle("Inhibition dynamics")
g[0].set_ylabel("copy numbers")
g[1].set_ylabel("copy numbers")
g[0].set_title("mRNA")
g[1].set_title("Product")
for i in range(2):
	g[i].set_xlabel("Time")

#now for the QSS solving using sympy
#set the various parameters as sympy symbols

for i in ar.keys():
	globals()[i]=Symbol(str(i))
for i in ary.keys():
	globals()[i]=Symbol(str(i))

#actual solving for dmldt and dmtdt null
sol_l=solve(k_mLb+k_mLa*(1/(1+((pT**n_T)/(theta_T**n_T))))-gamma_mL*mL,mL)
sol_t=solve(k_mTb+k_mTa*(1/(1+((pL**n_L)/(theta_L**n_L))))-gamma_mT*mT,mT)

stl=sol_l[0]
stt=sol_t[0]

#substitute the parmaters for their given value
for i in ar.keys():
	stl=stl.subs(i,ar[i])
	stt=stt.subs(i,ar[i])

#functions that returns the value of the QSS equations of mL and mT for any given product value
def sol_lf(pT,stl):
	n=int(pT)
	stl=stl.subs("pT",n)
	return eval(str(stl))
def sol_tf(pL,stt):
	n=int(pL)
	stt=stt.subs("pL",n)
	return eval(str(stt))

#create lists of values of mL* and mT*, using the sol functions defined previously, on the odeint data
mstarl=[]
for i in model[:,3]:
	ii=sol_lf(i,stl)
	mstarl.append(float(ii))
mstart=[]
for i in model[:,2]:
	ii=sol_tf(i,stt)
	mstart.append(float(ii))

#plot the mL* and mT* values on the mRNA graph
g[0].plot(t,mstarl,label="LacI (QSS)")
g[0].plot(t,mstart,label="TetR (QSS)")

#show our figure
g[0].legend()
g[1].legend()
plt.show()