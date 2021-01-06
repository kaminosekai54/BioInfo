import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint




# Helper function 1: ODE model for the Toggle Switch

def toggle_derivative (y, t, args):
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pT, gamma_pL, gamma_pT = args
    mL, mT, pL, pT = y

# Helper function 2: plot the behavior of the system in the protein state space

    dmLdt = k_mLb+k_mLa*(1/(1+((pT**n_L)/(theta_L**n_L))))-gamma_mL*mL
    dmTdt = k_mTb+k_mTa*(1/(1+((pL**n_T)/(theta_T**n_T))))-gamma_mT*mT
    dpLdt = k_pL*mL - gamma_pL*pL
    dpTdt = k_pT*mT - gamma_pT*pT
                         
    return  [dmLdt, dmTdt, dpLdt, dpTdt]

args = [0.0082, 0.0149, 1, 0.3865, 600, 500, 4,4, 0.04, 0.04, 0.1, 0.1, 0.002, 0.002]
y0 = [0, 0, 0, 500]                        #We randomly decide the initial conditions for the enzyme concentration and the substrate concentration. At t=0, we know that there is no enzyme/substrate complex and no product yet.
t = np.linspace(0,3000,3001)                   #set the timespan for the simulation 
model = odeint(toggle_derivative, y0, t, args = (args,)) # solve the differential equation 

plt.plot(t, model[:,0], 'r', label='mL concentration')
plt.plot(t, model[:,1], 'b', label='mT concentration')
plt.plot(t, model[:,2], 'k', label='pL concentration')
plt.plot(t, model[:,3], 'y', label='pT concentration')
plt.xlabel('time')
plt.ylabel('concentrations')
plt.title('evolution of the mRNA and protein concentrations')
plt.legend()
plt.show()  
    
    #we plot the graph presenting the concentrations depending on time 

#I compared to a 1000 and a 100 but you can choose your numbers there

if model[:,2][-1] >1000 and model[:,3][-1] < 1000  : # LacI wins # 
    trajectory_color= 'g' # green
    
elif model[:,3][-1] > 1000 and model[:,2][-1]< 1000: # TetR wins
    trajectory_color= 'b' # blue
        
else:
    trajectory_color= 'k' # black
            
plt.plot(model[:,2], model[:,3], color = trajectory_color)
plt.xticks(np.linspace(0,200,11)) #change the scale and the step
plt.yticks(np.linspace(0,200,11)) #change the scale and the step
plt.xlabel('LacI')
plt.ylabel('TetR')
plt.title('LacI vs TetR concentration')
plt.show()
