#Just testing the implementation of the stochastic test. 

#We decided to make a stochastic model to show that if a system is bistable and initial conditions
#are near separatrix (unstable equilibrium point) with the same initial conditions some trajectories will be
#blue and some will be green. This is due to the fact that small fluctuations are able to push our system
#over the line and switch to one of the two states.


#importing packages

import matplotlib.pyplot as plt
import numpy as np
import gillespy2

#Stochastic_toggle

class Toggle(gillespy2.Model):

    def __init__(self, parameter_values=None):

        # First call the gillespy2.Model initializer.

        gillespy2.Model.__init__(self, name='Toggle')



        # Define parameters for the rates of reactions
        
        #We chose to multiply K_mTa by 3 to make our system bistable (with the sliders it seemed to be the optimal parameter)
        
        k_mLb = gillespy2.Parameter(name='k_mLb', expression="0.0082")

        k_mTb = gillespy2.Parameter(name='k_mTb', expression="0.0149")
        
        k_mLa = gillespy2.Parameter(name='k_mLa', expression="1")
        
        k_mTa = gillespy2.Parameter(name='k_mTa', expression="0.3865*3")
        
        theta_L = gillespy2.Parameter(name='theta_L', expression = "600")
        
        theta_T = gillespy2.Parameter(name='theta_T', expression = "500")
        
        n_L = gillespy2.Parameter(name='n_L', expression = "4")
        
        n_T = gillespy2.Parameter(name='n_T', expression = "4")
        
        k_pL = gillespy2.Parameter(name='k_pL', expression = 0.1)
        
        k_pT = gillespy2.Parameter(name='k_pT', expression = 0.1)

        gamma_mL= gillespy2.Parameter(name='gamma_mL', expression = 0.04)
        
        gamma_mT = gillespy2.Parameter(name='gamma_mT', expression = 0.04)
        
        gamma_pT = gillespy2.Parameter(name='gamma_pT', expression = 0.002)
        
        gamma_pL = gillespy2.Parameter(name='gamma_pL', expression = 0.002)
        
        self.add_parameter([k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, k_pL, k_pT, gamma_mL, gamma_mT, gamma_pT, gamma_pL])



        # Define variables for mRNAs and proteins
        #Choosing initial values. 
        
        mL = gillespy2.Species(name='mL', initial_value=13)

        mT = gillespy2.Species(name='mT',   initial_value=9)
        
        pL = gillespy2.Species(name ='pL', initial_value=750)
        
        pT = gillespy2.Species(name ='pT', initial_value=500)

        self.add_species([mL, mT, pL, pT])



        # Define the reactions

        r_mL = gillespy2.Reaction(name="rmL_creation", reactants={}, products={mL:1}, propensity_function = "k_mLb+k_mLa*(1/(1+((pT**n_T)/(theta_T**n_T))))")

        r_mT = gillespy2.Reaction(name="rmT_creation", reactants={}, products={mT:1}, propensity_function = "k_mTb+k_mTa*(1/(1+((pL**n_L)/(theta_L**n_L))))")
        
        rd_mL = gillespy2.Reaction(name="rmL_degradation", rate=gamma_mL, reactants={mL:1}, products={})

        rd_mT = gillespy2.Reaction(name="rmT_degradation", rate=gamma_mT, reactants={mT:1}, products={})
        
        r_pL = gillespy2.Reaction(name="rpL_creation", rate=k_pL, reactants={mL:1}, products={mL:1, pL:1})

        r_pT = gillespy2.Reaction(name="rpT_creation", rate=k_pT, reactants={mT:1}, products={mT:1, pT:1})
        
        rd_pL = gillespy2.Reaction(name="rpL_degradation", rate=gamma_pL, reactants={pL:1}, products={})
        
        rd_pT = gillespy2.Reaction(name="rpT_degradation", rate=gamma_pT, reactants={pT:1}, products={})


        self.add_reaction([r_mL, r_mT, rd_mL, rd_mT, r_pL, r_pT, rd_pL, rd_pT])



        # Set the timespan for the simulation.

        self.timespan(np.linspace(0, 3000, 3001))



model = Toggle()

results = model.run(number_of_trajectories=10)

#checking the winning trajectory 

for index in range(0, 10):

    trajectory = results[index]
    
    if trajectory['pL'][-1] > 1000 and trajectory['pT'][-1] < 1000  : # LacI wins
        trajectory_color= 'g' # green
    
    elif trajectory['pT'][-1] > 1000 and trajectory['pL'][-1] < 1000: # TetR wins
        trajectory_color= 'b' # blue
        
    else:
        trajectory_color= 'k' # black
        
    plt.plot(trajectory['pL'], trajectory['pT'], color = trajectory_color)

#plotting

plt.xticks(np.linspace(0,2000,11))
plt.yticks(np.linspace(0,2000,11))
plt.xlabel('LacI')
plt.ylabel('TetR')
plt.title('LacI vs TetR concentration')
plt.show()