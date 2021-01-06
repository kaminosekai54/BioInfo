import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functions import *


def main():
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
    #  The "winner" is determined if one protein is over 1000 and the other protein is under 1000
    # model[:,2][-1], with that, we access to the 3rd list of the model, containing a sub list that contain all the value of concentration, so we get the last element of this sublist to get the final concentration
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

if __name__ == '__main__':
    main()