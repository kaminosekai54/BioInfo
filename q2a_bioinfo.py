# import of the package
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functions import *

def main():
	args = [0.0082, 0.0149, 1, 0.3865, 600, 500, 4,4, 0.04, 0.04, 0.1, 0.1, 0.002, 0.002]
	y0 = [0, 0, 0, 0]                        #We decide the initial conditions for the enzyme concentration and the substrate concentration. At t=0, we know that there is no enzyme/substrate complex and no product yet.
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
	plt.show()                            #we plot the graph presenting the concentrations depending on time 

if __name__ == '__main__':
	main()

	