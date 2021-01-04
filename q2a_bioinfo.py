# importation of the needed package

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


#  function toggle_derivative
# This function compute the derivative wanted to study the biologic circuit
# @param
# @y, a list containing the initial values for our reactances and product
# @t, a list of time range, mostly useful for the odeint function
# @args, a list containing all the different constant we need
def toggle_derivative (y, t, args):
	# initialising our value
	k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pT, gamma_pL, gamma_pT = args
	mL, mT, pL, pT = y

	# writing our derivative
	dmLdt = k_mLb+k_mTa*(1/(1+(pT*n_T)/(theta_T*n_T)))-gamma_mL*mL
	dmTdt = k_mTb+k_mLa*(1/(1+(pL*n_L)/(theta_L*n_L)))-gamma_mT*mT
	dpLdt = k_pL*mL - gamma_pL*pL
	dpTdt = k_pT*mT - gamma_pT*pT

                        #  just returning a list containing our different derivative expression
	return  [dmLdt, dmTdt, dpLdt, dpTdt]
# Main function
def main():
	#initialisation of our parameter
	args = [0.0082, 0.0149, 1, 0.3865, 600, 500, 4,4, 0.04, 0.04, 0.1, 0.1, 0.002, 0.002] # our different constant value
	y0 = [0, 0, 0, 0]                        # our initial concentrations
	t = np.linspace(0,3000,3001)                   #set the timespan for the simulation 

	# creating our model with odeint
	model = odeint(toggle_derivative, y0, t, args = (args,)) # solve the differential equation

# ploting our results
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