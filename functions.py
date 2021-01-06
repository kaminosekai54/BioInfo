# this file will contain all the function used for the different test

# importing the required packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


#  function toggle_derivative
# This function compute the derivative wanted to study the toggle switch circuit
# @param
# @y, a list containing the initial concentration values for our reactances and product
# @t, a list of time range, mostly useful for the odeint function
# @args, a list containing all the different constant we need
def toggle_derivative (y, t, args):
	# initialising our value
    # storing the different value of args into variables
    # the variable are named according to their name in the equation
    #they represent the different rates and constant
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pT, gamma_pL, gamma_pT = args
    # storing the different value of y into variables
    # they are named according to their name into the equation
    # they represent the different initial concentration of our mRNA and proteins
    mL, mT, pL, pT = y

	# writing our derivative
    dmLdt = k_mLb+k_mLa*(1/(1+((pT**n_T)/(theta_T**n_T))))-gamma_mL*mL
    dmTdt = k_mTb+k_mTa*(1/(1+((pL**n_L)/(theta_L**n_L))))-gamma_mT*mT
    dpLdt = k_pL*mL - gamma_pL*pL
    dpTdt = k_pT*mT - gamma_pT*pT

    #  return of the computed derivative in a list
    return [dmLdt, dmTdt, dpLdt, dpTdt]


# function sol_lf
#returns the value of the QSS equations of mL for any given product value
# @param,
# @pT,product value 
# @stl, expression of the equation
def sol_lf(pT, stl):
	n=int(pT)
	stl=stl.subs("pT",n)
	return eval(str(stl))
# function sol_lf
#returns the value of the QSS equations of mT for any given product value
# @param,
# @pL, product value 
# @stt, expression of the equation
def sol_tf(pL, stt):
	n=int(pL)
	stt=stt.subs("pL",n)
	return eval(str(stt))

def draw_phase_space (args):
    
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pT, gamma_pL, gamma_pT = args

    #first number = start point, second number = end point, third number how many samples to generate with equal differences among them 
    x = np.linspace(0,2000,35)
    pL_vector = np.linspace(0,2000,35)
    pT_vector = np.linspace(0,2000,35)

    # creating our vector (the size is 35 because of the linspace before)
    mL_vector = np.zeros((35))
    mT_vector = np.zeros((35))

    # Return coordinate matrices from coordinate vectors.
    LacI_vector, TetR_vector = np.meshgrid(pL_vector, pT_vector)

    # getting the value of our derivative 
    for i in range(len(pT_vector)):
        mL_vector[i]  = (k_mLb+k_mLa*(1/(1+(pT_vector[i]**n_T/theta_T**n_T))))/gamma_mL
    
    for i in range (len(pL_vector)):
        mT_vector[i]  = (k_mTb+k_mTa*(1/(1+(pL_vector[i]**n_L/theta_L**n_L))))/gamma_mT
        
    #create matrices for our quiver plot
    aux1= np.zeros((35, 35, 2))
    aux2 = np.zeros((35, 35, 2))

    #fill the quiver matrixes, with the created vectors of proteins.
    for i in range(len(pL_vector)):
        for j in range(len(pT_vector)):
            mL, mT, pL, pT = toggle_derivative([mL_vector[j], mT_vector[i], pL_vector[i], pT_vector[j]], 0, args)
            aux1[i][j] = [pL_vector[i],pT_vector[j]]
            aux2[i][j] = [pL, pT]

    #get the points of the nullclines, using nullified differential equations for product formation
    nullclineL=[]
    for i in mL_vector:
        pl=(k_pL*i)/gamma_pL
        nullclineL.append(pl)
    nullclineT=[]
    for i in mT_vector:
        pt=(k_pT*i)/gamma_pT
        nullclineT.append(pt)

    return     [pT_vector, nullclineT, nullclineL, pL_vector, aux1, aux2]