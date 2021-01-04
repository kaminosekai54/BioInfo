import numpy as np
import matplotlib.pyplot as plt

# Helper function 1: ODE model for the Toggle Switch

def toggle_derivative (y,t,args):
    
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pL, gamma_pL, gamma_pT = args
    
    mL, mT, pL, pT = y
    
    dmLdt = ?
    
    dmTdt = ?
    
    dpLdt = ?
    
    dpTdt = kTp*mT - gammaTp*pT
    
    return [dmLdt, dmTdt, dpLdt, dpTdt]

# Helper function 2: plot the behavior of the system in the protein state space

def draw_phase_space (args):
    
    
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pL, gamma_pL, gamma_pT = args
    
    
    ### Part 1: draw LacI/TetR vector field
    
    pL_vector = np.linspace(0,2000,35)
    pT_vector = np.linspace(0,2000,35)
    
    
    LacI_vector, TetR_vector = np.meshgrid(pL_vector, pT_vector)
    
    # expressions for steady-state concentrations of mRNAs: mL and mT
    
    mL_vector = ?
    
    mT_vector = ?
    
    
    # compute the value of the derivative of pL and pT for each pair of values in pL_vector and pT_vector
    aux1 = ?
    aux2 = ?
    counter = ?
    
    for i in range(?):
        for j in range(?):
            mL, mT, pL, pT = toggle_derivative([?, ?, ?,?],0,args)
            aux1[counter] = ?
            aux2[counter] = ?
            counter+=1


    # plot the vector field in the protein state space
    plt.quiver(LacI_vector, TetR_vector, aux1, aux2)

    ## Part 2: draw the LacI and TetR nullclines in the protein state space, using a steady state approximation for mRNAs
    
    LacI_vector= np.linspace(0,2000,300)
    TetR_vector= np.linspace(0,2000,300)

    # expression of LacI and TetR nullclines
    nLacI_vector = ?
    nTetR_vector = ?

    # plot nullclines
    plt.plot(nLacI_vector,TetR_vector,'g')
    plt.plot(LacI_vector,nTetR_vector,'c')


### Main part of project here:

    # Set toggle switch parameters here:

    # Q2: Perform numerical simulation here

    # Q2 a): Plot the protein and mRNA concentrations as a function of time here

    # Q2 b): Compute and plot steady-state concentrations of mRNA

    # Q2 c): Plot the system trajectory in the protein space

    # Q4: Plot the vector field in the protein space, and draw the protein nullclines, under the assumption that mRNAs are at steady-state
    draw_phase_space (?)

    # Q5: Plot in the state space system trajectories for which initial conditions are located on an 11x11 regular grid.
    LacI_vector= linspace(0,2000,11);
    TetR_vector= linspace(0,2000,11);
