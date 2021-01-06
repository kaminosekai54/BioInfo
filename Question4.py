#  import of the requiered package
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functions import *

def main():
    args = [0.0082, 0.0149, 1, 0.3865, 600, 500, 4, 4, 0.04, 0.04, 0.1, 0.1, 0.002, 0.002]
    pT_vector, nullclineT, nullclineL, pL_vector, aux1, aux2 = draw_phase_space(args)

    plt.plot(pT_vector, nullclineT,label="Nullcline TetR")
    plt.plot(nullclineL,pL_vector,label="Nullcline LacI")
    plt.xlabel("LacI")
    plt.ylabel("TetR")
    plt.quiver(aux1[:,:,0],aux1[:,:,1], aux2[:,:,0], aux2[:,:,1])
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()