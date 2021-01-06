import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functions import *
from matplotlib.widgets import Slider, Button, RadioButtons


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
    fig,g=plt.subplots()
    plt.subplots_adjust(bottom=0.35)

    axis_color = 'lightgoldenrodyellow'
    # Define an axes area and draw a slider in it
    pT_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.02], facecolor=axis_color)
    pT_slider = Slider(pT_slider_ax, 'pT', 0, 500000, valinit=0)

    pL_slider_ax  = fig.add_axes([0.25, 0.20, 0.65, 0.02], facecolor=axis_color)
    pL_slider = Slider(pL_slider_ax, 'pL', 0, 500000, valinit=0)

    t_slider_ax  = fig.add_axes([0.25, 0.25, 0.65, 0.02], facecolor=axis_color)
    t_slider = Slider(t_slider_ax, 't', 0, 10000, valinit=3000)

    KmTa_slider_ax  = fig.add_axes([0.25, 0.05, 0.65, 0.02], facecolor=axis_color)
    KmTa_slider = Slider(KmTa_slider_ax, 'Max. Inducible Rate TetR', 0.1, 1.5, valinit=0.3865)

    KmLa_slider_ax = fig.add_axes([0.25, 0.10, 0.65, 0.02], facecolor=axis_color)
    KmLa_slider = Slider(KmLa_slider_ax, 'Max. Inducible Rate LacI', 0.1, 1.5, valinit=1)

    l,=g.plot(model[:,2], model[:,3], color = trajectory_color)

    def sliders_on_changed(val):
        args[2]=KmLa_slider.val
        args[3]=KmTa_slider.val
        y0[2]=pL_slider.val
        y0[3]=pT_slider.val
        t = np.linspace(0,int(t_slider.val),int(t_slider.val)+1)
        model = odeint(toggle_derivative, y0, t, args = (args,)) # solve the differential equation 
        l.set_xdata(model[:,2])
        l.set_ydata(model[:,3])

        if model[:,2][-1] >1000 and model[:,3][-1] <100  : # LacI wins # 
            trajectory_color= 'g' # green
        
        elif model[:,3][-1] > 1000 and model[:,2][-1]<100: # TetR wins
            trajectory_color= 'b' # blue
            
        else:
            trajectory_color= 'k' # black

        l.set_color(trajectory_color)

    t_slider.on_changed(sliders_on_changed)
    pL_slider.on_changed(sliders_on_changed)
    pT_slider.on_changed(sliders_on_changed)
    KmTa_slider.on_changed(sliders_on_changed)
    KmLa_slider.on_changed(sliders_on_changed)

    reset_button_ax = fig.add_axes([0.05, 0.2, 0.1, 0.04])
    reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
    def reset_button_on_clicked(mouse_event):
        KmLa_slider.reset()
        KmTa_slider.reset()
        pL_slider.reset()
        pT_slider.reset()
        t_slider.reset()
    reset_button.on_clicked(reset_button_on_clicked)        
    g.set_xticks(np.linspace(0,2000,11)) #change the scale and the step
    g.set_yticks(np.linspace(0,2000,11)) #change the scale and the step
    g.set_xlabel('LacI')
    g.set_ylabel('TetR')
    g.set_title('LacI vs TetR concentration')
    plt.show()

if __name__ == '__main__':
    main()