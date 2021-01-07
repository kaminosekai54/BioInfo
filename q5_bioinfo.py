import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from matplotlib.widgets import Slider, Button, RadioButtons
from functions import*

#set parameters for integration and more
ar={'k_mLb': 0.0082, 'k_mTb': 0.0149, 'k_mLa': 1, 'k_mTa': 0.3865, 'theta_L': 600, 'theta_T': 500, 'n_L': 4, 'n_T': 4, 'gamma_mL': 0.04, 'gamma_mT': 0.04, 'k_pL': 0.1, 'k_pT': 0.1, 'gamma_pL': 0.002, 'gamma_pT': 0.002}
args = [0.0082, 0.0149, 1, 0.3865, 600, 500, 4,4, 0.04, 0.04, 0.1, 0.1, 0.002, 0.002]
t = np.linspace(0,3000,3001)

#general plotting commands
fig,g=plt.subplots(1,2)
plt.subplots_adjust(bottom=0.32)
plt.suptitle("Protein State Space")

#now for sliders!
axis_color = 'lightgoldenrodyellow'
# draw sliders for interesting parameters
KmTa_slider_ax  = fig.add_axes([0.25, 0.20, 0.65, 0.03], facecolor=axis_color)
KmTa_slider = Slider(KmTa_slider_ax, 'Max. Inducible Rate TetR', 0.1, 1.5, valinit=0.3865)

KmLa_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
KmLa_slider = Slider(KmLa_slider_ax, 'Max. Inducible Rate LacI', 0.1, 1.5, valinit=1)

gammap_slider_ax = fig.add_axes([0.25, 0.05, 0.65, 0.03], facecolor=axis_color)
gammap_slider = Slider(gammap_slider_ax, 'Degradation rate proteins', 0.0001, 0.05, valinit=0.002)

gammam_slider_ax = fig.add_axes([0.25, 0.10, 0.65, 0.03], facecolor=axis_color)
gammam_slider = Slider(gammam_slider_ax, 'Degradation rate mRNA', 0.01, 0.2, valinit=0.04)

#function to execute at any parameter change
#the recap boolean is used for the second figure plotted, which gathers everything on one subplot, with different quiver and colors for readability
def sliders_on_changed(val):
    if recap:
        g.clear()
        g.set_title("Trajectories on vector field")
        g.set_xlabel("LacI")
        g.set_ylabel("TetR")
    else:
        #clear both graphs
        g[0].clear()
        g[1].clear()
        #labeling
        g[0].set_title("Vector field and nullclines")
        g[1].set_title("Trajectories")
        for i in range(2):
            g[i].set_xlabel("LacI")
            g[i].set_ylabel("TetR")
    #important part, assigning slider value to appropriate arguments
    args[2]=KmLa_slider.val
    args[3]=KmTa_slider.val
    args[-1],args[-2]=gammap_slider.val,gammap_slider.val
    args[-5],args[-6]=gammam_slider.val,gammam_slider.val
    #fetching arguments from list for later use
    k_mLb, k_mTb, k_mLa, k_mTa, theta_L, theta_T, n_L, n_T, gamma_mL, gamma_mT, k_pL, k_pT, gamma_pL, gamma_pT = args
    #draw phase space for arguments set by sliders
    ps=draw_phase_space(args)
    #plot the phase space and nullclines
    if recap:
        g.quiver(ps[-2][:,:,0],ps[-2][:,:,1], ps[-1][:,:,0], ps[-1][:,:,1],minshaft=0.1,minlength=0.3,headlength=2,headaxislength=2,headwidth=3,alpha=0.4,width=0.002,linestyle='solid')
    else:
        g[0].quiver(ps[-2][:,:,0],ps[-2][:,:,1], ps[-1][:,:,0], ps[-1][:,:,1])
        g[0].plot(ps[2],ps[0],'g',label="LacI nullcline")
        g[0].plot(ps[3],ps[1],'c',label="TetR nullcline")

    #now for the second subplot, being protein trajectories
    #initiate vectors
    LacI_vector= np.linspace(0, 2000, 11)
    TetR_vector= np.linspace(0, 2000, 11)
    #fill those vectors using mRNA QSS equation
    for i in range(len(LacI_vector)):
        for j in range(len(TetR_vector)):
            mT=  (k_mTb+k_mTa*(1/(1+(LacI_vector[i]**n_L/theta_L**n_L))))/gamma_mT
            mL=  (k_mLb+k_mLa*(1/(1+(TetR_vector[j]**n_T/theta_T**n_T))))/gamma_mL
            
            y0 = [mL, mT,LacI_vector[i],TetR_vector[j]]

            #integrate with new arguments
            model = odeint(toggle_derivative, y0, t, args=(args,))

            #conditions defining curve color depending on final state of trajectory
            if model[:,2][-1] > 1000 and model[:,3][-1] < 250  : # LacI wins
                    trajectory_color= 'g' # green
    
            elif model[:,3][-1] > 1000 and model[:,2][-1] < 250: # TetR wins
                    trajectory_color= 'b' # blue
        
            else:
                    trajectory_color= 'k' # black

            #plot the new integrated values
            if recap:
                g.plot(model[:,2], model[:,3], color=trajectory_color,linewidth=0.75)
            else:
                g[1].plot(model[:,2], model[:,3], color=trajectory_color)
    if recap:
        g.plot(ps[2],ps[0],'g',label="LacI nullcline",linewidth=3,color="r")
        g.plot(ps[3],ps[1],'c',label="TetR nullcline",linewidth=3,color="#b938ff")
        g.legend()
    else:
        g[0].legend()
recap=False
#execute previously defined function when sliders are changed
KmTa_slider.on_changed(sliders_on_changed)
KmLa_slider.on_changed(sliders_on_changed)
gammam_slider.on_changed(sliders_on_changed)
gammap_slider.on_changed(sliders_on_changed)

# button for resetting the parameters
reset_button_ax = fig.add_axes([0.005, 0.025, 0.05, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    KmLa_slider.reset()
    KmTa_slider.reset()
    gammam_slider.reset()
    gammap_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)
#final pyplot commands
plt.legend()
plt.show()

#now for a big recap graph :
recap=True
fig,g=plt.subplots()
plt.subplots_adjust(bottom=0.32)
plt.title("Protein State Space")

#now for sliders!
axis_color = 'lightgoldenrodyellow'
# draw sliders for interesting parameters
KmTa_slider_ax  = fig.add_axes([0.25, 0.20, 0.65, 0.03], facecolor=axis_color)
KmTa_slider = Slider(KmTa_slider_ax, 'Max. Inducible Rate TetR', 0.1, 1.5, valinit=0.3865)

KmLa_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
KmLa_slider = Slider(KmLa_slider_ax, 'Max. Inducible Rate LacI', 0.1, 1.5, valinit=1)

gammap_slider_ax = fig.add_axes([0.25, 0.05, 0.65, 0.03], facecolor=axis_color)
gammap_slider = Slider(gammap_slider_ax, 'Degradation rate proteins', 0.0001, 0.05, valinit=0.002)

gammam_slider_ax = fig.add_axes([0.25, 0.10, 0.65, 0.03], facecolor=axis_color)
gammam_slider = Slider(gammam_slider_ax, 'Degradation rate mRNA', 0.01, 0.2, valinit=0.04)

#execute previously defined function when sliders are changed
KmTa_slider.on_changed(sliders_on_changed)
KmLa_slider.on_changed(sliders_on_changed)
gammam_slider.on_changed(sliders_on_changed)
gammap_slider.on_changed(sliders_on_changed)

# button for resetting the parameters
reset_button_ax = fig.add_axes([0.005, 0.025, 0.05, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    KmLa_slider.reset()
    KmTa_slider.reset()
    gammam_slider.reset()
    gammap_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)
#final pyplot commands
plt.legend()
plt.show()