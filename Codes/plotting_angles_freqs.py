import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from   plotting_fig2bovy import *

plt.ion()

fact = bovy_conversion.freq_in_Gyr(220.,8.)
dir = "/Users/anita/Documents/Grad_First_Year/Research/data_files/MWpot/"

value = aA_prog("MW")

# moved position of the progenitor in angles
val6 = (np.pi+(value[6]-np.median(value[6]))) % (2.*np.pi)
val8 = (np.pi+(value[8]-np.median(value[8]))) % (2.*np.pi)

def plotter(filename, ptype):

    data = np.loadtxt(dir+filename)
    
    # frequencies
    omegaR = data[:,3]
    omegaz = data[:,5]
    
    # angles
    thetaR = data[:,6]
    thetaz = data[:,8]
    
    thetaR_t = (np.pi+(thetaR-np.median(data[:,6]))) % (2.*np.pi)
    thetaz_t = (np.pi+(thetaz-np.median(data[:,8]))) % (2.*np.pi)
    
    colors = cm.rainbow(np.linspace(0, 1, len(omegaR)))
    
    if ptype == "freqs":
    
        plt.scatter(omegaR*fact,omegaz*fact, color=colors)
        
        # progenitor location
        plt.plot(value[3]*fact,value[5]*fact,"ko", ms=12)
        
        plt.xlabel(r"$\Omega_R \mathrm{(Gyr^{-1})}$",fontsize=20)
        plt.ylabel(r"$\Omega_z \mathrm{(Gyr^{-1})}$",fontsize=20)
        plt.grid(which='minor')
        plt.minorticks_on()
    
    
    elif ptype == "angles":
        # moved stream particles
        plt.scatter(thetaR_t,thetaz_t, color=colors)
    
        # progenitor location
        plt.plot(val6,val8,"ko", ms=12)
    
        plt.grid(which='minor')
        plt.minorticks_on()
        plt.xlabel(r"$\theta_R}$",fontsize=20)
        plt.ylabel(r"$\theta_z}$",fontsize=20)