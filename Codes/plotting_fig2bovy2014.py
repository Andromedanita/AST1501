# Anita Bahmanyar
# March 2016, U of T
# This code takes directory and filenames and plots fig. 2 in Bovy 2014 subplots
# So it basically returns the action angle coordinates

from nemo_funcs        import *
from galpy.actionAngle import actionAngleStaeckel, actionAngleIsochroneApprox
from galpy.potential   import LogarithmicHaloPotential, MWPotential2014
from galpy.util        import bovy_plot, bovy_coords, bovy_conversion



def aA_prog(ptype):
    """
    Parameters
    ----------------------------------------------------------------------------------
        ptype : potential type. It can be either "Log" for logarithmic potential or
                "MW" for Milky Way 2014 potential.
        
        
    Return
    ----------------------------------------------------------------------------------
        action-angle values of the progenitor for GD-1 stream given the potentia type
    """
    
    # setting up the action-angle instance
    if ptype == "Log":
        lp         = potential.LogarithmicHaloPotential(q=0.9, normalize=1.)
        aAS        = actionAngleIsochroneApprox(pot = lp, b=0.8)

    elif ptype == "MW":
        aAS        = actionAngleIsochroneApprox(pot = MWPotential2014, b=0.6, tintJ=1000,ntintJ=2000)
    
    # current position and velocity of the GD-1 stream progenitor in cartesian coordinates
    x, y, z    = np.array([12.4, 1.5, 7.1])
    vx, vy, vz = np.array([107., -243., -105.])

    # current position and velocity in cylindrical coordinate
    R, phi, zz = bovy_coords.rect_to_cyl(x, y, z)
    vR, vT, vz = bovy_coords.rect_to_cyl_vec(vx, vy, vz, x, y, z, cyl = False)
    
    ro = 8.0
    vo = 220.0

    # converting to galpy natural units
    R  /= ro
    zz /= ro
    vR /= vo
    vT /= vo
    vz /= vo

    # computing the action-angle coordinates
    val = aAS.actionsFreqsAngles(R, vR, vT, zz, vz, phi)

    return val



if __name__ == '__main__':
    
    # initializing the data
    
    # color lists
    clist_Log     = ["black", "blue", "green", "red", "yellow"]
    clist_MW2e4   = ["black", "blue", "green", "red"]

    # filename lists
    fnames_Log2e4 = ["concatenated_W2.5_M2e4.txt","concatenated_W3_M2e4.txt","concatenated_W6_M2e4.txt","concatenated_W9_M2e4.txt","concatenated_W12_M2e4.txt"]
    fnames_MW2e4  = ["concatenated_W2_M2e4_MWpot_orphanway.dat", "concatenated_W6_M2e4_MWpot_orphanway.dat", "concatenated_W9_M2e4_MWpot_orphanway.dat", "concatenated_W12_M2e4_MWpot_orphanway.dat"]
    fnames_Log2e6 = ["concatenated_W2_M2e6.txt","concatenated_W3_M2e6.txt","concatenated_W6_M2e6.txt","concatenated_W9_M2e6.txt","concatenated_W12_M2e6.txt"]
    fnames_Log2e7 = ["concatenated_W2_M2e7.txt","concatenated_W3_M2e7.txt","concatenated_W6_M2e7.txt","concatenated_W9_M2e7.txt","concatenated_W12_M2e7.txt"]

    # W0 values lists
    W0_Log2e4 = ["2.5","3.0","6.0","9.0","12.0"]
    W0_MW2e4  = ["2.0","6.0","9.0","12.0"]
    W0_Log2e6 = ["2.0","3.0","6.0","9.0","12.0"]
    W0_Log2e7 = ["2.0","3.0","6.0","9.0","12.0"]

    # directories
    dirname_Log = "/Users/anita/Documents/Grad_First_Year/Research/data_files/LOGpot/"
    dirname_MW  = "/Users/anita/Documents/Grad_First_Year/Research/data_files/MWpot/"


def plot_fig2(dirname, filename, plot_type):
    """
    Parameters
    ----------------------------------------------------------------------------------
        dirname : directory where action-angle files exist
        
        filename: action-angle filename
        
        plot_type: set of values returned by the function plot_fig2 to be plotted
        
    
    Return
    ----------------------------------------------------------------------------------
        set of values specified. It can be "OzOr" or "JzJr" for instance.
    """
    
    tail = np.loadtxt(dirname+filename)
    JR_t, Lz_t, Jz_t, omegaR_t, omega_phi_t, omegaz_t, thetaR_t, theta_phi_t, thetaz_t \
    = tail[:,0], tail[:,1], tail[:,2], tail[:,3], tail[:,4], tail[:,5], tail[:,6], tail[:,7], tail[:,8]
    indx_t = tail_cut(tail)
    JR_t, Lz_t, Jz_t                  = tail[indx_t,0], tail[indx_t,1], tail[indx_t,2]
    omegaR_t, omega_phi_t, omegaz_t   = tail[indx_t,3], tail[indx_t,4], tail[indx_t,5]
    thetaR_t, theta_phi_t, thetaz_t   = tail[indx_t,6], tail[indx_t,7], tail[indx_t,8]
    
    thetaR_t    = (np.pi+(thetaR_t-np.median(tail[:,6]))) % (2.*np.pi)
    theta_phi_t = (np.pi+(theta_phi_t-np.median(tail[:,7]))) % (2.*np.pi)
    thetaz_t    = (np.pi+(thetaz_t-np.median(tail[:,8]))) % (2.*np.pi)
    
    if plot_type == "OzOr":
        return omegaR_t, omegaz_t
    elif plot_type == "OphiOr":
        return omegaR_t, omega_phi_t

    elif plot_type == "JzJr":
        return JR_t, Jz_t

    elif plot_type == "LzJr":
        return JR_t, Lz_t

    elif plot_type == "AzAr":
        return thetaR_t, thetaz_t

    elif plot_type == "AphiAr":
        return thetaR_t, theta_phi_t


def plot_all(ptype):
    
    """
    Parameters
    ----------------------------------------------------------------------------------
        ptype : type of simulation. For example if it is Log2e4 (Logarithmic potential
                with mass of M = 2 X 10^4 or MW2e4 (Milky Way 2014 potential with the
                same mass)
        
        
    Return
    ----------------------------------------------------------------------------------
        Does not return anything. It only plots the specified set of variables set
        by plot_type argument.
    """
    
    plt.ion()

    fact  = bovy_conversion.freq_in_Gyr(220., 8.)

    if ptype == "Log2e4":
        val_prog = aA_prog("Log")
        plt.plot(val_prog[3]*fact, val_prog[5]*fact, "mo", ms=10, label="progenitor")
        i = 0
        for names in fnames_Log2e4:
            omegaR_t, omegaz_t = plot_fig2(dirname_Log,names,"OzOr")
            plt.plot(omegaR_t*fact, omegaz_t*fact, 'o', color=clist_Log[i],ms=2,label="$W_0= $" + W0_Log2e4[i])
            plt.title(r"Log potential, $M = 2 \times 10^4$")
            i += 1

    elif ptype == "MW2e4":
        val_prog = aA_prog("MW")
        plt.plot(val_prog[3]*fact, val_prog[5]*fact, "mo", ms=10, label="progenitor")
        i = 0
        for names in fnames_MW2e4:
            omegaR_t, omegaz_t = plot_fig2(dirname_MW,names, "OzOr")
            plt.plot(omegaR_t*fact, omegaz_t*fact,'o',  color=clist_MW2e4[i],label="$W_0= $" + W0_MW2e4[i])
            plt.title(r"MW2014 potential, $M = 2 \times 10^4$")
            i += 1

    elif ptype == "Log2e6":
        val_prog = aA_prog("Log")
        plt.plot(val_prog[3]*fact, val_prog[5]*fact, "mo", ms=10, label="progenitor")
        i = 0
        for names in fnames_Log2e6:
            omegaR_t, omegaz_t = plot_fig2(dirname_Log,names, "OzOr")
            plt.plot(omegaR_t*fact, omegaz_t*fact, 'o', ms=2, color=clist_Log[i],label="$W_0= $" + W0_Log2e6[i])
            plt.title(r"Log potential, $M = 2 \times 10^6$")
            i += 1

    elif ptype == "Log2e7":
        val_prog = aA_prog("Log")
        plt.plot(val_prog[3]*fact, val_prog[5]*fact, "mo", ms=10, label="progenitor")
        i = 0
        for names in fnames_Log2e7:
            omegaR_t, omegaz_t = plot_fig2(dirname_Log,names, "OzOr")
            plt.plot(omegaR_t*fact, omegaz_t*fact, 'o', ms = 2, color=clist_Log[i],label="$W_0= $" + W0_Log2e7[i])
            plt.title(r"Log potential, $M = 2 \times 10^7$")
            i += 1

    plt.xlabel("$\Omega_R (\mathrm{Gyr^{-1}})$", fontsize=20)
    plt.ylabel("$\Omega_z (\mathrm{Gyr^{-1}})$", fontsize=20)
    plt.legend(loc='best')

