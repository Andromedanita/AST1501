import numpy            as np
import matplotlib.pylab as plt
import os



def run_nemo(output_name,num_part, w0, mass, rt, wd_units, output_shifted, xs, ys, zs, vxs, vys, vzs, output_evol, tstop, eps, step, kmax, Nlev, fac, accname, accparse, output_final):

    """
    Parameters:
    -------------------------------------------------------
        output_name    : name of the output file
        
        num_part       : number of nbody particles
        
        w0             : dimensionless central potential in the King
                         profile
        
        mass           : total mass of the object
        
        rt             : tidal radius (where density -> 0)
        
        wd_units       : units of time
        
        output_shifted : name of the shifted output file
        
        xs, ys, zs     : shift in the position
        
        vxs, vys, vzs  : shift in the velocities
        
        output_evol    : name of evoluted output file
        
        tstop          : stop time
        
        eps            : softening length
        
        step           : time step
        
        kmax           : the longest timestep is taken to be 2^(-kmax) simulation time units
        
        Nlev           : 
        
        fac            : These factors control the average time step of a body to be the minimum
        
        accname        : if given, an external acceleration field wuth the specified name is used (multiple fields are
                         also accepted in the form of name1 + name2 + name3
        
        accparse       : parameters of accname (should be a list or an array). If multiple fields are given in accname,
                         the parameters should be given in the form of param1 ; param2 ; param3
        
        output_final   : name of the final output file
        
    
    Returns:
    -------------------------------------------------------
        Runs NEMO for King's profile using the specified 
        parameters. It then shifts the position and velocity 
        to its start time.
        It also flips y and z as well as vy and vz to be 
        able to use LogPot of gyrfalcon.
    """

    os.system('mkking' + ' ' +'out=' + output_name + ' ' + 'nbody=' + num_part + ' ' + 'W0=' + w0 + ' ' + 'mass=' + mass + ' ' + 'r_t=' + rt + ' ' + 'WD_units=' + wd_units)
    os.system('snapshift' + ' ' +  output_name + ' ' + output_shifted + 'rshift=' + xs + ',' + ys + ',' + zs + 'vshift=' + vxs + ',' + vys + ',' + vzs)
    os.system('gyrfalcON' + ' ' + 'in=' + output_shifted + ' ' + 'out=' + output_evol + ' ' + 'tstop=' + tstop + ' ' + 'eps=' + eps + ' ' + 'step=' + step + ' ' +'kmax=' + kmax + ' ' + 'Nlev=' + Nlev + ' ' + 'fac=' + fac + ' ' + 'accname=' + accname + ' ' + 'accpars=' + accpars[0] + ',' accpars[1] + ',' + accpars[2] + ',' + accpars[3] + ',' + accpars[4])
    os.system('s2a' + ' ' + output_evol + ' ' +  output_final)




def nemo_read_output(filename):
    
    """
    Parameter:
    -------------------------------------------------------
        filename: filename as a string
        
    Returns:
    -------------------------------------------------------
        mass, position and velocity (the two lather ones
        are in array format)
    """

    data = np.loadtxt(filename)
    mass = data.T[0]
    pos  = data.T[1:4]
    vel  = data.T[4:7]
 
    return mass, pos, vel


def nemo_plot(x,y):

    plt.ion()
    plt.plot(x,y,linewidth=2,color='blue')
    plt.xlabel(str(x),fontsize=20)
    plt.ylabel(str(y),fontsize=20)












