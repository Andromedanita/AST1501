import numpy             as     np
import matplotlib.pylab  as     plt
from   GD1_funcs         import *
from   galpy.actionAngle import actionAngleStaeckel
from   galpy.potential   import LogarithmicHaloPotential
import os

#def run_nemo(output_name,num_part, w0, mass, rt, wd_units, output_shifted, xs, ys, zs, vxs, vys, vzs, output_evol, tstop, eps, step, kmax, Nlev, fac, accname, accparse, output_final):

def run_nemo(num_part, w0, mass, rt, wd_units, xs, ys, zs, vxs, vys, vzs, tstop, eps, step, kmax, Nlev, fac, accname, accpars):
    """
    Parameters:
    -------------------------------------------------------
        output_name    : name of the output file
        
        num_part       : number of nbody particles
        
        w0             : dimensionless central potential in 
                         the King profile
        
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
        
        kmax           : the longest timestep is taken to be 2^(-kmax) 
                         simulation time units
        
        Nlev           : 
        
        fac            : These factors control the average time step 
                         of a body to be the minimum
        
        accname        : if given, an external acceleration field wuth 
                        the specified name is used (multiple fields are
                         also accepted in the form of name1 + name2 + name3
        
        accparse       : parameters of accname (should be a list or an array).
                         If multiple fields are given in accname, the params
                         should be given in the form of param1 ; param2 ; param3
        
        output_final   : name of the final output file
        
    
    Returns:
    -------------------------------------------------------
        Runs NEMO for King's profile using the specified 
        parameters. It then shifts the position and velocity 
        to its start time.
        It also flips y and z as well as vy and vz to be 
        able to use LogPot of gyrfalcon.
    """
    """
    os.system('mkking' + ' ' +'out=' + output_name + ' ' + 'nbody=' + num_part + ' ' + 'W0=' + w0 + ' ' + 'mass=' + mass + ' ' + 'r_t=' + rt + ' ' + 'WD_units=' + wd_units)
    print "Done first line"
    os.system('snapshift' + ' ' +  output_name + ' ' + output_shifted + 'rshift=' + xs + ',' + ys + ',' + zs + 'vshift=' + vxs + ',' + vys + ',' + vzs)
    print "Done second line"
    os.system('gyrfalcON' + ' ' + 'in=' + output_shifted + ' ' + 'out=' + output_evol + ' ' + 'tstop=' + tstop + ' ' + 'eps=' + eps + ' ' + 'step=' + step + ' ' +'kmax=' + kmax + ' ' + 'Nlev=' + Nlev + ' ' + 'fac=' + fac + ' ' + 'accname=' + accname + ' ' + 'accpars=' + accpars[0] + ',' + accpars[1] + ',' + accpars[2] + ',' + accpars[3] + ',' + accpars[4])
    print "Done third line"
    os.system('s2a' + ' ' + output_evol + ' ' +  output_final)
    """

    
    os.system('mkking' + ' ' +'out=gd1.nemo'  + ' ' + 'nbody=' + str(num_part) + ' ' + 'W0=' + str(w0) + ' ' + 'mass=' + str(mass) + ' ' + 'r_t=' + str(rt) + ' ' + 'WD_units=' + wd_units)
    print "Done first line"
    os.system('snapshift' + ' ' +  'gd1.nemo' + ' ' + 'gd1_shifted.nemo' + ' ' + 'rshift=' + str(xs) + ',' + str(ys) + ',' + str(zs) + ' ' + 'vshift=' + str(vxs) + ',' + str(vys) + ',' + str(vzs))
    #print 'snapshift' + ' ' +  'gd1.nemo' + ' ' + 'gd1_shifted.nemo' + ' ' + 'rshift=' + str(xs) + ',' + str(ys) + ',' + str(zs) + 'vshift=' + str(vxs) + ',' + str(vys) + ',' + str(vzs)
    print "Done second line"

    print 'gyrfalcON' + ' ' + 'in=gd1_shifted.nemo' + ' ' + 'out=gd1_evol.nemo'  + ' ' + 'tstop=' + str(tstop) + ' ' + 'eps=' + str(eps) + ' ' + 'step=' + str(step) + ' ' +'kmax=' + \
str(kmax) + ' ' + 'Nlev=' + str(Nlev) + ' ' + 'fac=' + str(fac) + ' ' + 'accname=' + str(accname) + ' ' + 'accpars=' + str(accpars[0]) + ',' + str(accpars[1]) + ',' + str(accpars[2]) + ',' + str(accpars[3]) + ','\
              + str(accpars[4])

    os.system('gyrfalcON' + ' ' + 'in=gd1_shifted.nemo' + ' ' + 'out=gd1_evol.nemo'  + ' ' + 'tstop=' + str(tstop) + ' ' + 'eps=' + str(eps) + ' ' + 'step=' + str(step) + ' ' +'kmax=' + \
str(kmax) + ' ' + 'Nlev=' + str(Nlev) + ' ' + 'fac=' + str(fac) + ' ' + 'accname=' + str(accname) + ' ' + 'accpars=' + str(accpars[0]) + ',' + str(accpars[1]) + ',' + str(accpars[2]) + ',' + str(accpars[3]) + ','\
              + str(accpars[4]))
    
    print "Done third line"
    os.system('s2a' + ' ' + 'gd1_evol.nemo' + ' ' +  'gd1_evol.dat')
    


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


def nemo_coord_convert(pos, vel, q, delta, C_use, ro, vo):

    """
    Parameter:
    -------------------------------------------------------
        pos   : x, y, z position as an array
        
        vel   : vx, vy, vz velocities as an array
        
        q     : flattening parameter for Logarithmic potential
        
        delta : focal distance
        
        C_use : True/False - to use C code or not
        
        ro    : radius for conversio to natural units (usually 8 kpc)
        
        vo    : velocity conversio to natural units(usually 220 km/s)
        
    Returns:
    -------------------------------------------------------
        action-angle coordinates and the omega (frequency) as an array
        array([jr,lz,jz,Omegar,Omegaphi,Omegaz,angler,anglephi,anglez])
        This is for the tail.
    
    """
    
    # Logarithmic potential and action-angle function initiated
    p   = LogarithmicHaloPotential(q = q, normalize = 1.)
    aAS = actionAngleStaeckel(pot = p, delta = delta, c = C_use)
    
    
    # position and velocity in cartesian coordinates
    # This flips y and z and vy and vz columns
    x, y, z    = pos[0], pos[2], pos[1]
    vx, vy, vz = vel[0], vel[2], vel[1]
    
    # make these functions vectorizable to use with arrays
    xyz_to_cyl_vect       = np.vectorize(xyz_to_cyl)
    vxvyvz_to_vrvtvz_vect = np.vectorize(vxvyvz_to_vrvtvz)
    
    # position and velocity in cylindrical coordinates
    R, zz , phi = xyz_to_cyl_vect(x, y, z)
    vR, vT, vz  = vxvyvz_to_vrvtvz_vect(x, y, z, vx, vy, vz)
    
    ro = ro
    vo = vo
    
    # convert to natural units for use in galpy
    R  /= ro
    zz /= ro
    vR /= vo
    vT /= vo
    vz /= vo
    
    # action-angle and omega values
    val = aAS.actionsFreqsAngles(R,vR,vT,zz,vz,phi)
    
    return val



def nemo_prog_action_angle(x, y, z, vx, vy, vz, R0, V0, q, end_time, delta, C_use):
    """
    Parameter:
    -------------------------------------------------------
        x, y, z : initial position of the progenitor
        
        vx, vy, vz : initial velocity of the progenitor
        
        R0, V0  : radius and circular velocity
        
        q       : flattening parameter
        
        end_time : time when the simulation ended
        
        delta    : focal distance
        
        C_use    : True/False - to use C code or not
        
       
    Returns:
    -------------------------------------------------------
        action-angle coordinates and the omega (frequency) as an array
        array([jr,lz,jz,Omegar,Omegaphi,Omegaz,angler,anglephi,anglez])
        This is for the progenitor.

    """
    
    p   = LogarithmicHaloPotential(q = q, normalize = 1.)
    xyz_to_cyl_vect       = np.vectorize(xyz_to_cyl)
    vxvyvz_to_vrvtvz_vect = np.vectorize(vxvyvz_to_vrvtvz)
    
    R, zz , phi = xyz_to_cyl_vect(x, y, z)
    vR, vT, vz  = vxvyvz_to_vrvtvz_vect(x, y, z, vx, vy, vz)
    
    
    # convert to natural units for use in galpy
    R  /= R0
    zz /= R0
    vR /= V0
    vT /= V0
    vz /= V0
    
    # initializing the orbit
    o = Orbit(vxvv=[R, vR, vT, zz, vz, phi], ro=R0, vo=V0)
    
    # to convert the time units to normal
    t = end_time * (V0/R0)
    time = np.linspace(0., t, 1e4)
    o.integrate(time, p)
                    
    Rf   = o.R(time)
    zzf  = o.z(time)
    vRf  = o.vR(time)
    vTf  = o.vT(time)
    vzf  = o.vz(time)
    phif = o.phi(time)
                    
    aAS = actionAngleStaeckel(pot = p, delta = delta, c = C_use)
    val = aAS.actionsFreqsAngles(Rf, vRf, vTf, zzf, vzf, phif)
    #o.plot()
    return val




def strip_time(filename, x, y, z, vx, vy, vz, R0, V0, q, end_time, delta, C_use):
    """
    Parameter:
    -------------------------------------------------------
    
        
        
    Returns:
    -------------------------------------------------------
        Stripping time (eq. 3 in Bovy 2014)
    
    """
    # properties of the tail
    mass, pos, vel = nemo_read_output(filename)
    
    val_prog = nemo_prog_action_angle(x, y, z, vx, vy, vz, R0, V0, q, end_time, delta, C_use)
    val_tail = nemo_coord_convert(pos, vel, q, delta, C_use, R0, V0)

    # progenitor frequency
    freq_r_prog   = val_prog[3]
    freq_phi_prog = val_prog[4]
    freq_z_prog   = val_prog[5]
    
    # tail frequency
    freq_r_tail   = val_tail[3]
    freq_phi_tail = val_tail[4]
    freq_z_tail   = val_tail[5]
    
    # progenitor angle
    theta_r_prog   = val_prog[6]
    theta_phi_prog = val_prog[7]
    theta_z_prog   = val_prog[8]
    
    # tail angle
    theta_r_tail   = val_tail[6]
    theta_phi_tail = val_tail[7]
    theta_z_tail   = val_tail[8]
    
    # frequency offset
    del_freq  = np.array([freq_r_tail-freq_r_prog, freq_phi_tail-freq_phi_prog, freq_z_tail-freq_z_prog])
    # angle offset
    del_theta = np.array([theta_r_tail-theta_r_prog, theta_phi_tail-theta_phi_prog, theta_z_tail-theta_z_prog])
    
    num   = np.dot(del_freq, del_theta)
    #denom =

    #ts = num/denom
    #return ts



def nemo_plot(x,y,xlabel,ylabel,label):
    
    plt.ion()
    plt.plot(x,y,linewidth=2,color='blue', label = label)
    plt.xlabel(xlabel,fontsize=20)
    plt.ylabel(ylabel,fontsize=20)




