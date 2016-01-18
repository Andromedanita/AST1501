import numpy             as     np
import matplotlib.pylab  as     plt
from   GD1_funcs         import *
from   galpy.actionAngle import actionAngleStaeckel, actionAngleIsochroneApprox
from   galpy.potential   import LogarithmicHaloPotential
from   galpy.util        import bovy_conversion
import os
import copy

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

    data = np.loadtxt(filename, delimiter=',')
    mass = data.T[0]
    pos  = data.T[1:4]
    vel  = data.T[4:7]
 
    return mass, pos, vel


def nemo_coord_convert(pos, vel, q, delta, C_use, ro, vo, m, n):

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
    aAS = actionAngleIsochroneApprox(pot=p, b=0.8)#actionAngleStaeckel(pot = p, delta = delta, c = C_use)
    
    
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
    val = aAS.actionsFreqsAngles(R[m:n],vR[m:n],vT[m:n],zz[m:n],vz[m:n],phi[m:n])
    
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
                    
    aAS = actionAngleIsochroneApprox(pot=p, b=0.8) #actionAngleStaeckel(pot = p, delta = delta, c = C_use)
    val = aAS.actionsFreqsAngles(Rf, vRf, vTf, zzf, vzf, phif)
    return val




def strip_time(filename_tail):
    """
    Parameter:
    -------------------------------------------------------
    
        
        
    Returns:
    -------------------------------------------------------
        Stripping time (eq. 3 in Bovy 2014)
    
    """
    
    data      = np.loadtxt(filename_tail)
    thetar    = data[:,6]
    thetar    = (np.pi+(thetar-np.median(thetar))) % (2.*np.pi)
    indx      = np.fabs(thetar-np.pi) > (5.*np.median(np.fabs(thetar-np.median(thetar))))
    thetar    = thetar[indx]
    thetap    = data[:,7]
    thetap    = (np.pi+(thetap-np.median(thetap))) % (2.*np.pi)
    thetap    = thetap[indx]
    thetaz    = data[:,8]
    thetaz    = (np.pi+(thetaz-np.median(thetaz))) % (2.*np.pi)
    thetaz    = thetaz[indx]
    # center around 0 (instead of pi)
    thetar   -= np.pi
    thetap   -= np.pi
    thetaz   -= np.pi
    # Frequencies
    Or        = data[:,3]
    Op        = data[:,4]
    Oz        = data[:,5]
    dOr       = Or[indx]-np.median(Or)
    dOp       = Op[indx]-np.median(Op)
    dOz       = Oz[indx]-np.median(Oz)
    # Times
    dangle    = np.vstack((thetar,thetap,thetaz))
    dO        = np.vstack((dOr,dOp,dOz))*bovy_conversion.freq_in_Gyr(220.,8.)
    ts        = np.sum(dO*dangle,axis=0)/np.sum(dO**2.,axis=0)
    del_freq  = np.sum(dO**2.,    axis=0)
    del_theta = np.sum(dangle**2.,axis=0)
    return  dO, dangle, del_freq, del_theta, ts



def output_cut(pos, vel, q, delta, C_use, ro, vo, N, var):
    
    m = 0
    
    while n<N:
        if  (N-n) < var:
            val = nemo_coord_convert(pos, vel, q, delta, C_use, ro, vo, m, N)
            np.savetxt("val_tail_{0}.txt".format(N))
        else:
            val = nemo_coord_convert(pos, vel, q, delta, C_use, ro, vo, m, n)
            np.savetxt("val_tail_{0}.txt".format(n))
            m += var
            n += var

    #return


def nemo_plot(x,y,xlabel,ylabel):
    
    plt.ion()
    plt.plot(x,y,linewidth=2,color='blue')
    plt.xlabel(xlabel,fontsize=20)
    plt.ylabel(ylabel,fontsize=20)


def tail_cut(data):
    
    thetar = data[:,6]
    thetar = (np.pi+(thetar-np.median(thetar))) % (2.*np.pi)
    indx   = np.fabs(thetar-np.pi) > (5.*np.median(np.fabs(thetar-np.median(thetar))))
    return indx



def hist_fig4(filename):
    data    = np.loadtxt(filename)
    thetar  = data[:,6]
    thetar  = (np.pi+(thetar-np.median(thetar))) % (2.*np.pi)
    indx    = np.fabs(thetar-np.pi) > (5.*np.median(np.fabs(thetar-np.median(thetar))))
    #Frequencies
    Or  = data[:,3]
    Op  = data[:,4]
    Oz  = data[:,5]
    dOr = Or[indx]-np.median(Or)
    dOp = Op[indx]-np.median(Op)
    dOz = Oz[indx]-np.median(Oz)
    dO  = np.vstack((dOr,dOp,dOz))*bovy_conversion.freq_in_Gyr(220.,8.)
    dO4dir = copy.copy(dO)
    dO4dir[:,dO4dir[:,0] < 0.]*= -1.
    dOdir  = np.median(dO4dir,axis=1)
    dOdir /= np.sqrt(np.sum(dOdir**2.))
    dO1d   = np.dot(dOdir,dO)
    #print "Misalignment:", numpy.arccos(numpy.sum(dOdir*progaa[-1,3:6])/numpy.sqrt(numpy.sum(dOdir**2.)*numpy.sum(progaa[-1,3:6]**2.)))/numpy.pi*180.-180.
    dO1d[dO1d < 0.]*= -1.
    return dO1d



def fig5(filename):
    dO, dangle, del_freq, del_theta, ts = strip_time(filename)
    #Direction in which the stream spreads
    dO4dir = copy.copy(dO)
    dO4dir[:,dO4dir[:,0] < 0.]*= -1.
    dOdir  = np.median(dO4dir,axis=1)
    dOdir /= np.sqrt(np.sum(dOdir**2.))
    #Times
    valx  = np.fabs(np.dot(dangle.T,dOdir))
    valy  = np.fabs(np.dot(dO.T,dOdir))
    return valx, valy





