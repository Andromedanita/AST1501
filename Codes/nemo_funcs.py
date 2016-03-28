# Anita Bahmanyar
# March 2016, U of T

import numpy             as     np
import matplotlib.pylab  as     plt
from   GD1_funcs         import *
from   galpy.actionAngle import actionAngleStaeckel, actionAngleIsochroneApprox
from   galpy.potential   import LogarithmicHaloPotential
from   galpy.util        import bovy_conversion
from   scipy             import optimize, special
from   galpy.util        import bovy_plot, bovy_coords, bovy_conversion
from   galpy.potential   import MWPotential2014, PowerSphericalPotentialwCutoff, MiyamotoNagaiPotential, NFWPotential
import os
import copy


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

    data = np.loadtxt(filename)#, delimiter=',')
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
    #data      = data.T
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
    """
    Parameter:
    -------------------------------------------------------
        
        
        
    Returns:
        -------------------------------------------------------
        
    """
    
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
    
    """
    Parameter:
    -------------------------------------------------------
        
        
        
    Returns:
    -------------------------------------------------------
        
    """
    
    plt.ion()
    plt.plot(x,y,linewidth=2,color='blue')
    plt.xlabel(xlabel,fontsize=20)
    plt.ylabel(ylabel,fontsize=20)


def tail_cut(data):
    
    """
    Parameter:
    -------------------------------------------------------
        data : concatenated simulation data
        
        
    
    Returns:
    -------------------------------------------------------
        index of the values that are not cut as in the
        tail cut process. This cuts the particles close to
        the stream out so that we get the progenitor at 
        the centre of the tails.
    """
    
    thetar = data[:,6]
    thetar = (np.pi+(thetar-np.median(thetar))) % (2.*np.pi)
    indx   = np.fabs(thetar-np.pi) > (5.*np.median(np.fabs(thetar-np.median(thetar))))
    return indx



def hist_fig4(filename):

    """
    Parameter:
    -------------------------------------------------------
        
        
        
    Returns:
    -------------------------------------------------------
        
    """
    
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
    dO4dir[:,dO4dir[:,0] < 0.] *= -1.
    dOdir  = np.median(dO4dir, axis=1)
    dOdir /= np.sqrt(np.sum(dOdir**2.))
    dO1d   = np.dot(dOdir,dO)
    dO1d[dO1d < 0.] *= -1.
    
    dOpar = np.dot(dOdir,dO)
    
    return dO1d, dOdir #dOpar



def fig5(filename):
    
    """
    Parameter:
    -------------------------------------------------------
        
        
        
    Returns:
    -------------------------------------------------------
    
    """
    
    dO, dangle, del_freq, del_theta, ts = strip_time(filename)
    
    #Direction in which the stream spreads
    dO4dir = copy.copy(dO)
    dO4dir[:,dO4dir[:,0] < 0.] *= -1.
    dOdir  = np.median(dO4dir,axis=1)
    dOdir /= np.sqrt(np.sum(dOdir**2.))
    
    # valx : Delta_theta_parallel
    # valy : Delta_omega_parallel
    valx  = np.fabs(np.dot(dangle.T,dOdir))
    valy  = np.fabs(np.dot(dO.T,dOdir))
    return valx, valy



  
def prob_Oparapar(tdisrupt, i, j, n, ptype):
    """
    Parameter:
    -------------------------------------------------------
        tdisrupt: disruption time in Gyr
        
        i, j : index of Omega and Theta, respectivle
               This is so that it can be used with 
               arrays and for loops easily
        
        n: number of points in omega and theta arrays
        
        ptype: type of the probability of ts
                p(ts) can be either of:
                1) Gaussian
                2) 1/ts
                3) ts
        
        
    Returns:
    -------------------------------------------------------
        Distribution of  parallel, the frequency
        offset between stream members and the progenitor 
        given the angle offset  theta parallel
    """
    
    Omega_par = np.linspace(0.1,0.4,n)
    theta_par = np.linspace(0.01,1.4,n)
    
    # these values are obtained from stream df improvements
    # from galpy. These are such that mu = meam/std = 6

    mean_O   = 0.192340929121
    sigma2_O = 0.0010276398059802503
    print "mu is:", mean_O/np.sqrt(sigma2_O)
  
    ts     = theta_par[j]/Omega_par[i]
    print "ts:", ts
    
    
    if (ts < tdisrupt and ts >= 0.):
        p_omega = np.exp(-0.5 * (Omega_par[i]-mean_O)**2./sigma2_O)/np.sqrt(sigma2_O)
    
    else:
        p_omega = 0.

    
    if ptype == "Gauss":
        p_omegapar_thetapar = (1./Omega_par[i]) * p_omega
    
    elif ptype == "one_over_ts":
        p_omegapar_thetapar = (1./theta_par[j]) * p_omega

    elif ptype == "ts":
        p_omegapar_thetapar = ts * p_omega
    
    return p_omegapar_thetapar




def gausstimesvalue(params,vals,nologsum=False):
    
    """
    Parameter:
    -------------------------------------------------------
        
        
        
    Returns:
    -------------------------------------------------------
    
    """
    
    tmean  = np.exp(params[0])
    tsig   = np.exp(params[1])
    norm   = tsig**2. * np.exp(-tmean**2./2./tsig**2.)+tsig*np.sqrt(np.pi/2.) * tmean * (1.+special.erf(tmean/np.sqrt(2.)/tsig))
    if nologsum:
        return np.fabs(vals)/norm * np.exp(-(vals-tmean)**2./2./tsig**2.)
    else:
        return -np.sum(np.log(np.fabs(vals)/norm * np.exp(-(vals-tmean)**2./2./tsig**2.)))


def plot_gauss(values):
    
    """
    Parameter:
    -------------------------------------------------------
        
        
        
    Returns:
    -------------------------------------------------------
    
    """

    import matplotlib.mlab as mlab
    mean     = np.mean(values)
    variance = np.var(values)
    sigma    = np.sqrt(variance)
    mu_sig   = mean/sigma
    xs       = np.linspace(0.05, np.max(values),200)
    plt.plot(xs, mlab.normpdf(xs,mean,sigma), 'r--', lw=2)
    
    bestfit= optimize.fmin_powell(gausstimesvalue, np.array([np.log(mean*2.), np.log(np.std(values))]), args=(values,))
    bovy_plot.bovy_plot(xs, gausstimesvalue(bestfit, xs, nologsum=True), '-', color='blue', overplot=True, lw=2., zorder=1)
    
    print
    print "Best fit of form output parameters:"
    print "mean is:"    ,np.exp(bestfit[0])
    print "sigma is:"   ,np.exp(bestfit[1])
    print "mu sigma is:", np.exp(bestfit[0])/np.exp(bestfit[1])

    print
    print "Gaussian output parameters:"
    print "mean  is:"   , mean
    print "sigma is:"   , sigma
    print "mu sigma is:", mu_sig
    print


def nemo_pot_params(duration, pot_type, Vo, Ro, q=None):
    
    from galpy.potential import nemo_accname, nemo_accpars
    from calc_shift_nemo import *
    
    '''
    bp = PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
    mp = MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
    np = NFWPotential(a=16/8.,normalize=.35)
    MWPotential2014 = [bp,mp,np]

    accname  = nemo_accname(MWPotential2014)
    accparse = nemo_accpars(MWPotential2014, Vo, Ro)
    '''
    data = calc_init_pos(duration, pot_type, Vo, Ro, q=q)
    return data



def mass_rt_interp(M):
    
    '''
    Parameters
    --------------------------------------------------------------------
        M: mass
        
    Return 
    --------------------------------------------------------------------
        tidal radius(rt) for the given M using the interpolation as well
        as rt(M) obtained from fitting a 3rd order polynomial
    '''
    
    from scipy import interpolate
    
    mass    = np.array([2e4,2e5,2e6,2e7,2e8,2e9])
    rt      = np.array([0.07,0.14,0.32,0.69,1.48,3.20])
    
    z = np.polyfit(mass, rt, 3)
    p = np.poly1d(z)
    
    def fit(mass):
        val = (p[3] * (mass**3)) + (p[2]*(mass**2)) + (p[1]*mass) + p[0]
        return val
    
    func_interp  = sp.interpolate.interp1d(mass, rt)
    
    return func_interp(M), fit(M)



def mass_eps_interp(M):

    '''
    Parameters
    --------------------------------------------------------------------
        M: mass
    
    Return
    --------------------------------------------------------------------
        softening length (epsilon) for the given M using the 
        interpolation as well as epsilon(M) obtained from fitting a 3rd 
        order polynomial
    '''
    
    from scipy import interpolate

    mass    = np.array([2e4,2e5,2e6,2e7,2e8,2e9])
    eps     = np.array([1.5,3.,6.,14.,30.,66.])

    z = np.polyfit(mass, eps, 3)
    p = np.poly1d(z)
    
    def fit(mass):
        val = (p[3] * (mass**3)) + (p[2]*(mass**2)) + (p[1]*mass) + p[0]
        return val

    func_interp  = sp.interpolate.interp1d(mass, eps)

    return func_interp(M), fit(M)



'''
How to run nemo simulation
mkking out=gd1.nemo nbody=1 W0=2. mass=20000 r_t=0.07 WD_units=t
    
snapshift gd1.nemo gd1_shifted.nemo rshift=13.95209126332265193,1.299800690425371164,10.41018639679868407 vshift=-100.6758702857176786,-242.0167773957586235,-17.19903761757237248


gyrfalcON in=gd1_shifted.nemo out=gd1_evol.nemo tstop=5.0 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=PowSphwCut+MiyamotoNagai+NFW accpars=0,1001.79126907,1.8,1.9#0,306770.418682,3.0,0.28#0,16.0,162.958241887

s2a gd1_evol.nemo gd1_evol.dat
'''




