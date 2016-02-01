from   GD1_funcs       import *
import mw_transform    as     mw
import scipy           as     sp
from   galpy.util      import bovy_conversion
from   galpy.potential import MWPotential2014


q  = 0.9

def actualtime_to_nemotime(t):
    '''
    t should be given in Gyrs.
    Converts actual time to nemo time units
    '''
    t     *= 1000.             # to convert to Myrs
    unit   = 977.7922212082034 # 1 nemo time unit is this many Myrs
    val    = t/unit
    return val


def nemotime_to_actualtime(tunit):
    '''
    Converts nemo time units to actual time in Gyr
    '''
    
    unit   = 977.7922212082034 # 1 nemo time unit is this many Myrs
    val    = unit * tunit
    val   /= 1000.             # to convert to Gyrs
    return val


def calc_init_pos(duration, pot_type, Vo, Ro, q=None):

    ro = Ro
    vo = Vo
    ts   = 1000 # number of timesteps
    # we need to convert to galpy time units since its not in Gyrs
    time = np.linspace(0., nemotime_to_actualtime(duration)/bovy_conversion.time_in_Gyr(vo,ro), ts)
    
    if pot_type == "Log":
        p    = potential.LogarithmicHaloPotential(q = q, normalize = 1)
    
    elif pot_type == "MW2014":
        from   galpy.potential   import MWPotential2014, PowerSphericalPotentialwCutoff, MiyamotoNagaiPotential, NFWPotential
        bp  = PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
        mp  = MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
        npp = NFWPotential(a=16/8.,normalize=.35)
        p   = MWPotential2014 = [bp,mp,npp]

    # the position and velocity of GD1 stream today in cartesian coordinates
    xnow, ynow, znow    = np.array([12.4,1.5,7.1])
    vxnow, vynow, vznow = np.array([107.0,-243.0,-105.0])

    # the position and velocity of GD1 stream today in cylindrical coordinates
    Ri,zcyli,phii  = xyz_to_cyl(xnow,ynow,znow)
    vri,vti,vzcyli = vxvyvz_to_vrvtvz(xnow,ynow,znow,-vxnow, -vynow, -vznow)

    # initializing the orbit
    o = Orbit(vxvv=[Ri/ro, vri/vo, vti/vo, zcyli/ro, vzcyli/vo, phii], ro=ro, vo=vo)
    o.integrate(time,p)

    x_init = o.x(time[-1], ro = ro, obs= [ro,0.,0.])
    y_init = o.y(time[-1], ro = ro, obs= [ro,0.,0.])
    z_init = o.z(time[-1], ro = ro, obs= [ro,0.,0.])
    
    vx_init = -o.vx(time[ts-1],ro = ro,obs= [ro,0.,0.])
    vy_init = -o.vy(time[ts-1],ro = ro,obs= [ro,0.,0.])
    vz_init = -o.vz(time[ts-1],ro = ro,obs= [ro,0.,0.])

    return np.array([x_init, y_init, z_init, vx_init, vy_init, vz_init])


