from   GD1_funcs    import *
import mw_transform as     mw
import scipy        as     sp
from galpy.util import bovy_conversion

ro = 8.
vo = 220.
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


ts   = 1000 # number of timesteps
# we need to convert to galpy time units since its not in Gyrs
time = np.linspace(0., nemotime_to_actualtime(5.125)/bovy_conversion.time_in_Gyr(vo,ro), ts)
p    = potential.LogarithmicHaloPotential(q = q, normalize = 1)

# the position and velocity of GD1 stream today in cartesian coordinates
xnow, ynow, znow    = np.array([12.4,1.5,7.1])
vxnow, vynow, vznow = np.array([107.0,-243.0,-105.0])

# the position and velocity of GD1 stream today in cylindrical coordinates
Ri,zcyli,phii  = xyz_to_cyl(xnow,ynow,znow)
vri,vti,vzcyli = vxvyvz_to_vrvtvz(xnow,ynow,znow,-vxnow, -vynow, -vznow)

# initializing the orbit
o = Orbit(vxvv=[Ri/ro, vri/vo, vti/vo, zcyli/ro, vzcyli/vo, phii], ro=ro, vo=vo)
o.integrate(time,p)


print
print "x start is:", o.x(time[-1], ro = ro, obs= [ro,0.,0.])
print "y start is:", o.y(time[-1], ro = ro, obs= [ro,0.,0.])
print "z start is:", o.z(time[-1], ro = ro, obs= [ro,0.,0.])

print "vx start is:", -o.vx(time[ts-1],ro = ro,obs= [ro,0.,0.])
print "vy start is:", -o.vy(time[ts-1],ro = ro,obs= [ro,0.,0.])
print "vz start is:", -o.vz(time[ts-1],ro = ro,obs= [ro,0.,0.])



