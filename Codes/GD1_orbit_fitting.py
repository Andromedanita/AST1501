from   GD1_funcs    import *
import mw_transform as     mw


# initial positions in cartesian coordinates in natural units
xi,yi,zi = np.array([-3.41,13.00,9.58])/distance  # (kpc) in natural units

# initial velocities in cartesian coordinates in natural units
vxi,vyi,vzi = np.array([-200.4,-162.6,13.9])/220. # (km/s) in natural units


# initial coordinates in cylindrical coordinates
Ri,zcyli,phii = xyz_to_cyl(xi,yi,zi)


# initial velocities in cylindrical coordinates
vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)

phii *= 180./np.pi  # in degrees to use in orbit initialization

# calling the potential
p    = potential.LogarithmicHaloPotential(q=0.9,normalize=1)

# initiating the orbit
ts   = 1000 # number of timesteps
time = np.linspace(0.,1e1,ts)
o    = Orbit(vxvv=[Ri,vri,vti,zcyli,vzcyli,phii],ro=8.5,vo=220.)

o.integrate(time,p)

RA  = o.ra(time,ro=distance,obs=[distance,0.])
DEC = o.dec(time,ro=distance,obs=[distance,0.])

func = np.vectorize(radec_to_phi12)

PHI1,PHI2 = func(RA,DEC,degree=True)
PHI1 *= 180./np.pi
PHI2 *= 180./np.pi

##### RA and Dec to lb --> phi1 and phi2 using Jo's code  #####
lb_test = mw.radec_to_lb(RA,DEC,degree=True,epoch=2000.0)
l_test = lb_test.T[0] # in degree
b_test = lb_test.T[1] # in degree

phi_test = mw.lb_to_phi12(l_test,b_test,degree=True)
phi1_test = phi_test.T[0]
phi2_test = phi_test.T[1]



###### paper data #####
phi1,phi2,phi2_err = table2_kop2010()

plt.ion()
plt.plot(PHI1,PHI2,linewidth=2,color='teal',label='galpy')
plt.plot(phi1,phi2,linewidth=2,color='red',label='GD data')
plt.legend(loc='best')
plt.xlabel("$\phi_1 \, [deg]$",fontsize=20)
plt.ylabel("$\phi_2 \, [deg]$",fontsize=20)
#plt.ylabel("$\ph2_1\,[deg]$",fontsize=20)


'''
############  testing ##################
    
# in R and z coordinates
    
phi1,phi2,phi2_err = table2_kop2010()    # getting Koposov values
dec                = np.zeros(len(phi1)) # initializing array
ra                 = np.zeros(len(phi1)) # ...
    
# converting phi1 and phi2 to RA and DEC
for i in range(len(phi1)):
ra[i],dec[i] = phi12_to_radec(phi1[i],phi2[i],degree=True)
    
# converting RA and DEC to xyz cartesian coordinates
x,y,z       = radec_to_xyz(dec,ra,distance,degree=False)
    
# converting xyz cartesian coordinates to cylindrical coordinates
R,z_cyl,phi = xyz_to_cyl(x,y,z)
    
# in phi1 and phi2 coordinates
Rvals   = o.R(time)
zvals   = o.z(time)
phivals = o.phi(time)
    

plt.figure(1)
plt.plot(Rvals,zvals,linewidth=2,color='blue')
plt.plot(R/distance,z_cyl/distance,'ro')
plt.title("in cylindrical coordinate")
'''


'''
xf,yf,zf = cyl_to_xyz(Rvals,zvals,phivals)
    
func        = np.vectorize(xyz_to_radex)
raf,decf,df = func(xf,yf,zf)
    
phi1f,phi2f = radec_to_phi12(decf,raf,df)
    
plt.figure(2)
plt.plot(phi1f,phi2f,linewidth=2,color='blue')
plt.plot(phi1,phi2,'ro')
plt.title("in stream coordinate")
'''

'''
# testing to compare coordinate transformation from phi12 to cylindrical
# and vice versa to see if they match. They only match if I add 180 degrees to
# the right ascension in phi12_to_radec function :(
    
xf,yf,zf     = cyl_to_xyz(R,z_cyl,phi)
    
#raf,decf,df = xyz_to_radex(xf,yf,zf)
    
func         = np.vectorize(xyz_to_radex)
raf,decf,df  = func(xf,yf,zf)
    
phi1f,phi2f  = radec_to_phi12(raf,decf,df)
'''


