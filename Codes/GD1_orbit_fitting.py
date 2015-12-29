from   GD1_funcs    import *
import mw_transform as     mw

#%pylab inline

#----------------------------------------------------------
#                   Koposov 2010 Fig. 1
#----------------------------------------------------------

# reference position and velocity
ro = 8.5
vo = 220.

# initial positions in cartesian coordinates in natural units
xi,yi,zi = np.array([3.41,13.00,9.58]) # (kpc) in physical units (same as Koposov 2010
                                       #  with x sign flipped
# initial velocities in cartesian coordinates in natural units
vxi,vyi,vzi = np.array([200.4,-162.6,13.9]) # (km/s) in physical units

# initial coordinates in cylindrical coordinates
Ri,zcyli,phii = xyz_to_cyl(xi,yi,zi) # phi is in radians

# initial velocities in cylindrical coordinates
vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)


# Initializing potential and orbit
p    = potential.LogarithmicHaloPotential(q=0.9,normalize=1)
ts   = 1000 # number of timesteps
time = np.linspace(0.,1e1,ts)
o    = Orbit(vxvv=[Ri/ro,vri/vo,vti/vo,zcyli/ro,vzcyli/vo,phii],ro=ro,vo=vo)
o.integrate(time,p)


# Orbit in Galactic coordinates
lval = o.ll(time,ro = ro,obs= [ro,0.,0.])
bval = o.bb(time,ro = ro,obs= [ro,0.,0.])

phi12 = mw.lb_to_phi12(lval,bval,degree=True)
phi12[phi12[:,0] > 180,0]-= 360.

# Koposov 2010 data
phi1,phi2,phi2_err = table2_kop2010()

x_err = np.ones(len(phi1))

# testing the likelihood

L_list = []
for i in range(len(phi1)):
    l = likelihood_single(phi12[:,0],phi1[i],x_err[i],phi12[:,1],phi2[i],phi2_err[i],time)
    L_list.append(l)

chi2 = -2. * np.log(L_list)


# plotting
plt.ion()
plt.figure(1,figsize=(15,8))
plt.subplot(121)
plt.plot(phi12[:,0],phi12[:,1],linewidth=2,color='black',marker='o',label='galpy fit')
plt.errorbar(phi1,phi2,xerr = x_err,yerr=phi2_err,marker='o',linestyle='',ecolor='r',color='red',label='GD-1 data')
plt.legend(loc='best')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\phi_2 \, \mathrm{[deg]}$",fontsize=20)
plt.xlim(-80,20)
plt.ylim(-4,2)


plt.subplot(122)
plt.plot(phi1,chi2,linewidth=2,color='teal',label='likelihood')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)
'''
#----------------------------------------------------------
#                   Koposov 2010 Fig. 2
#----------------------------------------------------------

phi1,mu1,mu2,sigma_mu = table4_kop2010()

#vll  = o.vll(time,ro = ro, obs=[ro,0.,0.,0.,vo,0.])
#vbb  = o.vbb(time,ro = ro, obs=[ro,0.,0.,0.,vo,0.])

pmll = o.pmll(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])
pmbb = o.pmbb(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])

galpy_vel = mw.pmllpmbb_to_pmphi12(pmll,pmbb,lval,bval,degree=True)


#x_err = np.random.ranf(len(phi1))
x_err = np.ones(len(phi1))
L_mu1 = []
L_mu2 = []
for i in range(len(phi1)):
    l1 = likelihood(phi12[:,0],phi1[i],x_err[i],galpy_vel.T[0],mu1[i],sigma_mu[i],time)
    l2 = likelihood(phi12[:,0],phi1[i],x_err[i],galpy_vel.T[1],mu2[i],sigma_mu[i],time)
    L_mu1.append(l1)
    L_mu2.append(l2)

chi2_mu1 = -2. * np.log(L_mu1)
chi2_mu2 = -2. * np.log(L_mu2)

plt.figure(2,figsize=(15,8))
plt.subplot(121)
plt.plot(phi12[:,0],galpy_vel.T[0],linewidth=2,color='black',label='galpy fit $\phi_1$')
plt.plot(phi12[:,0],galpy_vel.T[1],linewidth=2,color='green',label='galpy fit $\phi_2$')
plt.errorbar(phi1,mu1,xerr = x_err,yerr=sigma_mu,marker='o',color='red',label='$\mu_{\phi_1}$')
plt.errorbar(phi1,mu2,xerr = x_err,yerr=sigma_mu,marker='o',color='blue',label='$\mu_{\phi_2}$')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\mu \, \mathrm{[mas/yr]}$",fontsize=20)
plt.xlim(-80,20)
plt.ylim(-15,-1)
plt.legend(loc='best')


plt.subplot(122)
plt.plot(phi1,chi2_mu1,linewidth=2,color='teal',label='likelihood $\mu_{\phi_1}$')
plt.plot(phi1,chi2_mu2,linewidth=2,color='red',label='likelihood $\mu_{\phi_2}$')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)
plt.legend(loc='best')

#----------------------------------------------------------
#                   Koposov 2010 Fig. 3
#----------------------------------------------------------

phi1,dist,dist_err = table3_kop2010()

#x_err  = np.random.ranf(len(phi1))
x_err = np.ones(len(phi1))

L_list = []
for i in range(len(phi1)):
    l = likelihood(phi12[:,0],phi1[i],x_err[i],o.dist(time),dist[i],dist_err[i],time)
    L_list.append(l)

chi2_dist = -2. * np.log(L_list)

plt.figure(3,figsize=(15,8))
plt.subplot(121)
plt.plot(phi12[:,0],o.dist(time),linewidth=2,color='black',label='galpy fit')
plt.errorbar(phi1,dist,xerr = x_err,yerr = dist_err,marker='o',linestyle='',color='red',label='GD-1 data')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("distance [kpc]",fontsize=15)
plt.title("Distance vs. $\phi_1$")
plt.legend(loc='best')
plt.xlim(-80,20)
plt.ylim(6,14)


plt.subplot(122)
plt.plot(phi1,chi2_dist,linewidth=2,color='teal',label='likelihood')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)


#----------------------------------------------------------
#                   Koposov 2010 Fig. 4
#----------------------------------------------------------

phi1,Vrad,V_err = table1_kop2010()
vrad_galpy      = o.vlos(time)


x_err = np.ones(len(phi1))
L_vrad = []
for i in range(len(phi1)):
    l = likelihood(phi12[:,0],phi1[i],x_err[i],vrad_galpy,Vrad[i],V_err[i],time)
    L_vrad.append(l)

chi2_vrad = -2. * np.log(L_vrad)


plt.figure(4,figsize=(15,8))
plt.subplot(121)
plt.errorbar(phi1,Vrad,xerr = x_err,yerr=V_err,marker='o',linestyle='',color='red',label='GD-1 data')
plt.plot(phi12[:,0],vrad_galpy,linewidth=2,color='black',label='galpy fit')
plt.xlim(-80,20)
plt.ylim(-350,150)
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$V_{rad} \, \mathrm{[km/s]}$",fontsize=20)
plt.legend(loc='best')


plt.subplot(122)
plt.plot(phi1,chi2_vrad,linewidth=2,color='teal',label='likelihood')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)



'''


'''
#### galpy values ####
RA  = o.ra(time,ro=distance,obs=[distance,0.,0.02])   # in degrees
DEC = o.dec(time,ro=distance,obs=[distance,0.,0.02])  # in degrees

func      = np.vectorize(radec_to_phi12)
PHI1,PHI2 = func(RA,DEC,degree=True) # phi1 and phi2 are in radians

PHI1 *= 180./np.pi # in degrees
PHI2 *= 180./np.pi # in degrees

##### RA and Dec to lb --> phi1 and phi2 using Jo's code  #####
lb_test = bovy_coords.radec_to_lb(RA,DEC,degree=True,epoch=2000.0)
l_test  = lb_test.T[0] # in degree
b_test  = lb_test.T[1] # in degree

phi_test  = mw.lb_to_phi12(l_test,b_test,degree=True)
phi1_test = (phi_test.T[0]) #* 180./np.pi
phi2_test = (phi_test.T[1]) #* 180./np.pi
'''



