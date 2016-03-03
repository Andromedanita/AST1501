from   GD1_funcs    import *
import mw_transform as     mw
from   utils        import *

#%pylab inline

#----------------------------------------------------------
#                   Koposov 2010 Fig. 1
#----------------------------------------------------------

# reference position and velocity
ro = 8.5
vo = 220.

def orb_fit(params, find_posvel):

    if find_posvel:
        pos, vel    = conversion_likelihood_fig(params)
        xi,yi,zi    = pos[0], pos[1], pos[2]
        vxi,vyi,vzi = vel[0], vel[1], vel[2]
    
    else:
        xi,yi,zi    = np.array([3.41,13.00,9.58])
        vxi,vyi,vzi = np.array([200.4,-162.6,13.9])

    
    Ri,zcyli,phii  = xyz_to_cyl(xi,yi,zi)
    vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)

    p    = potential.LogarithmicHaloPotential(q=0.9,normalize=1)
    ts   = 1000 # number of timesteps
    time = np.linspace(0.,1e1,ts)
    o    = Orbit(vxvv=[Ri/ro,vri/vo,vti/vo,zcyli/ro,vzcyli/vo,phii],ro=ro,vo=vo)
    o.integrate(time,p)

    lval = o.ll(time,ro = ro,obs= [ro,0.,0.])
    bval = o.bb(time,ro = ro,obs= [ro,0.,0.])

    phi12 = mw.lb_to_phi12(lval,bval,degree=True)
    phi12[phi12[:,0] > 180,0]-= 360.

    pmll = o.pmll(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])
    pmbb = o.pmbb(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])

    galpy_vel = mw.pmllpmbb_to_pmphi12(pmll,pmbb,lval,bval,degree=True)

    distvals = o.dist(time)
    Vrad     = o.vlos(time)
    mu1      = galpy_vel.T[0]
    mu2      = galpy_vel.T[1]

    return phi12, distvals, Vrad, mu1, mu2

# max and second max parameters when using 20 bins only
#params_max     = np.array([-3.00012631, 16.40515978, 0.97098013, -2.98959206, -78.66099582])
#params_sec_max = np.array([-1.99474384, 16.5094308, 0.97975578, -2.9905328, -78.59489525])

# max and second max parameters when using 100 bins
params_max        = np.array([ -1.99636521,  16.50468115,   0.97936593,  -2.99043471, -78.60036725])
params_sec_max    = np.array([ -2.13626669,  16.48876609,   0.97821678,  -2.98980611, -78.60974939])
params_third_max  = np.array([ -1.99636521,  16.50468115,   0.97936593,  -2.99043471, -78.60036725])
params_fourth_maz = np.array([ -2.26596714,  16.47531928,   0.97710049,  -2.98962618, -78.61708341])

phi12_max, dist_max,Vrad_max, mu1_max, mu2_max                                     = orb_fit(params_max,        find_posvel = True)
phi12_sec_max, dist_sec_max, Vrad_second_max, mu1_second_max, mu2_second_max       = orb_fit(params_sec_max,    find_posvel = True)
phi12_third_max, dist_third_max, Vrad_third_max, mu1_third_max, mu2_third_max      = orb_fit(params_third_max,  find_posvel = True)
phi12_fourth_max, dist_fourth_max, Vrad_fourth_max, mu1_fourth_max, mu2_fourth_max = orb_fit(params_fourth_maz, find_posvel = True)

phi12_galpy, dist_galpy,vrad_galpy2, mu1_galpy, mu2_galpy = orb_fit(params_max, find_posvel = False)

# Koposov 2010 data
phi1,phi2,phi2_err = table2_kop2010()

x_err = np.ones(len(phi1))


'''
# testing the likelihood
L_list = []
for i in range(len(phi1)):
    l = likelihood_single(phi12[:,0],phi1[i],x_err[i],phi12[:,1],phi2[i],phi2_err[i],time)
    L_list.append(l)

chi2 = -2. * np.log(L_list)
'''

# plotting
plt.ion()
#plt.figure(1,figsize=(15,8))
#plt.subplot(121)
plt.figure(1,figsize=(8,7))
#plt.plot(phi12[:,0],phi12[:,1],linewidth=2,color='black',label='optimization fit')
#plt.plot(phi12[:,0],phi12[:,1],linewidth=2,color='green', marker='o', label='galpy fit')
plt.errorbar(phi1,phi2,xerr = x_err,yerr=phi2_err,marker='o',linestyle='',ecolor='r',color='red',label='GD-1 data')
plt.legend(loc='best')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\phi_2 \, \mathrm{[deg]}$",fontsize=20)
plt.xlim(-80,20)
plt.ylim(-4,2)


plt.plot(phi12_galpy[:,0],phi12_galpy[:,1],linewidth=2,color='green', marker='o', label='galpy fit')
plt.plot(phi12_max[:,0],phi12_max[:,1],linewidth=2,color='black',label='optimization fit 1st max')
plt.plot(phi12_sec_max[:,0],phi12_sec_max[:,1],linewidth=2,color='dodgerblue',label='optimization fit 2nd max')
plt.plot(phi12_third_max[:,0],phi12_third_max[:,1],linewidth=2,color='orange',label='optimization fit 3rd max')
plt.plot(phi12_fourth_max[:,0],phi12_fourth_max[:,1],linewidth=2,color='red',label='optimization fit 4th max')

plt.legend(loc='best')

#plt.subplot(122)
#plt.plot(phi1,chi2,linewidth=2,color='teal',label='likelihood')
#plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
#plt.ylabel("$\chi^2$",fontsize=20)


#----------------------------------------------------------
#                   Koposov 2010 Fig. 2
#----------------------------------------------------------

phi1,mu1,mu2,sigma_mu = table4_kop2010()

#vll  = o.vll(time,ro = ro, obs=[ro,0.,0.,0.,vo,0.])
#vbb  = o.vbb(time,ro = ro, obs=[ro,0.,0.,0.,vo,0.])

#pmll = o.pmll(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])
#pmbb = o.pmbb(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])

#galpy_vel = mw.pmllpmbb_to_pmphi12(pmll,pmbb,lval,bval,degree=True)


#x_err = np.random.ranf(len(phi1))
x_err = np.ones(len(phi1))
'''
L_mu1 = []
L_mu2 = []
for i in range(len(phi1)):
    l1 = likelihood(phi12[:,0],phi1[i],x_err[i],galpy_vel.T[0],mu1[i],sigma_mu[i],time)
    l2 = likelihood(phi12[:,0],phi1[i],x_err[i],galpy_vel.T[1],mu2[i],sigma_mu[i],time)
    L_mu1.append(l1)
    L_mu2.append(l2)

chi2_mu1 = -2. * np.log(L_mu1)
chi2_mu2 = -2. * np.log(L_mu2)
'''

#plt.figure(2,figsize=(15,8))
#plt.subplot(121)
#plt.plot(phi12[:,0],galpy_vel.T[0],linewidth=2,color='black',label='galpy fit $\phi_1$')
#plt.plot(phi12[:,0],galpy_vel.T[1],linewidth=2,color='green',label='galpy fit $\phi_2$')
#plt.errorbar(phi1,mu1,xerr = x_err,yerr=sigma_mu,marker='o',color='red',label='$\mu_{\phi_1}$')
#plt.errorbar(phi1,mu2,xerr = x_err,yerr=sigma_mu,marker='o',color='blue',label='$\mu_{\phi_2}$')
#plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
#plt.ylabel("$\mu \, \mathrm{[mas/yr]}$",fontsize=20)
#plt.xlim(-80,20)
#plt.ylim(-15,-1)
#plt.legend(loc='best')

plt.figure(2)

plt.errorbar(phi1,mu1,xerr = x_err,yerr=sigma_mu,marker='o',linestyle='',ecolor='r',color='red',label=r'GD-1 data $\mu_1$')
plt.errorbar(phi1,mu2,xerr = x_err,yerr=sigma_mu,marker='o',linestyle='',ecolor='purple',color='purple',label=r'GD-1 data $\mu_2$')

plt.plot(phi12_galpy[:,0],
         mu1_galpy,     linewidth=2,color='green',     marker='o', label='galpy fit')
plt.plot(phi12_max[:,0],       mu1_max,       linewidth=2,color='black',     label='optimization fit 4th max')
plt.plot(phi12_sec_max[:,0],   mu1_second_max,linewidth=2,color='dodgerblue',label='optimization fit 1st max')
plt.plot(phi12_third_max[:,0], mu1_third_max, linewidth=2,color='orange',    label='optimization fit 2nd max')
plt.plot(phi12_fourth_max[:,0],mu1_fourth_max,linewidth=2,color='red',       label='optimization fit 3rd max')


plt.plot(phi12_galpy[:,0],     mu2_galpy,     linewidth=2,color='green',     marker='o')
plt.plot(phi12_max[:,0],       mu2_max,       linewidth=2,color='black')
plt.plot(phi12_sec_max[:,0],   mu2_second_max,linewidth=2,color='dodgerblue')
plt.plot(phi12_third_max[:,0], mu2_third_max, linewidth=2,color='orange')
plt.plot(phi12_fourth_max[:,0],mu2_fourth_max,linewidth=2,color='red')

plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\mu \, \mathrm{[mas/yr]}$",fontsize=20)
plt.xlim(-80,20)
plt.ylim(-15,-1)
plt.legend(loc='best')

'''
plt.subplot(122)
plt.plot(phi1,chi2_mu1,linewidth=2,color='teal',label='likelihood $\mu_{\phi_1}$')
plt.plot(phi1,chi2_mu2,linewidth=2,color='red',label='likelihood $\mu_{\phi_2}$')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)
plt.legend(loc='best')
'''
#----------------------------------------------------------
#                   Koposov 2010 Fig. 3
#----------------------------------------------------------

phi1,dist,dist_err = table3_kop2010()

x_err = np.ones(len(phi1))

'''
L_list = []
for i in range(len(phi1)):
    l = likelihood(phi12[:,0],phi1[i],x_err[i],o.dist(time),dist[i],dist_err[i],time)
    L_list.append(l)

chi2_dist = -2. * np.log(L_list)
'''

plt.figure(3)
#plt.figure(3,figsize=(15,8))
#plt.subplot(121)
#plt.plot(phi12[:,0],o.dist(time),linewidth=2,color='black',label='galpy fit')
plt.errorbar(phi1,dist,xerr = x_err,yerr = dist_err,marker='o',linestyle='',color='red',label='GD-1 data')
#plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
#plt.ylabel("distance [kpc]",fontsize=15)
#plt.title("Distance vs. $\phi_1$")
#plt.legend(loc='best')
#plt.xlim(-80,20)
#plt.ylim(6,14)


plt.plot(phi12_galpy[:,0],dist_galpy,linewidth=2,color='green', marker='o', label='galpy fit')
plt.plot(phi12_max[:,0],dist_max,linewidth=2,color='black',label='optimization fit 4th max')
plt.plot(phi12_sec_max[:,0],dist_sec_max,linewidth=2,color='dodgerblue',label='optimization fit 1st max')
plt.plot(phi12_third_max[:,0],dist_third_max,linewidth=2,color='orange',label='optimization fit 2nd max')
plt.plot(phi12_fourth_max[:,0],dist_fourth_max,linewidth=2,color='red',label='optimization fit 3rd max')

plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("distance [kpc]",fontsize=15)
plt.xlim(-80,20)
plt.ylim(6,14)
plt.legend(loc='best')



'''
plt.subplot(122)
plt.plot(phi1,chi2_dist,linewidth=2,color='teal',label='likelihood')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)
'''


#----------------------------------------------------------
#                   Koposov 2010 Fig. 4
#----------------------------------------------------------

phi1,Vrad,V_err = table1_kop2010()
#vrad_galpy      = o.vlos(time)


x_err = np.ones(len(phi1))
'''
L_vrad = []
for i in range(len(phi1)):
    l = likelihood(phi12[:,0],phi1[i],x_err[i],vrad_galpy,Vrad[i],V_err[i],time)
    L_vrad.append(l)

chi2_vrad = -2. * np.log(L_vrad)
'''

'''
plt.figure(4,figsize=(15,8))
plt.subplot(121)
plt.errorbar(phi1,Vrad,xerr = x_err,yerr=V_err,marker='o',linestyle='',color='red',label='GD-1 data')
plt.plot(phi12[:,0],vrad_galpy,linewidth=2,color='black',label='galpy fit')
plt.xlim(-80,20)
plt.ylim(-350,150)
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$V_{rad} \, \mathrm{[km/s]}$",fontsize=20)
plt.legend(loc='best')
'''
'''
plt.subplot(122)
plt.plot(phi1,chi2_vrad,linewidth=2,color='teal',label='likelihood')
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$\chi^2$",fontsize=20)
'''

plt.figure(4)
plt.errorbar(phi1,Vrad,xerr = x_err,yerr=V_err,marker='o',linestyle='',ecolor='r',color='red',label='GD-1 data')

plt.plot(phi12_galpy[:,0],     vrad_galpy2,    linewidth=2,color='green',     marker='o', label='galpy fit')
plt.plot(phi12_max[:,0],       Vrad_max,       linewidth=2,color='black',     label='optimization fit 4th max')
plt.plot(phi12_sec_max[:,0],   Vrad_second_max,linewidth=2,color='dodgerblue',label='optimization fit 1st max')
plt.plot(phi12_third_max[:,0], Vrad_third_max, linewidth=2,color='orange',    label='optimization fit 2nd max')
plt.plot(phi12_fourth_max[:,0],Vrad_fourth_max,linewidth=2,color='red',       label='optimization fit 3rd max')

plt.xlim(-80,20)
plt.ylim(-350,150)
plt.xlabel("$\phi_1 \, \mathrm{[deg]}$",fontsize=20)
plt.ylabel("$V_{rad} \, \mathrm{[km/s]}$",fontsize=20)
plt.legend(loc='best')

