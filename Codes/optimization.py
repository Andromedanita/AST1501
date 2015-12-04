from   GD1_funcs    import *
import mw_transform as     mw
import scipy        as     sp

# converting initial condition to be in phi1 and phi2 coordinate

### initial position in cartesian coordinate
xi,yi,zi = np.array([3.41,13.00,9.58])

### convert from cartesian to lbd coordinate
li,bi,di = bovy_coords.XYZ_to_lbd(xi,yi,zi,degree = True)

### convert lb to phi1 and phi2
phi12i = mw.lb_to_phi12(li,bi,degree=True)
    
    
# initial velocities in cartesian coordinates in natural units
vxi,vyi,vzi = np.array([200.4,-162.6,13.9]) # (km/s) in physical units
    
# initial velocities in cylindrical coordinates
vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)

ts   = 1000 # number of timesteps
time = np.linspace(0.,1e1,ts)

phi1_pos,phi2_pos,phi2_err = table2_kop2010()
phi1_dist,dist,dist_err    = table3_kop2010()
phi1_vrad,Vrad,V_err       = table1_kop2010()
phi1_mu,mu1,mu2,sigma_mu   = table4_kop2010()


x_err_pos  = np.ones(len(phi1_pos))
x_err_dist = np.ones(len(phi1_dist))
x_err_vrad = np.ones(len(phi1_vrad))
x_err_mu   = np.ones(len(phi1_mu))


dict = {"Vc":220.,"q":0.9}


def optimize(phi2,D,mu_phi1,mu_phi2,Vrad):
    
    Vc = dict['Vc']
    q  = dict['q']
    
    # choose a value of phi2 for the given phi1 (phi1 should stay constant from
    # initial conditions obtained above
    
    phi1i = phi12i[0]
    phi2i = phi2
    #phi2i = phi12i[1]
    
    # convert the phi1 and phi2 to be in cylindrical coordinate
    lf, bf        = mw.phi12_to_lb(phi1i,phi2i,degree=True)
    xf, yf, zf    = bovy_coords.lbd_to_XYZ(lf,bf,D, degree = True)
    Rf,zcylf,phif = xyz_to_cyl(xf,yf,zf)


    # initial guess for the velocities
    # need to convert from cylindrical velocities to muphi1 and muphi2
    # to get an idea of the initial guess and then convert it to
    # be in cylindrical coordinate again to be used in intializing the
    # orbit of the stream

    # use the above initial positions to calculate the potential and orbit
    ro = 8.5
    vo = Vc
    p  = potential.LogarithmicHaloPotential(q=q,normalize=1)
    o  = Orbit(vxvv=[Rf/ro,vri/vo,vti/vo,zcylf/ro,vzcyli/vo,phif],ro=ro,vo=vo)
    o.integrate(time,p)
    
    time = np.linspace(0.,1e1,1e4)

    # compute the likelihood with the above initial position
    lval = o.ll(time,ro = ro,obs= [ro,0.,0.])
    bval = o.bb(time,ro = ro,obs= [ro,0.,0.])
    
    pmll = o.pmll(time, ro = ro, obs=[ro,0.,0.,-10.,vo+5.25,7.17])
    pmbb = o.pmbb(time, ro = ro, obs=[ro,0.,0.,-10.,vo+5.25,7.17])
    
    galpy_vel  = mw.pmllpmbb_to_pmphi12(pmll,pmbb,lval,bval,degree=True)
    vrad_galpy = o.vlos(time,ro = ro, obs=[ro,0.,0.,-10.,vo+5.25,7.17])
    
    
    phi12 = mw.lb_to_phi12(lval,bval,degree=True)
    phi12[phi12[:,0] > 180,0]-= 360.
    
    
    L_pos  = likelihood_all_test(phi12[:,0],phi1_pos,x_err_pos,phi12[:,1],phi2_pos,phi2_err,time)
    L_dist = likelihood_all_test(phi12[:,0],phi1_dist,x_err_dist,o.dist(time),dist,dist_err,time)
    L_vrad = likelihood_all_test(phi12[:,0],phi1_vrad,x_err_vrad,vrad_galpy,Vrad,V_err,time)
    L_mu1  = likelihood_all_test(phi12[:,0],phi1_mu,x_err_mu,galpy_vel.T[0],mu1,sigma_mu,time)
    L_mu2  = likelihood_all_test(phi12[:,0],phi1_mu,x_err_mu,galpy_vel.T[1],mu2,sigma_mu,time)
    
    # likelihood value in log unit
    L_total = L_pos + L_dist + L_vrad + L_mu1 + L_mu2

    
    # optimize the initial paramters (maximize the ln(likeliood) = L_total)
    # return the computed initial position to use it for marginalization


    #return



