from   GD1_funcs    import *
import mw_transform as     mw
import scipy        as     sp


#----------------------------------------------
# converting initial positions to be in phi1
#            and phi2 coordinate
#----------------------------------------------


# initial position in cartesian coordinate
xi_kop, yi_kop, zi_kop = np.array([3.41,13.00,9.58])

# convert from cartesian to lbd coordinate
li_kop, bi_kop, di_kop = bovy_coords.XYZ_to_lbd(xi_kop, yi_kop, zi_kop, degree = True)

# convert lb to phi1 and phi2
phi12i_kop = mw.lb_to_phi12(li_kop, bi_kop, degree=True)
    

#----------------------------------------------
# converting initial velocities to be in phi1
#            and phi2 coordinate
#----------------------------------------------

# initial velocities in cartesian coordinates in natural units
vxi_kop, vyi_kop, vzi_kop = np.array([200.4,-162.6,13.9]) # (km/s) in physical units
    
# initial velocities in cylindrical coordinates
vri_kop, vti_kop, vzcyli_kop = vxvyvz_to_vrvtvz(xi_kop, yi_kop, zi_kop, vxi_kop, vyi_kop, vzi_kop)


#----------------------------------------------
#  Initializing the Koposov 2010 table values
#----------------------------------------------

phi1_pos,phi2_pos,phi2_err = table2_kop2010()
phi1_dist,dist,dist_err    = table3_kop2010()
phi1_vrad,Vrad_kop,V_err   = table1_kop2010()
phi1_mu,mu1,mu2,sigma_mu   = table4_kop2010()

#----------------------------------------------
#  erros on phi1 for all 4 tables in Koposov
#----------------------------------------------

x_err_pos  = np.ones(len(phi1_pos))
x_err_dist = np.ones(len(phi1_dist))
x_err_vrad = np.ones(len(phi1_vrad))
x_err_mu   = np.ones(len(phi1_mu))

# default values for Vc and q - this can be
# changed to use in the contour plotting
dict = {"Vc":220.,"q":0.9}

ts   = 1000 # number of timesteps
time = np.linspace(0.,1e1,ts)

def optimizer_func(phi2,D,mu_phi1,mu_phi2,Vrad):
    
    '''
    Parameters
    ----------------------------------------------------
        phi2:
            phi2 in degrees
        
        D:
            distance in kpc
        
        
        mu_phi1, mu_phi2:
            proper motion in phi1 and phi2 coordinates
        
        Vrad:
            radial velocity
        
    Functionality
    ----------------------------------------------------
    
        
    Return
    ----------------------------------------------------
        Log likelihood using the specified Vc,q, phi2,
        D, proper motions in stream coordinates and 
        the radial velocity.
    '''
    
    Vc = dict['Vc']
    q  = dict['q']
    
    # choose a value of phi2 for the given phi1 (phi1 should stay constant from
    # initial conditions obtained above in degrees
    
    phi1i = phi12i[0]
    phi2i = phi2
    #phi2i = phi12i[1]
    
    # convert the phi1 and phi2 to be in cylindrical coordinate
    lf, bf        = mw.phi12_to_lb(phi1i, phi2i, degree=True)
    xf, yf, zf    = bovy_coords.lbd_to_XYZ(lf, bf, D, degree = True)
    # guessed initial position in cylindrical coordinate
    Rf,zcylf,phif = xyz_to_cyl(xf, yf, zf)
    
    # convert proper motion in phi1 and phi2 coordinates to
    # proper motion in Galactic (vl,vb) coordinate
    vl,vb = pmphi12_to_pmllpmbb(mu_phi1, mu_phi2, phi1i, phi2i, degree = True)
    
    # convert the vl, vb to proper motion in cartesian coordinate
    vxf, vyf, vzf  = bovy_coords.vrpmllpmbb_to_vxvyvz(Vrad, vl, vb, lf, bf, D, degree = True)
    
    # convert vx,vy,vz to be in cylindrical coordinates
    vrf, vtf, vzff = vxvyvz_to_vrvtvz(xf, yf, zf, vxf, vyf, vzf)
    
    
    # use the above initial positions to calculate the potential and orbit
    ro = 8.5
    vo = Vc
    p  = potential.LogarithmicHaloPotential(q=q,normalize=1)
    o  = Orbit(vxvv=[Rf/ro, vrf/vo, vtf/vo, zcylf/ro, vzff/vo, phif], ro=ro, vo=vo)
    o.integrate(time,p)
    
    #time = np.linspace(0.,1e1,1e4)

    # compute the likelihood with the above initial position
    lval = o.ll(time, ro = ro, obs = [ro,0.,0.])
    bval = o.bb(time, ro = ro, obs = [ro,0.,0.])
    
    pmll = o.pmll(time, ro = ro, obs = [ro, 0., 0., -10., vo+5.25, 7.17])
    pmbb = o.pmbb(time, ro = ro, obs = [ro, 0., 0., -10., vo+5.25, 7.17])
    
    galpy_vel  = mw.pmllpmbb_to_pmphi12(pmll,pmbb,lval,bval,degree = True)
    vrad_galpy = o.vlos(time, ro = ro, obs = [ro, 0., 0., -10., vo+5.25, 7.17])
    
    
    phi12 = mw.lb_to_phi12(lval, bval, degree = True)
    phi12[phi12[:,0] > 180,0]-= 360.
    
    
    L_pos  = likelihood_all_test(phi12[:,0], phi1_pos,  x_err_pos,  phi12[:,1],     phi2_pos, phi2_err, time)
    L_dist = likelihood_all_test(phi12[:,0], phi1_dist, x_err_dist, o.dist(time),   dist,     dist_err, time)
    L_vrad = likelihood_all_test(phi12[:,0], phi1_vrad, x_err_vrad, vrad_galpy,     Vrad_kop, V_err,    time)
    L_mu1  = likelihood_all_test(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[0], mu1,      sigma_mu, time)
    L_mu2  = likelihood_all_test(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[1], mu2,      sigma_mu, time)
    
    # likelihood value in log unit
    L_total = L_pos + L_dist + L_vrad + L_mu1 + L_mu2

    return L_total



def optimize(optimizer_func):




    return



