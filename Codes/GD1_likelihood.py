from   GD1_funcs    import *
import mw_transform as     mw


# initial positions in cartesian coordinates in natural units
xi,yi,zi = np.array([3.41,13.00,9.58]) # (kpc) in physical units (same as Koposov 2010)
                                        #  with x sign flipped
# initial velocities in cartesian coordinates in natural units
vxi,vyi,vzi = np.array([200.4,-162.6,13.9]) # (km/s) in physical units

# initial coordinates in cylindrical coordinates
Ri,zcyli,phii = xyz_to_cyl(xi,yi,zi) # phi is in radians

# initial velocities in cylindrical coordinates
vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)

ts   = 1000 # number of timesteps
time = np.linspace(0.,1e1,ts)

phi1_pos,phi2_pos,phi2_err = table2_kop2010()
phi1_dist,dist,dist_err    = table3_kop2010()
phi1_vrad,Vrad,V_err       = table1_kop2010()
phi1_mu,mu1,mu2,sigma_mu   = table4_kop2010()


x_err = np.ones(len(phi1))

def contour_singlebox(Vc,q):
    ro = 8.5
    vo = Vc
    p  = potential.LogarithmicHaloPotential(q=q,normalize=1)
    o    = Orbit(vxvv=[Ri/ro,vri/vo,vti/vo,zcyli/ro,vzcyli/vo,phii],ro=ro,vo=vo)
    o.integrate(time,p)
    
    # Orbit in Galactic coordinates
    lval = o.ll(time,ro = ro,obs= [ro,0.,0.])
    bval = o.bb(time,ro = ro,obs= [ro,0.,0.])
    
    
    pmll = o.pmll(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])
    pmbb = o.pmbb(time, ro = ro, obs=[ro,0.,0.,0.,vo,0.])

    galpy_vel  = mw.pmllpmbb_to_pmphi12(pmll,pmbb,lval,bval,degree=True)
    vrad_galpy = o.vlos(time)
    

    phi12 = mw.lb_to_phi12(lval,bval,degree=True)
    phi12[phi12[:,0] > 180,0]-= 360.

    L_pos  = likelihood_all(phi12[:,0],phi1_pos,x_err,phi12[:,1],phi2_pos,phi2_err,time)
    L_dist = likelihood_all(phi12[:,0],phi1_dist,x_err,o.dist(time),dist,dist_err,time)
    L_vrad = likelihood_all(phi12[:,0],phi1_vrad,x_err,vrad_galpy,Vrad,V_err,time)
    L_mu1  = likelihood_all(phi12[:,0],phi1_mu,x_err,galpy_vel.T[0],mu1,sigma_mu,time)
    L_mu2  = likelihood_all(phi12[:,0],phi1_mu,x_err,galpy_vel.T[1],mu2,sigma_mu,time)

    L_total = L_pos * L_dist * L_vrad * L_mu1 * L_mu2
    return L_total



