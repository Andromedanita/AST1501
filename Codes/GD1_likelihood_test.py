from   GD1_funcs    import *
import mw_transform as     mw
import scipy        as     sp


#----------------------------------------------
#  Initializing the Koposov 2010 table values
#----------------------------------------------
phi1_pos,phi2_pos,phi2_err = table2_kop2010()
phi1_dist,dist,dist_err    = table3_kop2010()
phi1_vrad,Vrad,V_err       = table1_kop2010()
phi1_mu,mu1,mu2,sigma_mu   = table4_kop2010()


x_err_pos  = np.ones(len(phi1_pos))
x_err_dist = np.ones(len(phi1_dist))
x_err_vrad = np.ones(len(phi1_vrad))
x_err_mu   = np.ones(len(phi1_mu))



ts   = 1000 # number of timesteps
time = np.linspace(0.,1e1,ts)

def contour_singlebox(phi2, d, mu1, mu2, Vrad, Vc, q):
    ro = 8.5
    vo = Vc
    p  = potential.LogarithmicHaloPotential(q = q, normalize = 1)
    
    
    #---------------------------------------------------------------
    # evaluating the x,y,z and vx,vy,vz from the function variables
    #---------------------------------------------------------------
    phi1          = 34.776756565525091 # degrees
    lf, bf        = mw.phi12_to_lb(phi1, phi2, degree=True)
    xf, yf, zf    = bovy_coords.lbd_to_XYZ(lf, bf, d, degree = True)
    # guessed initial position in cylindrical coordinate
    Rf,zcylf,phif = xyz_to_cyl(xf, yf, zf)
    
    # convert proper motion in phi1 and phi2 coordinates to
    # proper motion in Galactic (vl,vb) coordinate
    vl,vb = mw.pmphi12_to_pmllpmbb(mu1, mu2, phi1, phi2, degree = True)
    
    # convert the vl, vb to proper motion in cartesian coordinate
    vxf, vyf, vzf  = bovy_coords.vrpmllpmbb_to_vxvyvz(Vrad, vl, vb, lf, bf, d, degree = True)
    
    # convert vx,vy,vz to be in cylindrical coordinates
    vrf, vtf, vzff = vxvyvz_to_vrvtvz(xf, yf, zf, vxf, vyf, vzf)
    
    

    o  = Orbit(vxvv=[Rf/ro, vrf/vo, vtf/vo, zcylf/ro, vzff/vo, phif], ro = ro, vo = vo)
    global time
    o.integrate(time,p)
    
    time = np.linspace(0.,1e1,1e4)
    
    # Orbit in Galactic coordinates
    lval = o.ll(time, ro = ro, obs= [ro,0.,0.])
    bval = o.bb(time, ro = ro, obs= [ro,0.,0.])
    
    
    pmll = o.pmll(time, ro = ro, obs=[ro,0.,0.,-10.,vo+5.25,7.17])
    pmbb = o.pmbb(time, ro = ro, obs=[ro,0.,0.,-10.,vo+5.25,7.17])
    
    galpy_vel  = mw.pmllpmbb_to_pmphi12(pmll, pmbb, lval, bval, degree = True)
    vrad_galpy = o.vlos(time, ro = ro, obs = [ro,0.,0.,-10.,vo+5.25,7.17])
    
    
    phi12 = mw.lb_to_phi12(lval, bval, degree = True)
    phi12[phi12[:,0] > 180,0]-= 360.
    
    
    L_pos  = likelihood_all_test_sum(phi12[:,0], phi1_pos,  x_err_pos,  phi12[:,1],     phi2_pos, phi2_err, time)
    L_dist = likelihood_all_test_sum(phi12[:,0], phi1_dist, x_err_dist, o.dist(time),   dist,     dist_err, time)
    L_vrad = likelihood_all_test_sum(phi12[:,0], phi1_vrad, x_err_vrad, vrad_galpy,     Vrad,     V_err,    time)
    L_mu1  = likelihood_all_test_sum(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[0], mu1,      sigma_mu, time)
    L_mu2  = likelihood_all_test_sum(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[1], mu2,      sigma_mu, time)
    
    L_total = L_pos + L_dist + L_vrad + L_mu1 + L_mu2

    return L_total


#Vc_list = np.linspace(160.,300.,20)
#q_list  = np.linspace(0.4,1.6,20)
#Vc_list = np.arange(160,300,15) 
#q_list  = np.arange(0.5,1.5,0.1)
#q_list  = np.arange(0.5,1.5,0.05)
#Vc_list = np.arange(160,280,6)
q_list  = np.arange(0.5,1.5,0.01)
Vc_list = np.arange(151.,300,1.5)

table = np.zeros([len(Vc_list),len(q_list)])

import pickle
#fileObject = open("TNCmethod_table_optimization.dat",'r')
#fileObject = open("TNCmethod_table_optimization_finergrid.dat",'r')
fileObject = open("TNCmethod_table_optimization_100bins_phi1_minus30_phi2_minus_2.dat",'r')
params     = pickle.load(fileObject)

for i in range(len(Vc_list)):
    print i
    for j in range(len(q_list)):
        table[j][i] = contour_singlebox(params[j][i][0], params[j][i][1], params[j][i][2], params[j][i][3], params[j][i][4], Vc_list[i],q_list[j])

'''
Vc_list = np.linspace(160.,300.,20)
q_list  = np.linspace(0.4,1.6,20)


    table = np.zeros([len(Vc_list),len(q_list)])
    
    for i in range(len(Vc_list)):
    print i
    for j in range(len(q_list)):
    table[j][i] = contour_singlebox(Vc_list[i],q_list[j])
    
'''
'''
def P_R0(R0,R0_val,R0_err):
    
    num   = -(R0 - R0_val)**2.
    denom = 2. * (R0_err**2.)
    val   = (1./2./np.pi) * ( np.exp(num/denom) )
    return val


def marginalize_func(R0,Vc,q,R0_val,R0_err):
    
    val = contour_singlebox(Vc,q) * P_R0(R0,R0_val,R0_err)
    return val


def marginalize_int(Vc,q,R0_val,R0_err,R0_min,R0_max):
    
    #int_val,int_err = sp.integrate.quad(marginalize_func,R0_min,R0_max,args=(Vc,q,R0_val,R0_err))
    
    return int_val

'''





