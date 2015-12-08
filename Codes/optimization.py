from   GD1_funcs           import *
import mw_transform        as     mw
import scipy               as     sp
import GD1_likelihood_test as     lt


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

def optimizer_func(input):#(phi2,D,mu_phi1,mu_phi2,Vrad):
    
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
    
        
    Return
    ----------------------------------------------------
        Log likelihood using the specified Vc,q, phi2,
        D, proper motions in stream coordinates and 
        the radial velocity.
    '''
    phi2    = input[0]
    D       = input[1]
    mu_phi1 = input[2]
    mu_phi2 = input[3]
    Vrad    = input[4]
    
    Vc = dict['Vc']
    q  = dict['q']
    
    print "vc in function:", Vc
    print "q in function:" , q
    
    # choose a value of phi2 for the given phi1 (phi1 should stay constant from
    # initial conditions obtained above in degrees
    
    phi1i = phi12i_kop[0]
    phi2i = phi2
    
    # convert the phi1 and phi2 to be in cylindrical coordinate
    lf, bf        = mw.phi12_to_lb(phi1i, phi2i, degree=True)
    xf, yf, zf    = bovy_coords.lbd_to_XYZ(lf, bf, D, degree = True)
    # guessed initial position in cylindrical coordinate
    Rf,zcylf,phif = xyz_to_cyl(xf, yf, zf)
    
    # convert proper motion in phi1 and phi2 coordinates to
    # proper motion in Galactic (vl,vb) coordinate
    vl,vb = mw.pmphi12_to_pmllpmbb(mu_phi1, mu_phi2, phi1i, phi2i, degree = True)
    
    # convert the vl, vb to proper motion in cartesian coordinate
    vxf, vyf, vzf  = bovy_coords.vrpmllpmbb_to_vxvyvz(Vrad, vl, vb, lf, bf, D, degree = True)
    
    # convert vx,vy,vz to be in cylindrical coordinates
    vrf, vtf, vzff = vxvyvz_to_vrvtvz(xf, yf, zf, vxf, vyf, vzf)
    
    
    # use the above initial positions to calculate the potential and orbit
    ro = 8.5
    vo = Vc
    p  = potential.LogarithmicHaloPotential(q = q,normalize = 1)
    o  = Orbit(vxvv=[Rf/ro, vrf/vo, vtf/vo, zcylf/ro, vzff/vo, phif], ro = ro, vo = vo)
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
    
    # calculating likelihood for each set of parameters
    L_pos  = likelihood_all_test(phi12[:,0], phi1_pos,  x_err_pos,  phi12[:,1],     phi2_pos, phi2_err, time)
    L_dist = likelihood_all_test(phi12[:,0], phi1_dist, x_err_dist, o.dist(time),   dist,     dist_err, time)
    L_vrad = likelihood_all_test(phi12[:,0], phi1_vrad, x_err_vrad, vrad_galpy,     Vrad_kop, V_err,    time)
    L_mu1  = likelihood_all_test(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[0], mu1,      sigma_mu, time)
    L_mu2  = likelihood_all_test(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[1], mu2,      sigma_mu, time)
    
    # likelihood value in log unit
    L_total = L_pos + L_dist + L_vrad + L_mu1 + L_mu2
    
    # chi^2 = -2ln(L) where ln(L) = L_total
    chi2 = -2. * L_total

    return chi2




def check_phi12_to_cylind(phi1,phi2, d, degree = False):
    
    '''
    Parameters
    ----------------------------------------------------
        phi1 and phi2:
            position in stream coordinates
        
        d:
            distance in kpc
        
        
    Return
    ----------------------------------------------------
        position in cylindrical coordinates
    '''
    
    l, b          = mw.phi12_to_lb(phi1, phi2, degree)
    x, y, z       = bovy_coords.lbd_to_XYZ(l, b, d, degree)
    R, zcyl , phi = xyz_to_cyl(xf, yf, zf)
    return np.array([R, zcyl, phi])



def vxvyvz_to_pmphi12(x, y, z, vx, vy, vz, degree):
    
    '''
    Parameters
    ----------------------------------------------------
        x, y, z : 
            position in cartesian coordinate
            
        vx, vy, vz :
            velocity in cartesian coordinate
        
    Return
    ----------------------------------------------------
        proper motion in stream coordinate as an array
        and the radial velocity as a number

    '''

    l, b, d      = bovy_coords.XYZ_to_lbd(x, y, z, degree = degree)
    vr, vl, vb   = bovy_coords.vxvyvz_to_vrpmllpmbb(vx,vy,vz, l, b, d, degree = degree)
    vphi1, vphi2 = mw.pmllpmbb_to_pmphi12(vl, vb, l, b, degree = degree)

    return np.array([vphi1,vphi2]), vr

Nfeval = 1

def callbackF(Xi):
    global Nfeval
    print '{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}'.format(Nfeval, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], optimizer_func(Xi))
    Nfeval += 1

print  '{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}'.format('Iter', ' X1', ' X2', ' X3', 'X4', 'X5', 'f(X)')


def optimize():
    
    '''
    Parameters
    ----------------------------------------------------
        
        
    Return
    ----------------------------------------------------
        returns an array with shape (5,1) that includes
        the optimized parameters from optimizer_func()
        function. The parameters are phi2, D, mu_phi1,
        mu_phi2 and Vrad, respectively.

    '''
    
    bnds       = ((-90., 90.), (0., None), (None,None), (None,None), (None,None))
    init_guess = (phi12i_kop[1], 8.5, mu_array[0], mu_array[1], Vrad)
    val       = sp.optimize.minimize(optimizer_func, init_guess, method = 'BFGS', bounds = bnds, callback=callbackF)
    #val       = sp.optimize.fmin_l_bfgs_b(optimizer_func, init_guess, bounds = bnds)
    return val.x


# initial guess for proper motion and radial velocity
mu_array, Vrad = vxvyvz_to_pmphi12(xi_kop, yi_kop, zi_kop, vxi_kop, vyi_kop, vzi_kop, True)

# Vc and q arrays
Vc_list = np.linspace(160.,300.,20)
q_list  = np.linspace(0.4,1.6,20)

table   = [[0] * len(Vc_list) for i in range(1)]#range(len(q_list))]
table_contour = np.zeros([len(Vc_list),len(q_list)])



#----------------------------------------------
# computing the optimization parameters for
# each Vc and q in the table and returning a
# 2D list where each cell includes the optimized
# values for phi2, d, mu1, mu2 and Vrad,
# respectively
#----------------------------------------------

for i in range(len(Vc_list)):
    dict['Vc'] = Vc_list[i]
    print
    print "Vc is", dict['Vc']
    for j in range(1):#(len(q_list)):
        dict['q'] = q_list[j]
        print "going through", j
        print "q is:", dict['q']
        table[j][i] = optimize()
        print "done one q"
        print


#----------------------------------------------
# computing the log-likelihood values for
# each Vc and q in the table and returning a
# 2D table
#----------------------------------------------

for i in range(len(Vc_list)):
    print i
    for j in range((len(q_list)):
        print j
        table_contour[j][i] = contour_singlebox(Vc_list[i],q_list[j])
  
                   
                
#----------------------------------------------
#       Plotting log-likelihood contour
#----------------------------------------------
                   
plt.ion()
plt.contourf(Vc_list,q_list,table_contour,100)
plt.xlabel("$V_c \, [kms^{-1}]$",fontsize=20)
plt.ylabel("$q_{\phi}$",fontsize=20)
plt.title("Log-Likelihood")
plt.colorbar()




'''
### test #####

def func(x,y,z):
    return x+y+z


init_guess = (2,3,4)
bnds       = ((1,None),(0.5,None),(None,None))
sp.optimize.minimize(func, init_guess, method = 'BFGS', bounds = bnds)
'''






