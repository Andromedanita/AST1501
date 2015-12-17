from   GD1_funcs           import *
import mw_transform        as     mw
import scipy               as     sp
import pickle
#import GD1_likelihood_test as     lt


#----------------------------------------------
# converting initial positions to be in phi1
#            and phi2 coordinate
#----------------------------------------------


# initial position in cartesian coordinate
xi_kop, yi_kop, zi_kop = np.array([3.41,13.00,9.58])

# convert from cartesian to lbd coordinate and returns l and b in degrees
li_kop, bi_kop, di_kop = bovy_coords.XYZ_to_lbd(xi_kop, yi_kop, zi_kop, degree = True)

# convert lb to phi1 and phi2 in degrees
phi12i_kop = mw.lb_to_phi12(li_kop, bi_kop, degree=True)
#phi12i_kop[phi12i_kop[0] > 180,0]-= 360.
    

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

ts        = 1000 # number of timesteps
time_glob = np.linspace(0.,1e1,ts)

def optimizer_func(input,Vc,q):
    
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
    o.integrate(time_glob,p)
    
    time = np.linspace(0.,1e1,1e4)

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
    L_pos  = likelihood_all_test_sum(phi12[:,0], phi1_pos,  x_err_pos,  phi12[:,1],     phi2_pos, phi2_err, time)
    L_dist = likelihood_all_test_sum(phi12[:,0], phi1_dist, x_err_dist, o.dist(time),   dist,     dist_err, time)
    L_vrad = likelihood_all_test_sum(phi12[:,0], phi1_vrad, x_err_vrad, vrad_galpy,     Vrad_kop, V_err,    time)
    L_mu1  = likelihood_all_test_sum(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[0], mu1,      sigma_mu, time)
    L_mu2  = likelihood_all_test_sum(phi12[:,0], phi1_mu,   x_err_mu,   galpy_vel.T[1], mu2,      sigma_mu, time)
    
    print "L position:", L_pos
    print "L dist:"    , L_dist
    print "L vrad:"    , L_vrad
    print "L mu1:"     , L_mu1
    print "L mu2:"     , L_mu2
    print
    
    L_total = L_pos + L_dist + L_vrad + L_mu1 + L_mu2
    
  
    # chi^2 = -2ln(L) where ln(L) = L_total
    chi2 = -2. * L_total


    print
    print "chi2: ",chi2
    print
    
    print "input:", input
    
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
    Nfeval += 1

#print  '{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}'.format('Iter', ' X1', ' X2', ' X3', 'X4', 'X5', 'f(X)')


Vc_inf_list = []
q_inf_list  = []

def optimize(Vc,q):
    
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
    init_guess = (phi12i_kop[1], di_kop, mu_array[0], mu_array[1], Vrad)
    
    val = sp.optimize.minimize(optimizer_func, init_guess, args=(Vc,q), method = 'BFGS', bounds = bnds, options={'maxiter':5,'disp': True,'maxfun':5}, callback=callbackF)
    
    return val.x



# initial guess for proper motion and radial velocity
mu_array, Vrad = vxvyvz_to_pmphi12(xi_kop, yi_kop, zi_kop, vxi_kop, vyi_kop, vzi_kop, True)

# Vc and q arrays
Vc_list = np.linspace(160.,300.,20)
q_list  = np.linspace(0.4,1.6,20)


table   = [[0] * len(Vc_list) for i in range(len(q_list))]
table_contour = np.zeros([len(Vc_list),len(q_list)])


#----------------------------------------------
# computing the optimization parameters for
# each Vc and q in the table and returning a
# 2D list where each cell includes the optimized
# values for phi2, d, mu1, mu2 and Vrad,
# respectively
#----------------------------------------------

for i in range(len(Vc_list)):
    for j in range(len(q_list)):
        print "i = ",i, "j = ", j
        table[j][i] = optimize(Vc_list[i],q_list[j])
        print



file_Name  = "corrected_table_optimization.dat"
fileObject = open(file_Name,'wb')
pickle.dump(table,fileObject)
fileObject.close()

#b = pickle.load(fileObject)



#----------------------------------------------
# computing the log-likelihood values for
# each Vc and q in the table and returning a
# 2D table
#----------------------------------------------
'''
for i in range(len(Vc_list)):
    print i
    for j in range((len(q_list)):
        print j
        table_contour[j][i] = contour_singlebox(Vc_list[i],q_list[j])
  
  '''
                
#----------------------------------------------
#       Plotting log-likelihood contour
#----------------------------------------------
'''
plt.ion()
plt.contourf(Vc_list,q_list,table_contour,100)
plt.xlabel("$V_c \, [kms^{-1}]$",fontsize=20)
plt.ylabel("$q_{\phi}$",fontsize=20)
plt.title("Log-Likelihood")
plt.colorbar()
'''



'''
### test #####

def func(input):
    x = input[0]
    y = input[1]
    z = input[2]
    return x**2 - np.sqrt(y) - z


init_guess = (2,3,4)
bnds       = ((1,None),(0.5,None),(None,None))
sp.optimize.minimize(func, init_guess, method = 'BFGS', bounds = bnds)
'''






