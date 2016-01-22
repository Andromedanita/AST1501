from   GD1_funcs    import *
import mw_transform as     mw
import scipy        as     sp


ts   = 1000 # number of timesteps
time = np.linspace(0.,5.,ts)


ro = 8.
vo = 220.
q  = 0.9
p  = potential.LogarithmicHaloPotential(q = q, normalize = 1)

xnow, ynow, znow = np.array([-11.63337239, -20.76235661, -10.631736273934635])
vxnow, vynow, vznow = np.array([-128.8281653, 42.88727925, 79.172383882274971])
    
#Ri,zcyli,phii = xyz_to_cyl(xi,yi,zi)
#vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)
    
#o    = Orbit(vxvv=[Ri/ro,vri/vo,vti/vo,zcyli/ro,vzcyli/vo,phii],ro=ro,vo=vo)
#o.integrate(time,p)

def optimizer_func(input):

    xi  = input[0]
    yi  = input[1]
    zi  = input[2]
    vxi = input[3]
    vyi = input[4]
    vzi = input[5]

    Ri,zcyli,phii = xyz_to_cyl(xi,yi,zi)
    vri,vti,vzcyli = vxvyvz_to_vrvtvz(xi,yi,zi,vxi,vyi,vzi)
    
    o = Orbit(vxvv=[Ri/ro,vri/vo,vti/vo,zcyli/ro,vzcyli/vo,phii],ro=ro,vo=vo)
    o.integrate(time,p)

    var1 = np.fabs(np.fabs(o.x(time[999],ro = ro,obs= [ro,0.,0.]))- 12.4)
    var2 = np.fabs(np.fabs(o.y(time[999],ro = ro,obs= [ro,0.,0.]))- 1.5)
    var3 = np.fabs(np.fabs(o.z(time[999],ro = ro,obs= [ro,0.,0.]))- 7.1)

    var4 = np.fabs(np.fabs(o.x(time[999],ro = ro,obs= [ro,0.,0.]))- 128.8281653)
    var5 = np.fabs(np.fabs(o.x(time[999],ro = ro,obs= [ro,0.,0.]))- 42.88727925)
    var6 = np.fabs(np.fabs(o.x(time[999],ro = ro,obs= [ro,0.,0.]))- 79.172383882274971)
    
    fin_val = var1 + var2 + var3 + var4 + var5 + var6
    
    print
    print "x diff :",  var1
    print "y diff :",  var2
    print "z diff :",  var3
    print "vx diff :", var4
    print "vy diff :", var5
    print "vz diff :", var6
    print "final value:",fin_val
    print
    
    return fin_val



Nfeval = 1

def callbackF(Xi):
    global Nfeval
    Nfeval += 1


def optimize():
    
    bnds       = ((-15.,-9), (-30., -10.), (-20.,0.), (-140.,-115.), (35.,50.), (65.,90.))
    init_guess = (-11., -20., -10., -128., 42., 79.)
    
    val = sp.optimize.minimize(optimizer_func, init_guess, method = 'L-BFGS-B', bounds = bnds, options={'maxiter':5,'disp': True,'maxfun':5}, callback=callbackF)
    
    return val.x




table = optimize()


