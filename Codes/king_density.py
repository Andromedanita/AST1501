import numpy            as     np
import matplotlib.pylab as     plt
from   scipy            import integrate, special, interpolate

rho1 = 1.
G    = 1. #6.673e-11
sig  = 1.
sig2 = sig**2.
W0   = 3.

rrange = np.arange(0.01,100.,0.01)

def solvr(psi, r):

    x1   = psi[0]
    x2   = psi[1]
    val  = (np.exp(x1/sig2) * special.erf(np.sqrt(x1)/sig)) - (np.sqrt((4. * x1)/(np.pi*sig2)) * (1. + ((2*x1)/(3*sig2))))
    val_fin = (-4 * np.pi * G * rho1 * val) - ((2. * x2)/r)
    
    print
    print "x1:   ", x1
    print "x2:   ", x2
    print "exp:  ", np.exp(x1/sig2)
    print "efr:  ", special.erf(np.sqrt(x1)/sig)
    print "sqrt: ", np.sqrt((4. * x1)/(np.pi*sig2))
    print "last: ", 1. + ((2*x1)/(3*sig2))
    print "2/r:", (2. * x2)/r
    print

    return np.array([x2, val_fin], float)



def dens():
    '''
    Return
    -------------------------------------------
        solution of the ode by solving the 
        integral numerically. 
        rsol is an array with 2 columns
        rsol[:,0] = psi(r) values
        rsol[:,1] = psi'(r) (derivative)
    
    '''
    
    rsol = integrate.odeint(solvr, [W0 * sig2,0.01], rrange)
    w    = np.column_stack((rrange, rsol[:,0], rsol[:,1]))
    
    if np.isnan(rsol[:,0]):
        np.pop(rsol)
    np.savetxt("test_sols.txt", w, fmt="%.9e")
    return rsol


#if __name__ == '__main__':
#    q = dens()

# Manual 4th-order Runge-Kutta
a = 0.
b = 100.
N = 1000
h = (b-a)/N

rpoints  = np.arange(a,b,h)
x1points = []
x2points = []

t = np.array([W0*sig2, 0.],float)
for r in rpoints:
    x1points.append(t[0])
    x2points.append(t[1])
    k1 = h*solvr(t,r)
    k2 = h*solvr(t+0.5*k1, r+0.5*h)
    k3 = h*solvr(t+0.5*k2, r+0.5*h)
    k4 = h*solvr(t+k3, r+h)
    t += (k1+2*k2+2*k3+k4)/6

plt.ion()
plt.plot(rpoints,x1points)
plt.plot(rpoints,x2points)



def psi_r(r):
    '''
    Parameter
    -------------------------------------------
        r: radius
    
    Return
    -------------------------------------------
        interpolated value of psi(r) obtained
        from solving the ode numerically.
    '''
    val = interpolate.interp1d(q[:,0],rrange)
    return val


def rho_r(r):
    x1 = psi_r(r)
    val = rho1 * ((np.exp(x1/sig2) * special.erf(np.sqrt(x1)/sig)) - (np.sqrt((4. * x1)/(np.pi*sig2)) * (1. + ((2*x1)/(3*sig2)))))
    return val


'''
    data = np.loadtxt("test_sols.txt")
    plt.ion()
    plt.plot(data[:,0], data[:,1], 'o')
    plt.xlabel("$r$", fontsize=20)
    plt.ylabel("$\psi(r)$", fontsize=20)
    '''
