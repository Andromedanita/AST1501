import numpy as np
import matplotlib.pylab as plt

def log_datafit(x, y, deg):
    z = np.polyfit(np.log10(x), np.log10(y), deg)
    p = np.poly1d(z)
    A = np.zeros(np.shape(p)[0])
    for i in range(np.shape(p)[0]):
        A[::-1][i] = p[i]

    yvals = 0.    
    for j in range(np.shape(p)[0]):
        yvals += (((np.log10(x))**j)*A[::-1][j])

    plt.ion()
    plt.loglog(x, y, 'bo', label='Data')
    plt.loglog(x, 10**(yvals), 'g--', lw=2, label='Best Fit')
    plt.legend(loc='best')
    plt.grid(which='minor')
    plt.minorticks_on()

    print "Ax+B"
    print "A = ", A[0]
    print "B =",  A[1]



