import numpy            as np
import matplotlib.pylab as plt
from   galpy.util  import bovy_conversion, bovy_coords
import copy



#dir      = "/Users/anita/Documents/Grad_First_Year/Research/data_files/MWpot/"
dir = "/Users/anita/Documents/Grad_First_Year/Research/data_files/LOGpot/"



def cart_spher(x,y,z):
    R         = np.sqrt((x**2.) + (y**2.) + (z**2.))
    theta_sph = np.arctan2(y,x)
    phi_sph   = np.arccos(z/R)

    return np.array([R, theta_sph, phi_sph])

def eig(filename):

    # compute direction in cylindrical
    data    = np.loadtxt(dir+filename)
    thetar  = data[:,6]
    thetar  = (np.pi+(thetar-np.median(thetar))) % (2.*np.pi)
    indx    = np.fabs(thetar-np.pi) > (5.*np.median(np.fabs(thetar-np.median(thetar))))
    
    #Frequencies
    Or      = data[:,3]
    Op      = data[:,4]
    Oz      = data[:,5]
    dOr     = Or[indx]-np.median(Or)
    dOp     = Op[indx]-np.median(Op)
    dOz     = Oz[indx]-np.median(Oz)
    dO      = np.vstack((dOr,dOp,dOz))*bovy_conversion.freq_in_Gyr(220.,8.)
    dO4dir  = copy.copy(dO)
    dO4dir[:,dO4dir[:,0] < 0.]*= -1.
    dOdir   = np.median(dO4dir,axis=1)
    dOdir  /= np.sqrt(np.sum(dOdir**2.))

    # direction in cylindrical
    Or, Ophi, Oz = dOdir[0], dOdir[1], dOdir[2]

    # convert cylindrical to cartesian
    x,y,z = bovy_coords.cyl_to_rect(Or, Ophi, Oz)
    
    # convert cartesian to spherical
    R, theta_sph, phi_sph = cart_spher(x,y,z)

    val = np.array([R, theta_sph, phi_sph])
    return val


#filename = ["concatenated_W2_M2e4_MWpot_orphanway.dat","concatenated_W6_M2e4_MWpot_orphanway.dat", "concatenated_W9_M2e4_MWpot_orphanway.dat","concatenated_W12_M2e4_MWpot_orphanway.dat"]

#filename = ["concatenated_W2.5_M2e4.txt","concatenated_W3_M2e4.txt","concatenated_W6_M2e4.txt","concatenated_W9_M2e4.txt", "concatenated_W12_M2e4.txt"]


filename = ["concatenated_W2_M2e6.txt","concatenated_W3_M2e6.txt","concatenated_W6_M2e6.txt","concatenated_W12_M2e6.txt"]

val_list = []

for names in filename:
    value = eig(names)
    val_list.append(value)

val_array = np.array(val_list)





