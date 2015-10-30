def get_epoch_angles(epoch=2000.0):
    """
    NAME:
       get_epoch_angles
    PURPOSE:
       get the angles relevant for the transformation from ra, dec to l,b for the given epoch
    INPUT:
       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported
    OUTPUT:
       set of angles
    HISTORY:
       2010-04-07 - Written - Bovy (NYU)
    """
    import math as m
    import numpy as nu
    import scipy as sc
    
    if epoch == 2000.0:
        theta= 122.932/180.*sc.pi
        dec_ngp= 27.12825/180.*sc.pi
        ra_ngp= 192.85948/180.*sc.pi
    elif epoch == 1950.0:
        theta= 123./180.*sc.pi
        dec_ngp= 27.4/180.*sc.pi
        ra_ngp= 192.25/180.*sc.pi
    else:
        raise IOError("Only epochs 1950 and 2000 are supported")
    return (theta,dec_ngp,ra_ngp)



def radec_to_lb(ra,dec,degree=False,epoch=2000.0):
    """
    NAME:
       radec_to_lb
    PURPOSE:
       transform from equatorial coordinates to Galactic coordinates
    INPUT:
       ra - right ascension
       dec - declination
       degree - (Bool) if True, ra and dec are given in degree and l and b will be as well
       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)
    OUTPUT:
       l,b
       For vector inputs [:,2]
    HISTORY:
       2009-11-12 - Written - Bovy (NYU)
       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)
    """
    import math as m
    import numpy as nu
    import scipy as sc
    #First calculate the transformation matrix T
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    T= sc.dot(sc.array([[sc.cos(theta),sc.sin(theta),0.],[sc.sin(theta),-sc.cos(theta),0.],[0.,0.,1.]]),sc.dot(sc.array([[-sc.sin(dec_ngp),0.,sc.cos(dec_ngp)],[0.,1.,0.],[sc.cos(dec_ngp),0.,sc.sin(dec_ngp)]]),sc.array([[sc.cos(ra_ngp),sc.sin(ra_ngp),0.],[-sc.sin(ra_ngp),sc.cos(ra_ngp),0.],[0.,0.,1.]])))
    #Whether to use degrees and scalar input is handled by decorators
    XYZ= nu.array([nu.cos(dec)*nu.cos(ra),
                   nu.cos(dec)*nu.sin(ra),
                   nu.sin(dec)])
    galXYZ= nu.dot(T,XYZ)
    b= nu.arcsin(galXYZ[2])
    l= nu.arctan2(galXYZ[1]/sc.cos(b),galXYZ[0]/sc.cos(b))
    l[l<0.]+= 2.*nu.pi
    out= nu.array([l,b])
    return out.T



def lb_to_phi12(l,b,degree=False):
    """
    NAME:
       lb_to_phi12
    PURPOSE:
       Transform Galactic coordinates (l,b) to (phi1,phi2)
    INPUT:
       l - Galactic longitude (rad or degree)
       b - Galactic latitude (rad or degree)
       degree= (False) if True, input and output are in degrees
    OUTPUT:
       (phi1,phi2) for scalar input
       [:,2] array for vector input
    HISTORY:
        2014-11-04 - Written - Bovy (IAS)
    """
    import numpy
    from   galpy.util import bovy_coords
    
    _TKOP= numpy.zeros((3,3))
    _TKOP[0,:]= [-0.4776303088,-0.1738432154,0.8611897727]
    _TKOP[1,:]= [0.510844589,-0.8524449229,0.111245042]
    _TKOP[2,:]= [0.7147776536,0.4930681392,0.4959603976]
    #First convert to ra and dec
    radec= bovy_coords.lb_to_radec(l,b)
    ra= radec[:,0]
    dec= radec[:,1]
    XYZ= numpy.array([numpy.cos(dec)*numpy.cos(ra),
                      numpy.cos(dec)*numpy.sin(ra),
                      numpy.sin(dec)])
    phiXYZ= numpy.dot(_TKOP,XYZ)
    phi2= numpy.arcsin(phiXYZ[2])
    phi1= numpy.arctan2(phiXYZ[1],phiXYZ[0])
    phi1[phi1<0.]+= 2.*numpy.pi
    return numpy.array([phi1,phi2]).T


