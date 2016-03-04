# Test consistency between galpy and NEMO
from __future__    import print_function, division
import os
import numpy
import subprocess
from   galpy.orbit import Orbit
from   galpy       import potential
from   galpy.util  import bovy_conversion, bovy_coords

def test_nemo_MN3ExponentialDiskPotential():
    mn= potential.MN3ExponentialDiskPotential(normalize=1.,hr=0.5,hz=0.1)
    tmax= 3.
    vo,ro= 215., 8.75
    o= Orbit([1.,0.1,1.1,0.3,0.1,0.4],ro=ro,vo=vo)
    run_orbitIntegration_comparison(o,mn,tmax,vo,ro)
    return None

def test_nemo_MiyamotoNagaiPotential():
    mp= potential.MiyamotoNagaiPotential(normalize=1.,a=0.5,b=0.1)
    tmax= 4.
    vo,ro= 220., 8.
    o= Orbit([1.,0.1,1.1,0.3,0.1,0.4],ro=ro,vo=vo)
    run_orbitIntegration_comparison(o,mp,tmax,vo,ro)
    return None

def test_nemo_NFWPotential():
    np= potential.NFWPotential(normalize=1.,a=3.)
    tmax= 3.
    vo,ro= 200., 7.
    o= Orbit([1.,0.5,1.3,0.3,0.1,0.4],ro=ro,vo=vo)
    run_orbitIntegration_comparison(o,np,tmax,vo,ro)
    return None

def test_nemo_PowerSphericalPotentialwCutoffPotential():
    pp= potential.PowerSphericalPotentialwCutoff(normalize=1.,alpha=1.,rc=0.4)
    tmax= 2.
    vo,ro= 180., 9.
    o= Orbit([1.,0.03,1.03,0.2,0.1,0.4],ro=ro,vo=vo)
    run_orbitIntegration_comparison(o,pp,tmax,vo,ro)
    return None

def test_nemo_LogarithmicHaloPotential():
    lp= potential.LogarithmicHaloPotential(normalize=1.)
    tmax= 2.
    vo,ro= 210., 8.5
    o= Orbit([1.,0.1,1.1,0.3,0.1,0.4],ro=ro,vo=vo)
    run_orbitIntegration_comparison(o,lp,tmax,vo,ro,tol=0.03)
    return None

def test_nemo_PlummerPotential():
    pp= potential.PlummerPotential(normalize=1.,b=2.)
    tmax= 3.
    vo,ro= 213., 8.23
    xnow, ynow, znow    = numpy.array([12.4, 1.5, 7.1])
    vxnow, vynow, vznow = numpy.array([107., -243., -105.])
    Ri, phii, zcyl  = bovy_coords.rect_to_cyl(xnow, ynow, znow)
    vri, vti, vzcyl = bovy_coords.rect_to_cyl_vec(-vxnow,-vynow,-vznow, xnow, ynow, znow, cyl = False)
    o = Orbit([Ri/ro, vri/vo, vti/vo, zcyl/ro, vzcyl/vo, phii], ro = ro, vo = vo)
    run_orbitIntegration_comparison(o,pp,tmax,vo,ro,tol=0.03)
    return None

def test_nemo_MWPotential2014(filename,ntot):
    mp    = potential.MWPotential2014
    tmax  = 4.75   # time when its at apocenter
    vo,ro = 220., 8.
    #fact = bovy_conversion.velocity_in_kpcGyr(vo, ro)/vo

    # current GD1 position and velocity in cartesian coordinates
    xnow, ynow, znow    = numpy.array([12.4, 1.5, 7.1])
    vxnow, vynow, vznow = numpy.array([107., -243., -105.])

    # current position and velocity of GD1 in cylindrical coordinates
    Ri, phii, zcyl  = bovy_coords.rect_to_cyl(xnow, ynow, znow)
    vri, vti, vzcyl = bovy_coords.rect_to_cyl_vec(-vxnow,-vynow,-vznow, xnow, ynow, znow, cyl = False)

    # Initializing the orbit
    o = Orbit([Ri/ro, vri/vo, vti/vo, zcyl/ro, vzcyl/vo, phii], ro = ro, vo = vo)
    run_orbitIntegration_comparison(o,mp,tmax,vo,ro,filename,ntot,isList=True)
    return None

def run_orbitIntegration_comparison(orb,pot,tmax,vo,ro,filename,ntot,isList=False,
                                    tol=0.01):
    # time in Gyr
    ts = numpy.linspace(0.,tmax/bovy_conversion.time_in_Gyr(vo,ro),1001)
    orb.integrate(ts,pot)
    
    # initial position and velocity of GD1 based on the given time passed in units
    # of [x] = kpc and [v] = kpc/Gyr
    xval  = orb.x(ts[-1])
    yval  = orb.y(ts[-1])
    zval  = orb.z(ts[-1])
    vxval = -orb.vx(ts[-1], use_physical=False)*bovy_conversion.velocity_in_kpcGyr(vo,ro) 
    vyval = -orb.vy(ts[-1], use_physical=False)*bovy_conversion.velocity_in_kpcGyr(vo,ro) 
    vzval = -orb.vz(ts[-1], use_physical=False)*bovy_conversion.velocity_in_kpcGyr(vo,ro) 
    
    # saving the initial position and velocity
    numpy.savetxt('orb.dat',
                  numpy.array([[10.**-6.,orb.x(ts[-1]),orb.y(ts[-1]),orb.z(ts[-1]),
                                -orb.vx(ts[-1], use_physical=False)\
                                    *bovy_conversion.velocity_in_kpcGyr(vo,ro),
                                -orb.vy(ts[-1], use_physical=False)\
                                    *bovy_conversion.velocity_in_kpcGyr(vo,ro),
                                -orb.vz(ts[-1], use_physical=False)\
                                    *bovy_conversion.velocity_in_kpcGyr(vo,ro)]]))

    # initial position and velocities from mkking
    convert_from_nemo('/home/bahmanyar/' + filename + '.nemo','/home/bahmanyar/' + filename + '.dat')
    data_init = numpy.loadtxt('/home/bahmanyar/' + filename + '.dat')
    M_init    = data_init[:,0]
    x_init    = data_init[:,1]
    y_init    = data_init[:,2]
    z_init    = data_init[:,3]
    vx_init   = data_init[:,4]
    vy_init   = data_init[:,5]
    vz_init   = data_init[:,6]

    # shifted position and velocity from mkking file
    # using the initial values obtained
    xshift  = xval  + x_init
    yshift  = yval  + y_init
    zshift  = zval  + z_init
    vxshift = vxval + vx_init
    vyshift = vyval + vy_init
    vzshift = vzval + vz_init

    # stacking the shifted values 
    data_shift = numpy.column_stack((M_init, xshift, yshift, zshift, vxshift, vyshift, vzshift))

    # formatting the shifted file to be the same as the format for gyrfalcON reading
    file = open("orb_shifted.dat",'w')
    file.write("#\n# a2s orb_shifted.dat orb_shifted.nemo\n#\n# run Wed Feb 17 10:43:02 2016\n#  by user  bahmanyar\n#  on host  yngve\n#  with pid 128259\n#\n# time: 0.000000\n# Ntot: %d, Nsink:0, Ngas: 0, Nstd: %d\n#\n# 'mass ' 'pos  ' 'vel  '\n#\n" %(ntot,ntot))

    for i in range(ntot):
    	s1 = str('%.6E' %(data_shift[i][0]))
    	s2 = str('%.6E' %(data_shift[i][1]))
    	s3 = str('%.6E' %(data_shift[i][2]))
    	s4 = str('%.6E' %(data_shift[i][3]))
    	s5 = str('%.6E' %(data_shift[i][4]))
    	s6 = str('%.6E' %(data_shift[i][5]))
    	s7 = str('%.6E' %(data_shift[i][6]))

    	file.write('%15s%15s%15s%15s%15s%15s%15s\n' %(s1,s2,s3,s4,s5,s6,s7))

    file.close()

    # converting the shifted output to nemo format
    os.system('a2s in=orb_shifted.dat out=orb_shifted.nemo N=%d read=mxv' %(ntot))

    # running gyrfalcON
    os.system('gyrfalcON in=orb_shifted.nemo out=gd1_evol_python_test.nemo tstop=' + str(tmax) + ' eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=PowSphwCut+MiyamotoNagai+NFW accpars=0,1001.79126907,1.8,1.9#0,306770.418682,3.0,0.28#0,16.0,162.958241887')
	
    # converting nemo output to dat file
    os.system('s2a gd1_evol_python_test.nemo gd1_evol_python_test.dat')

    

def convert_to_nemo(infile,outfile):
    subprocess.check_call(['a2s','in=%s'% infile,'out=%s' % outfile,'N=10000',
                           'read=mxv'])
    
def convert_from_nemo(infile,outfile):
    subprocess.check_call(['s2a','in=%s' % infile,'out=%s' % outfile])
    
def integrate_gyrfalcon(infile,outfile,tmax,nemo_accname,nemo_accpars):
    """Integrate a snapshot in infile until tmax in Gyr, save to outfile"""
    '''
    with open('gyrfalcON.log','w') as f:
        subprocess.check_call(['gyrfalcON',
                               'in=%s' % infile,
                               'out=%s' % outfile,
                               'tstop=%g' % tmax,
                               'eps=0.0015',
                               'step=0.01',
                               'kmax=10',
                               'Nlev=8',
                               'fac=0.01',
                               'accname=%s' % nemo_accname,
                               'accpars=%s' % nemo_accpars],
                              stdout=f)
    return None
    '''
    os.system('gyrfalcON' + ' ' + 'in=orb.nemo' + ' ' + 'out=orb_evol.nemo'  + ' ' + 'tstop=' + str(tmax) + ' ' + 'eps=0.0015' + ' ' + 'step=0.01' + ' ' +'kmax=10' + ' ' + 'Nlev=8' + ' ' + 'fac=0.01' + ' ' + 'accname=' + str(nemo_accname) + ' ' + 'accpars=' + str(nemo_accpars))

    os.system('s2a' + ' ' + 'orb_evol.nemo' + ' ' + 'orb_evol.dat')



