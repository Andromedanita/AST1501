import numpy            as    np
import matplotlib.pylab as    plt
import galpy
from   galpy.util      import bovy_conversion
from   nemo_funcs      import tail_cut

fact  = bovy_conversion.freq_in_Gyr(220., 8.)
fact2 = bovy_conversion.freq_in_kmskpc(220., 8.)

prog = np.loadtxt("/Users/anita/dyn-modeling-streams-2014/sim/test.dat", delimiter=',')
tail = np.loadtxt("/Users/anita/Desktop/Meeting_2/tail_concatenated_final.txt")

JR_p, Lz_p, Jz_p, omegaR_p, omega_phi_p, omegaz_p, thetaR_p, theta_phi_p, thetaz_p = prog[:,0], prog[:,1], prog[:,2], prog[:,3], prog[:,4], prog[:,5], prog[:,6], prog[:,7], prog[:,8]
JR_t, Lz_t, Jz_t, omegaR_t, omega_phi_t, omegaz_t, thetaR_t, theta_phi_t, thetaz_t = tail[:,0], tail[:,1], tail[:,2], tail[:,3], tail[:,4], tail[:,5], tail[:,6], tail[:,7], tail[:,8]

omegaR_t, omegaz_t = tail_cut(tail)

plt.ion()
plt.plot(omegaR_t*fact, omegaz_t * fact, '.', ms=2)

''''
thetar = thetaR_t

thetar = (np.pi+(thetar-np.median(thetar))) % (2.*np.pi)
indx   = np.fabs(thetar-np.pi) > (5.*np.median(np.fabs(thetar-np.median(thetar))))
plotx  = tail[indx,3]*bovy_conversion.freq_in_Gyr(220.,8.)
ploty  = tail[indx,5]*bovy_conversion.freq_in_Gyr(220.,8.)
'''




plt.ion()
plt.figure(1)

plt.subplot(231)
plt.plot(tail[:,0], tail[:,2], 'k.', ms=2)
plt.plot(prog[:,0], prog[:,2], 'ro')
plt.xlabel("$J_R$")
plt.ylabel("$J_z$")

plt.subplot(232)
plt.plot(omegaR_t*fact, omegaz_t * fact, 'k.', ms=2)
plt.plot(prog[:,3]*fact, prog[:,5]*fact, 'ro')
plt.xlabel("$\Omega_R \mathrm{(Gyr^{-1})}$")
plt.ylabel("$\Omega_z \mathrm{(Gyr^{-1})}$")
plt.xlim(15.45,16.)
plt.ylim(11.70,12.05)

plt.subplot(233)
plt.plot(tail[:,6], tail[:,8], 'k.', ms=2)
plt.plot(prog[:,6], prog[:,8], 'ro')
plt.xlabel(r"$\theta_R$")
plt.ylabel(r"$\theta_z$")

plt.subplot(234)
plt.plot(tail[:,0], tail[:,1], 'k.', ms=2)
plt.plot(prog[:,0], prog[:,1], 'ro')
plt.xlabel("$J_R$")
plt.ylabel("$L_z$")

plt.subplot(235)
plt.plot(tail[:,3]*fact, tail[:,4]*fact, 'k.', ms=2)
plt.plot(prog[:,3]*fact, prog[:,4]*fact, 'ro')
plt.xlabel("$\Omega_R \mathrm{(Gyr^{-1})}$")
plt.ylabel("$\Omega_{\phi} \mathrm{(Gyr^{-1})}$")
plt.xlim(15.45,16.)
plt.ylim(-11.,-10.65)

plt.subplot(236)
plt.plot(tail[:,6], tail[:,7], 'k.', ms=2)
plt.plot(prog[:,6], prog[:,7], 'ro')
plt.xlabel(r"$\theta_R$")
plt.ylabel(r"$\theta_{\phi}$")


plt.suptitle("$W = 2.0$")


