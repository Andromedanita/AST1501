import numpy            as    np
import matplotlib.pylab as    plt
import galpy
from   galpy.util      import bovy_conversion
from   nemo_funcs      import tail_cut

R0 = 8.
V0 = 220.

fact  = bovy_conversion.freq_in_Gyr(220., 8.)
fact2 = bovy_conversion.freq_in_kmskpc(220., 8.)

prog = np.loadtxt("/Users/anita/dyn-modeling-streams-2014/sim/test.dat", delimiter=',')
tail = np.loadtxt("/Users/anita/Desktop/Meeting_2/tail_concatenated_final.txt")

JR_p, Lz_p, Jz_p, omegaR_p, omega_phi_p, omegaz_p, thetaR_p, theta_phi_p, thetaz_p = prog[:,0], prog[:,1], prog[:,2], prog[:,3], prog[:,4], prog[:,5], prog[:,6], prog[:,7], prog[:,8]
JR_t, Lz_t, Jz_t, omegaR_t, omega_phi_t, omegaz_t, thetaR_t, theta_phi_t, thetaz_t = tail[:,0], tail[:,1], tail[:,2], tail[:,3], tail[:,4], tail[:,5], tail[:,6], tail[:,7], tail[:,8]

indx_t = tail_cut(tail)
indx_p = tail_cut(prog)

# Tail
omegaR_t, omega_phi_t, omegaz_t,  = tail[indx_t,3], tail[indx_t,4], tail[indx_t,5]
JR_t, Lz_t, Jz_t                  = tail[indx_t,0], tail[indx_t,1], tail[indx_t,2]
thetaR_t, theta_phi_t, thetaz_t   = tail[indx_t,6], tail[indx_t,7], tail[indx_t,8]

thetaR_t    = (np.pi+(thetaR_t-np.median(tail[:,6]))) % (2.*np.pi)
theta_phi_t = (np.pi+(theta_phi_t-np.median(tail[:,7]))) % (2.*np.pi)
thetaz_t    = (np.pi+(thetaz_t-np.median(tail[:,8]))) % (2.*np.pi)

'''
# Progenitor
omegaR_p, omega_phi_p, omegaz_p,  = prog[indx_p,3], prog[indx_p,4], prog[indx_p,5]
JR_p, Lz_p, Jz_p                  = prog[indx_p,0], prog[indx_p,1], prog[indx_p,2]
thetaR_p, theta_phi_p, thetaz_p   = prog[indx_p,6], prog[indx_p,7], prog[indx_p,8]

thetaR_p    = (np.pi+(thetaR_p-np.median(prog[:,6]))) % (2.*np.pi)
theta_phi_p = (np.pi+(theta_phi_p-np.median(prog[:,7]))) % (2.*np.pi)
thetaz_p    = (np.pi+(thetaz_p-np.median(prog[:,8]))) % (2.*np.pi)
'''
JR_p, Lz_p, Jz_p                = prog[:,0][len(prog)-1], prog[:,1][len(prog)-1], prog[:,2][len(prog)-1]
omegaR_p, omega_phi_p, omegaz_p = prog[:,3][len(prog)-1], prog[:,4][len(prog)-1], prog[:,5][len(prog)-1]
thetaR_p, theta_phi_p, thetaz_p = prog[:,6][len(prog)-1]+np.pi, prog[:,7][len(prog)-1]+np.pi, prog[:,8][len(prog)-1]+np.pi


plt.ion()
plt.plot(omegaR_t*fact, omegaz_t * fact, '.', ms=2)


plt.ion()
plt.figure(1)

plt.subplot(231)
plt.plot(JR_t * R0, Jz_t * R0, 'k.', ms=2)
plt.plot(JR_p * R0, Jz_p * R0, 'ro')
plt.xlabel("$J_R$")
plt.ylabel("$J_z$")

plt.subplot(232)
plt.plot(omegaR_t * fact, omegaz_t * fact, 'k.', ms=2)
plt.plot(omegaR_p * fact, omegaz_p *fact, 'ro')
plt.xlabel("$\Omega_R \mathrm{(Gyr^{-1})}$")
plt.ylabel("$\Omega_z \mathrm{(Gyr^{-1})}$")
plt.xlim(15.45,16.)
plt.ylim(11.70,12.05)

plt.subplot(233)
plt.plot(thetaR_t, thetaz_t, 'k.', ms=2)
#plt.plot(thetaR_p, thetaz_p, 'ro')
plt.xlabel(r"$\theta_R$")
plt.ylabel(r"$\theta_z$")

plt.subplot(234)
plt.plot(JR_t * R0, Lz_t * R0, 'k.', ms=2)
plt.plot(JR_p * R0, Lz_p * R0, 'ro')
plt.xlabel("$J_R$")
plt.ylabel("$L_z$")

plt.subplot(235)
plt.plot(omegaR_t * fact, omega_phi_t * fact, 'k.', ms=2)
plt.plot(omegaR_p * fact, omega_phi_p * fact, 'ro')
plt.xlabel("$\Omega_R \mathrm{(Gyr^{-1})}$")
plt.ylabel("$\Omega_{\phi} \mathrm{(Gyr^{-1})}$")
plt.xlim(15.45,16.)
plt.ylim(-11.,-10.65)

plt.subplot(236)
plt.plot(thetaR_t, theta_phi_t, 'k.', ms=2)
#plt.plot(thetaR_p, theta_phi_p, 'ro')
plt.xlabel(r"$\theta_R$")
plt.ylabel(r"$\theta_{\phi}$")


plt.suptitle("$W = 2.0$")


