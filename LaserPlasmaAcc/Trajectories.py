import math as m
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy import constants as cs
from scipy.integrate import odeint
import pristinifier as ps

[tableau20, tableau20Edge] = ps.rgb_array()
font = {'family' : 'Tahoma',
        'weight' : 'bold',
        'size'   : 11}

matplotlib.rc('font', **font)
##############################################
#        INPUT                               #
##############################################

n_e_cm = 1e19 # in cm-3
a0 = 3.
zeta_s = 4.
zeta_0 = 0.0
wavelength = 800*1e-9

##############################################
#    CONVERTED AND PHYS. PARAM               #
##############################################

n_e = n_e_cm*1e6 # converted to m-3
omega_p = m.sqrt(n_e*cs.e*cs.e/(cs.m_e*cs.epsilon_0))
omega_0 = 2*3.14*cs.c/wavelength
beta_p = cs.c # underdense plasma
k_p = omega_p/cs.c
v_p = cs.c*(1+0.5*m.pow(omega_p/omega_0, 2))

print 1/k_p
##############################################
#    FUNCTIONS                               #
##############################################

# Hamiltonian
def H(_gamma, _beta_p, _phi):
    return _gamma-_beta_p*m.sqrt(_gamma*_gamma-1)-_phi

# Electrostatic potential
def phi_def(_zeta):
    return m.cos(_zeta)

# Field
def a_zeta(_zeta):
    return a0*m.exp(-((_zeta-zeta_0)/zeta_s)*((_zeta-zeta_0)/zeta_s))
    #return a0


##############################################
#    SOLVE for electrostatic potential       #
##############################################

# Make time array for solution
tStop = 20.
tInc = tStop/1000
t = np.arange(-10.0, tStop, tInc)


def f(y, t):
    phi, x = y      # unpack current values of y
    derivs = [x,      # list of dy/dt=f functions
              0.5*((1. + (a_zeta(t)*a_zeta(t))/2.)/((1. + phi)*(1. + phi)) - 1.)]
    return derivs


# Initial values
phi0 = 0.0     # initial field
x0 = 0.0       # initial field derivative wrt zeta


# Bundle initial conditions for ODE solver
y0 = [phi0, x0]

# Call the ODE solver
psoln = odeint(f, y0, t, hmax=2*tInc)


# Plot results
fig = plt.figure(1, figsize=(8,6), dpi=80, facecolor='w', edgecolor='k', linewidth=5)


# Plot theta as a function of time
ax1 = fig.add_subplot(211)
ax1.plot(t, psoln[:, 0], linewidth=2, color=tableau20[0])
ax1.plot(t, [a_zeta(zeta) for zeta in t], linewidth=2, color=tableau20[4])
ax1.set_xlabel('$k_p\zeta$', fontsize=15)
ax1.set_ylabel('$\phi(k_p\zeta)$', fontsize=15)

# Plot derivative as a function of zetakp
ax2 = fig.add_subplot(212)
ax2.plot(t, -np.gradient(psoln[:, 0]), linewidth=2, color=tableau20[2])
ax2.set_xlabel('$k_p\zeta$', fontsize=15)
ax2.set_ylabel('$E(k_p\zeta )$', fontsize=15)

plt.tight_layout()
plt.show()

