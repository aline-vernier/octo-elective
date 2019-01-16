import math as m
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy import constants as cs
from scipy.integrate import odeint
import pristinifier as ps

##############################################
#             PATHS AND SETTINGS             #
##############################################

[tableau20, tableau20Edge] = ps.rgb_array()
font = {'family' : 'Tahoma',
        'weight' : 'bold',
        'size'   : 11}

matplotlib.rc('font', **font)

my_path = "C:\Users\\vernier\PycharmProjects\octo-elective\\"
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

gamma_p = omega_0/omega_p
beta_p = 1 - 0.5*m.pow(omega_p/omega_0, 2)

k_p = omega_p/cs.c
v_p = cs.c*(1+0.5*m.pow(omega_p/omega_0, 2))

##############################################
#    FUNCTIONS                               #
##############################################

# Hamiltonian
def H(_gamma, _beta_p, _phi):
    return _gamma-_beta_p*m.sqrt(_gamma*_gamma-1)-_phi

# Solution for gamma
def gamma (_H0, _phi):
    if m.pow(gamma_p, 2)*m.pow(_H0 + _phi, 2) - 1 < 0:
        print "H0 = " + str(_H0)
        print "phi = " + str(_phi)
        return None, None
    else:
        g0 = gamma_p*gamma_p*(_H0 + _phi)
        g1 = m.sqrt(m.pow(gamma_p, 2)-1) * m.sqrt(m.pow(gamma_p, 2)*m.pow(_H0 + _phi, 2) - 1)
        sol1 = g0 - g1
        sol2 = g0 + g1
        return  sol1, sol2

# Electrostatic potential
def phi_def(_zeta):
    return m.cos(_zeta)

# Field
def a_zeta(_zeta):
    return a0*m.exp(-((_zeta-zeta_0)/zeta_s)*((_zeta-zeta_0)/zeta_s))


##############################################
#    SOLVE for electrostatic potential       #
##############################################

# Make time array for solution
tStart = -10.0
tStop = 20.
tNum = 1000
tInc = tStop/tNum
t = np.arange(tStart, tStop, tInc)


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

##############################################
#    PLOT      electrostatic potential       #
##############################################

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
plt.savefig(my_path + '1D_Field_and_Potential.png')

##############################################
#    CREATE TRAJECTORIES IN STATIC FIELD     #
##############################################


phi_zeta_kp_sol = psoln[:, 0]
zeta_kp = t

i0_min = 1
i0_max = tNum - 1
i0_num = 10
i0_step = int((i0_max - i0_min)/i0_num)

gamma_min = 5.
gamma_max = 15.
gamma_num = 10


gamma_step = (gamma_max - gamma_min)/gamma_num


# A few trajectories
fig2 = plt.figure(2, figsize=(8,6), dpi=80, facecolor='w', edgecolor='k', linewidth=5)
ax1 = fig2.add_subplot(111)

for gamma_0 in np.arange(gamma_min, gamma_max, gamma_step):
    for i0 in np.arange(i0_min, i0_max, i0_step):
        phi_0 = phi_zeta_kp_sol[i0]

        H0 = H(gamma_0, beta_p, phi_0)
        s1 = []
        s2 = []
        for i in np.arange(i0, len(zeta_kp)):
            phi = phi_zeta_kp_sol[i]
            a1, a2 = gamma(H0, phi)
            s1.append(a1)
            s2.append(a2)

        ax1.plot(t[i0:], s1, linewidth=1, color=tableau20[0])
        ax1.plot(t[i0:], s2, linewidth=1, color=tableau20[2])


plt.show()

