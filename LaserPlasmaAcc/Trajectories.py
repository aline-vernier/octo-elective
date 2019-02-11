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
 #       'weight' : 'bold',
        'size'   : 15}

matplotlib.rc('font', **font)

my_path = "C:\Users\\vernier\PycharmProjects\octo-elective\\"
##############################################
#        INPUT                               #
##############################################

n_e_cm = 1e19 # in cm-3
a0 = .1
zeta_s = 2.
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

print "Gamma p = " + str(gamma_p)
print "Beta p = " + str(beta_p)
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
        return None, None
    else:
        g0 = gamma_p*gamma_p*(_H0 + _phi)
        g1 = m.sqrt(m.pow(gamma_p, 2)-1) * m.sqrt(m.pow(gamma_p, 2)*m.pow(_H0 + _phi, 2) - 1)
        sol1 = g0 - g1
        sol2 = g0 + g1
        return  sol1, sol2


def uz (_H0, _phi, _a):
    gamma_t_sq = 1 + _a*_a
    if m.pow(gamma_p, 2)*pow((_H0 + _phi), 2) < gamma_t_sq :
        return None, None
    g0 = gamma_p * gamma_p * beta_p*(_H0 + _phi)
    g1 = gamma_p*m.sqrt(
             m.pow(gamma_p, 2)*pow((_H0 + _phi), 2) - gamma_t_sq
         )
    sol1 = g0 - g1
    sol2 = g0 + g1
    return sol1 +1., sol2 + 1.

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
tStart = -5.0
tStop = 20.
tNum = 500000
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
fig = plt.figure(1, figsize=(10,6), dpi=80, facecolor='w', edgecolor='k', linewidth=5)

# Plot theta as a function of time
ax1 = fig.add_subplot(211)
ax1.plot(-t/(2*3.14), psoln[:, 0], linewidth=2, color=tableau20[0], label = '$\phi(k_p\zeta)$')
ax1.plot(-t/(2*3.14), [a_zeta(zeta)*a_zeta(zeta) for zeta in t], dashes=[6, 2], linewidth=1, color=tableau20[4], label='$a^2(\zeta)$')

#ax1.set_ylabel('$\phi(k_p\zeta)$', fontsize=15, color=tableau20[0])


ax1.set_xlim(-2.5, 0.5)
ax1.set_ylim(-0.01, 0.011)

# Plot derivative as a function of zetakp
ax2 = fig.add_subplot(212)
ax2.plot(-t/(2*3.14), np.gradient(psoln[:, 0])/max(np.gradient(psoln[:, 0])), linewidth=2, color=tableau20[2], label='$E(k_p\zeta )/E_{Max}$')
ax2.set_xlabel('$k_p\zeta/2\pi$', fontsize=15)
#ax2.set_ylabel('$E(k_p\zeta )/E_{Max}$', fontsize=15, color=tableau20[2])
ax2.set_xlim(-1.7, 0.5)
ax2.set_ylim(-1.5, 1.1)

ax1.legend(loc='upper left')
ax2.legend(loc='upper left')
plt.tight_layout()
plt.savefig(my_path + '1D_Field_and_Potential_linear.png')
#plt.show()

##############################################
#    CREATE TRAJECTORIES IN STATIC FIELD     #
##############################################


phi_zeta_kp_sol = psoln[:, 0]
zeta_kp = t

i0_min = 1
i0_max = tNum - 1
i0_num = 1
i0_step = int((i0_max - i0_min)/i0_num)

gamma_min = 1.
gamma_max = 2.8
gamma_num = 10.


gamma_step = (gamma_max)/gamma_num


# # A few trajectories
# fig2 = plt.figure(2, figsize=(8,6), dpi=80, facecolor='w', edgecolor='k', linewidth=5)
# ax1 = fig2.add_subplot(111)
#
# for gamma_0 in np.arange(gamma_min, gamma_max, gamma_step):
#     for i0 in np.arange(i0_min, i0_max, i0_step):
#         phi_0 = phi_zeta_kp_sol[i0]
#
#         H0 = H(gamma_0, beta_p, phi_0)
#
#
#         s1 = []
#         s2 = []
#
#         for i in np.arange(i0, len(zeta_kp)):
#             phi = phi_zeta_kp_sol[i]
#             a = a_zeta(i)
#             a1, a2 = uz(H0, phi, a)
#
#             s1.append(a1)
#             s2.append(a2)
#
#         ax1.semilogy(-t[i0:]/(2*3.14), s1, dashes=[6, 2], linewidth=1, color=tableau20[2])
#         ax1.semilogy(-t[i0:]/(2*3.14), s2, dashes=[6, 2], linewidth=1, color=tableau20[2])
#
#
# Hsep = 1/gamma_p - min(phi_zeta_kp_sol)
# print "Hsep = " + str(Hsep)
# s1 = []
# s2 = []
# for i in np.arange(i0_min, len(zeta_kp)):
#     phi = phi_zeta_kp_sol[i]
#     a = a_zeta(i)
#     a1, a2 = uz(Hsep, phi, a)
#
#     s1.append(a1)
#     s2.append(a2)
#
# ax1.semilogy(-t[i0:]/(2*3.14), s1, linewidth=2, color=tableau20[6])
# ax1.semilogy(-t[i0:]/(2*3.14), s2, linewidth=2, color=tableau20[6])
#
# Hin = [Hsep/i for i in [2, 3, 5, 50]]
#
# for h in Hin:
#     for i0 in np.arange(i0_min, i0_max, i0_step):
#         phi_0 = phi_zeta_kp_sol[i0]
#
#         H0 = h
#         s1 = []
#         s2 = []
#
#         for i in np.arange(i0, len(zeta_kp)):
#             phi = phi_zeta_kp_sol[i]
#             a = a_zeta(i)
#             a1, a2 = uz(H0, phi, a)
#
#             s1.append(a1)
#             s2.append(a2)
#
#         ax1.semilogy(-t[i0:]/(2*3.14), s1, dashes=[6, 2], linewidth=1, color=tableau20[2])
#         ax1.semilogy(-t[i0:]/(2*3.14), s2, dashes=[6, 2], linewidth=1, color=tableau20[2])
#
# ax1.set_xlim(-1.7, 0.5)
# ax1.set_xlabel('$k_p\zeta/2\pi$', fontsize=15)
# ax1.set_ylabel('$u_z+1$', fontsize=15)
#
#
# #plt.savefig(my_path + 'trajectories.png')

##############################################
#    PLOT TRAJECTORIES WITH SEPERATRIX       #
##############################################
#
# fig3 = plt.figure(2, figsize=(8,6), dpi=80, facecolor='w', edgecolor='k', linewidth=5)
# ax1 = fig3.add_subplot(111)
#
#
# Hsep = 1/gamma_p - min(phi_zeta_kp_sol)
# print "Hsep = " + str(Hsep)
#
#
# s1 = []
# s2 = []
# for i in np.arange(i0_min, len(zeta_kp)):
#     phi = phi_zeta_kp_sol[i]
#     a = a_zeta(i)
#     a1, a2 = uz(Hsep, phi, a)
#
#     s1.append(a1)
#     s2.append(a2)
#
# ax1.semilogy(-t[i0_min:]/(2*3.14), s1, linewidth=2, color=tableau20[6])
# ax1.semilogy(-t[i0_min:]/(2*3.14), s2, linewidth=2, color=tableau20[6])
#
#
# Hin = [Hsep/i for i in [1.1, 1.6, 2, 5, 10, 200]]
#
# for h in Hin:
#     for i0 in np.arange(i0_min, i0_max, i0_step):
#
#         H0 = h
#         s1 = []
#         s2 = []
#
#         for i in np.arange(i0, len(zeta_kp)):
#             phi = phi_zeta_kp_sol[i]
#             a = a_zeta(i)
#             a1, a2 = uz(H0, phi, a)
#
#             s1.append(a1)
#             s2.append(a2)
#
#         ax1.semilogy(-t[i0:]/(2*3.14), s1, dashes=[4, 2], linewidth=1.2, color=tableau20[0])
#         ax1.semilogy(-t[i0:]/(2*3.14), s2, dashes=[4, 2], linewidth=1.2, color=tableau20[0])
#
# ax1.set_xlim(-1.7, 0.5)
# ax1.set_xlabel('$k_p\zeta/2\pi$', fontsize=15)
# ax1.set_ylabel('$u_z+1$', fontsize=15)
#
# Hout = [Hsep*i for i in [1.1, 1.5]]
#
# for h in Hout:
#     for i0 in np.arange(i0_min, i0_max, i0_step):
#
#         H0 = h
#         s1 = []
#         s2 = []
#
#         for i in np.arange(i0, len(zeta_kp)):
#             phi = phi_zeta_kp_sol[i]
#             a = a_zeta(i)
#             a1, a2 = uz(H0, phi, a)
#
#             s1.append(a1)
#             s2.append(a2)
#
#         ax1.semilogy(-t[i0:]/(2*3.14), s1, dashes=[4, 2],linewidth=2, color=tableau20[8])
#   #      ax1.semilogy(-t[i0:]/(2*3.14), s2, dashes=[6, 2], linewidth=1, color=tableau20[8])
#
# for i0 in np.arange(i0_min, i0_max, i0_step):
#
#     H0 = H(1, beta_p,0)
#     s1 = []
#     s2 = []
#
#     for i in np.arange(i0, len(zeta_kp)):
#         phi = phi_zeta_kp_sol[i]
#         a = a_zeta(i)
#         a1, a2 = uz(H0, phi, a)
#
#         s1.append(a1)
#         s2.append(a2)
#
#     ax1.semilogy(-t[i0:]/(2*3.14), s1, linewidth=2, color='black')
#
# ax1.set_xlim(-2.5, 0.5)
# ax1.set_xlabel('$k_p\zeta/2\pi$', fontsize=15)
# ax1.set_ylabel('$u_z+1$', fontsize=15)
# plt.savefig(my_path + 'trajectories.png')