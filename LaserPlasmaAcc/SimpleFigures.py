import math as m
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy import constants as cs
from scipy.integrate import odeint
import pristinifier as ps
import cmath as cm
from matplotlib.widgets import Cursor


##############################################
#             PATHS AND SETTINGS             #
##############################################

[tableau20, tableau20Edge] = ps.rgb_array()
font = {'family' : 'Tahoma',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

my_path = "C:\Users\\vernier\PycharmProjects\octo-elective\\"




##############################################
#       STABILITY DIAGRAM LINEAR CAVITY      #
##############################################

# def g2_f(_g1):
#     return 1/_g1
#
# g1_t = np.append(np.arange(-10, -0.001, 0.001), np.arange(0.001, 10, 0.001))
# g2_t = [g2_f(g1) for g1 in g1_t]
#
# fig = plt.figure(1, figsize=(10,10), dpi=80, facecolor='w', edgecolor='k', frameon=False)
# plt.xticks([])
# plt.yticks([])
# # Plot theta as a function of time
# ax1 = fig.add_subplot(111)
# ax1.plot(g1_t, g2_t, linewidth=2, color=tableau20[2])
# ax1.fill_between(g1_t, 0, g2_t, linewidth=2, color=tableau20[3])
#
# ax1.set_xlabel('$g1$', fontsize=30)
# ax1.set_ylabel('$g2$', fontsize=30)
# ax1.set_ylim(-1., 1.)
# ax1.axhline(0, color='black')
# ax1.axvline(0, color='black')
# plt.tight_layout()
# #plt.show()
# plt.savefig(my_path + 'stability.png')
#


##############################################
#       OSCILLATOR IN PHASE SPACE            #
##############################################

# tmin = 0
# tmax = 10.
# tnum = 1000
# tstep = (tmax - tmin)/tnum
#
# t_array = np.arange(tmin, tmax, tstep)
#
#
# fig = plt.figure(1, figsize=(10,10), dpi=80, facecolor='w', edgecolor='k')
#
# # Plot theta as a function of time
# ax1 = fig.add_subplot(111)
# for a0 in np.arange(1, 10, 1):
#     x = [a0*m.cos(t) for t in t_array]
#     y = [a0*m.sin(t) for t in t_array]
#     ax1.plot(x, y, linewidth=2, color=tableau20[a0])
#
# ax1.set_xlabel('$x$', fontsize=30)
# ax1.set_ylabel('$\\frac{p^2(t)}{m\omega_0}$', fontsize=30)
# plt.tight_layout()
# plt.savefig(my_path + 'simple_phase_space.png')

##############################################
#    DAMPED OSCILLATOR IN PHASE SPACE, BLANK #
# ##############################################
#
#
# tmin = 0
# tmax = 10.
# tnum = 1000
# tstep = (tmax - tmin)/tnum
# #
# t_array = np.arange(tmin, tmax, tstep)
#
#
# fig = plt.figure(1, figsize=(11,10), dpi=80, facecolor='w', edgecolor='k')
#
# # Plot theta as a function of time
# ax1 = fig.add_subplot(111)
#
# # for a0 in np.arange(1, 10, 1):
# #     x = [a0*m.cos(t) for t in t_array]
# #     y = [a0*m.sin(t) for t in t_array]
# #     ax1.plot(x, y, linewidth=2, color=tableau20[a0])
# #
# ax1.set_xlim(-1, 1)
# ax1.set_ylim(-1, 1)
# ax1.set_xlabel('$x$', fontsize=30)
# ax1.set_ylabel('$\\frac{p^2(t)}{m\omega_0}$', fontsize=30)
# plt.tight_layout()
# plt.savefig(my_path + 'empty_phase_space.png')

##############################################
#       SPECTRUM OF LONGITUDINAL MODES       #
##############################################
# c = 3e8
#
#
# wlambda = 632.8e-9
# f_mid = c/wlambda
# fmin = f_mid*(1.+6.e-6)
# fmax = f_mid*(1.-6.e-6)
#
# lw = f_mid*1e-8
#
# fnum = 10000
# fstep = (fmax - fmin)/fnum
#
# f_array = np.arange(fmin, fmax, fstep)
#
#
# R1 = 10.
# R2 = 20.
# L = 0.1
#
# p_mid = int(2*L/wlambda)
#
# A = 1-2L/R1
# B = 2L*(1-L/R2)
# z = L*(L-R2)/(R1+R2-2*L)
# zR = m.sqrt(L*(R1+R2-L)*(R1-L)*(R2-L)/(R1+R2-2*L)/(R1+R2-2*L))
#
# def w_sq(_z):
#         return zR*wlambda/3.14*(1+m.pow(_z/zR, 2))
#
#
# def dw(_p, _m, _n):
#         q_inv = 1 / R1 - 1j * wlambda / (3.14 * w_sq(z))
#         return c/(2*L)*(_p-(_n+_m+1)/(2*3.14)*cm.phase(A+B/q_inv))
#
#
# def lor(_f, _f0):
#         return 1/(1+pow((_f-_f0)/lw, 2))
#
# spectrum = np.zeros(fnum)
# p_lim = 4
# n_m_lim = 4
# for p_i in range(p_mid-p_lim, p_mid+p_lim):
#         for m_i in range(0, n_m_lim):
#                 for n_i in range(0,1):
#                         for ii in range(0, fnum):
#                                 spectrum[ii] += lor(f_array[ii], dw(p_i, m_i, n_i))\
#                                                 *m.pow(1/(3.14*m.factorial(m_i)*m.factorial(n_i)), 0.5)\
#                                                 *m.pow(2, -(n_i+m_i)/2)
#
#
# fig = plt.figure(1, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
# ax1 = fig.add_subplot(111)
# ax1.plot(f_array*1e-12, spectrum/max(spectrum), linewidth=2, color=tableau20[0])
# #cursor = Cursor(ax1, useblit=True, color='red', linewidth=1)
#
# ax1.set_xlabel('Frequence (THz)', fontsize=20)
# ax1.set_ylabel('Amplitude du mode (a.u.)', fontsize=20)
# ax1.set_ylim(0, 1.05)
# ax1.ticklabel_format(useOffset=False)
# ax1.set_xlim(fmin*1e-12, fmax*1e-12)
# #plt.show()
# plt.savefig(my_path + 'ModeAmplitude.png')
#
##############################################
#    COULOMB POTENTIAL                       #
# ############################################

rmax = 10.
rmin = 0.6
rnum = 1000
rstep = rmax/rnum

alpha_1 = 10.
r_array = np.append(np.arange(-rmax, -rmin, rstep), np.arange(rmin, rmax, rstep))


fig = plt.figure(1, figsize=(17, 7), dpi=80, facecolor='w', edgecolor='k')

# Plot theta as a function of time
ax1 = fig.add_subplot(111)

def coulomb(r, alpha_2):
    return -alpha_1/m.sqrt(m.pow(r, 2)) - alpha_2*r

a2_arr = [0., 0.5, 1., 1.5, 2., 2.5]
for ii in range(0, 4):
    a2 = a2_arr[ii]
    c_p = [coulomb(r, a2) for r in r_array]
    ax1.plot(r_array, c_p, linewidth=4, color=tableau20[ii], label='eE = '+str(a2)+'$\,$(u.a.)')

ax1.set_xlim(-rmax, rmax)
ax1.set_ylim(-13, 0)
ax1.set_xlabel('$r\,\mathrm{(u.a.)}$', fontsize=25)
ax1.set_ylabel('$V(r)\,\mathrm{(u.a.)}$', fontsize=25)
#ax1.legend(loc='lower left')
#plt.show()
#plt.tight_layout()
plt.savefig(my_path + 'suppression_barriere.png')

