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
#       OBJECT AND IMAGE SPACE               #
##############################################

fig = plt.figure(1, figsize=(10,10), dpi=80, facecolor='w', edgecolor='k', frameon=False)
plt.xticks([])
plt.yticks([])
# Plot theta as a function of time
ax1 = fig.add_subplot(111)

ax1.fill_between(x, 0, 1,  facecolor='green', alpha=0.5)

ax1.set_xlabel('$g1$', fontsize=30)
ax1.set_ylabel('$g2$', fontsize=30)
ax1.set_ylim(-1., 1.)
ax1.set_xlim(-1., 1.)
ax1.axhline(0, color='black')
ax1.axvline(0, color='black')
plt.tight_layout()
#plt.show()
#plt.savefig(my_path + 'espaces.png')


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
#       SPECTRUM OF LONGITUDINAL MODES       #
##############################################



