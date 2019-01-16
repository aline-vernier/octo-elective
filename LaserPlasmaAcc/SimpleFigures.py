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
#       OSCILLATOR IN PHASE SPACE            #
##############################################

tmin = 0
tmax = 10.
tnum = 1000
tstep = (tmax - tmin)/tnum

t_array = np.arange(tmin, tmax, tstep)
x = [m.cos(t) for t in t_array]
y = [m.sin(t) for t in t_array]

fig = plt.figure(1, figsize=(10,10), dpi=80, facecolor='w', edgecolor='k')

# Plot theta as a function of time
ax1 = fig.add_subplot(111)
ax1.plot(x, y, linewidth=2, color=tableau20[0])

ax1.set_xlabel('$x$', fontsize=30)
ax1.set_ylabel('$\\frac{p^2(t)}{m\omega_0}$', fontsize=30)
plt.tight_layout()
plt.savefig(my_path + 'simple_phase_space.png')