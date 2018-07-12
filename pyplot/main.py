#!/home/vmagda/.anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from plot import *

# -----------------------------------------------------------------------------------------------------------
# load data
# -----------------------------------------------------------------------------------------------------------
x       = np.loadtxt('../sorties/x.dat')
w       = np.loadtxt('../sorties/w.dat')

rho     = np.loadtxt('../sorties/rho.dat')
u       = np.loadtxt('../sorties/u.dat')
P       = np.loadtxt('../sorties/P.dat')

nor     = np.loadtxt('../sorties/nor.dat')
div_nor = np.loadtxt('../sorties/div_nor.dat')

grad_P  = np.loadtxt('../sorties/grad_P.dat')
fts     = np.loadtxt('../sorties/fts.dat')
dmu     = np.loadtxt('../sorties/dmu.dat')



# -----------------------------------------------------------------------------------------------------------
# tracés
# -----------------------------------------------------------------------------------------------------------
fig, ax = plot_particles_2D(x, P)
ax.quiver(x[:, 0], x[:, 1], dmu[:, 0], dmu[:, 1],
    angles='xy', scale=0.05, label='Mvt quantity', color='black')
ax.legend(loc='best', fontsize='x-large', fancybox=True, framealpha=0.5)
fig.tight_layout()



fig, ax = plot_particles_2D(x, P)
ax.quiver(x[:, 0], x[:, 1], nor[:, 0], nor[:, 1],
    angles='xy', scale=50, label='Mvt quantity', color='black')
ax.legend(loc='best', fontsize='x-large', fancybox=True, framealpha=0.5)
fig.tight_layout()



# -----------------------------------------------------------------------------------------------------------
# show
# -----------------------------------------------------------------------------------------------------------
plt.show()
