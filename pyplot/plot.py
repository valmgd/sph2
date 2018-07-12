import numpy as np
import matplotlib.pyplot as plt

def plot_particles_2D(x, P, title='', xlab='$x$ [m]', ylab='$y$ [m]') :
    fig, ax = plt.subplots()

    im = ax.scatter(x[:, 0], x[:, 1], s=20, c=P, cmap='jet', label='Particules')
    cax = fig.colorbar(im, ax=ax)
    cax.set_label('Pressure [Pa]')

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_aspect('equal', 'datalim')
    ax.set_facecolor('lavender')
    ax.set_title(title)

    ax.grid(color='gray', linestyle='--')

    return(fig, ax)
#}
