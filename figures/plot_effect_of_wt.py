# -*- coding: utf-8 -*-
"""
plot_effect_of_wt.py -- paper "effect of w_t" panel (file: figure_3.eps).

Three steady-state occupation maps on a 20x20 square lattice at
w_a = 60 w_r, phi = 0.75 and increasing translational diffusion:
    (a) F_vel.dat  w_t = 0
    (b) I_vel.dat  w_t = 5  w_r
    (c) G_vel.dat  w_t = 10 w_r
Requires a LaTeX installation with the kpfonts package.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
plt.rcParams.update({
    # "axes.formatter.use_locale": True,
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["kpfonts"],
    "text.latex.preamble": r'\usepackage{kpfonts} \usepackage[T1]{fontenc}',
    'font.size': 20,
    'axes.labelsize': 18,
    'axes.titlesize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 11,
    'axes.labelpad' : 1
}) 
files = ["F_vel.dat", "I_vel.dat", "G_vel.dat"]
labels = ['(a)', '(b)', '(c)']

fig, axes = plt.subplots(1, 3, figsize=(10, 4))  # narrower figure

images = []

for ax, fname, lab in zip(axes, files, labels):
    data = np.loadtxt(os.path.join(DATA_DIR, fname))

    i = data[:, 0].astype(int) - 1
    j = data[:, 1].astype(int) - 1
    rho = data[:, 2]

    grid = np.zeros((20, 20))
    grid[j, i] = rho

    im = ax.imshow(grid, origin='lower', cmap='Greys', vmin=0, vmax=1)
    images.append(im)

    # ax.set_xticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=14)
    # ax.set_yticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.text(0.95, 0.95, lab,
            transform=ax.transAxes,
            ha='right', va='top',
            color='white', fontsize=20)

# tighten spacing between panels
plt.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.20, wspace=0.08)

cbar = fig.colorbar(images[0], ax=axes,
                    orientation='horizontal',
                    fraction=0.06, pad=0.12)

cbar.set_label(r'$\rho$', fontsize=20)
cbar.ax.tick_params(labelsize=18)
plt.savefig('figure_3.eps',dpi=300)
plt.show()