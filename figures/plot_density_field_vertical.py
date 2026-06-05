# -*- coding: utf-8 -*-
"""
plot_density_field_vertical.py -- paper Figure 6 (file: figure_6.eps).

Two stacked panels for the locally stable inhomogeneous state OUTSIDE the
spinodal surface (w_a = 20 w_r, w_t = 0, phi = 0.5):
    (a) C_vel.dat  density + persistence velocity
    (b) C_pc.dat   density + probability current
Data files are read from ./data/ . Requires a LaTeX install with kpfonts.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize

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
# ***** Fig. 6 (file figure_6.eps): density + persistence velocity (a)
#       and density + probability current (b) for the same state.
files = ['C_vel.dat', 'C_pc.dat']

fig = plt.figure(figsize=(5, 8))
gs = GridSpec(2, 1, hspace=0.05)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])

axes = [ax1, ax2]
images = []

for idx, f in enumerate(files):

    # -----------------------------
    # LOAD DATA: i, j, rho, vx, vy
    # -----------------------------
    data = np.loadtxt(os.path.join(DATA_DIR, f))

    # Original coordinates: 1..20
    i = data[:, 0].astype(int) - 1
    j = data[:, 1].astype(int) - 1
    rho = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]

    # Lattice size (should now be 20)
    Nx = i.max() + 1
    Ny = j.max() + 1

    # -----------------------------------
    # BUILD 2D MATRICES FROM POINT DATA
    # -----------------------------------
    density = np.zeros((Ny, Nx))
    vx_matrix = np.zeros((Ny, Nx))
    vy_matrix = np.zeros((Ny, Nx))

    density[j, i] = rho
    vx_matrix[j, i] = vx
    vy_matrix[j, i] = vy

    # -----------------------------------
    # VECTOR FIELD PROPERTIES
    # -----------------------------------
    modulus = np.sqrt(vx_matrix**2 + vy_matrix**2)

    # Avoid division by zero
    u_dir = np.zeros_like(vx_matrix)
    v_dir = np.zeros_like(vy_matrix)

    mask_nonzero = modulus > 0
    u_dir[mask_nonzero] = vx_matrix[mask_nonzero] / modulus[mask_nonzero]
    v_dir[mask_nonzero] = vy_matrix[mask_nonzero] / modulus[mask_nonzero]

    # -----------------------------------
    ax = axes[idx]

    # Inverted Greys heatmap
    vmax = np.max(density)
    im = ax.imshow(density,
                   cmap='Greys',
                   vmin=0,
                   vmax=1,
                   origin='lower')

    images.append(im)

    # Grid for quiver
    X, Y = np.meshgrid(np.arange(Nx), np.arange(Ny))

    # -----------------------------------
    # INVERTED COLORMAP ARROW LOGIC
    # -----------------------------------
    norm = Normalize(vmin=0, vmax=vmax)
    density_norm = norm(density)

    threshold = 0.5
    dark_regions = density_norm < threshold
    light_regions = density_norm >= threshold

    # Flatten arrays
    X_f = X.flatten()
    Y_f = Y.flatten()
    U_f = u_dir.flatten()
    V_f = v_dir.flatten()
    modulus_f = modulus.flatten()

    # Width scaling ∝ magnitude
    max_mod = np.max(modulus_f) if np.max(modulus_f) > 0 else 1.0
    desired_widths = 0.00 + 0.02 * (modulus_f / max_mod)

    # Bin widths
    n_bins = 30
    bins = np.linspace(np.min(desired_widths),
                       np.max(desired_widths),
                       n_bins)

    # ---- DARK REGIONS → BLACK ARROWS ----
    for k in range(n_bins - 1):
        mask_width = (desired_widths >= bins[k]) & (desired_widths < bins[k+1])
        mask = mask_width & dark_regions.flatten()

        if np.any(mask):
            ax.quiver(X_f[mask], Y_f[mask],
                      U_f[mask], V_f[mask],
                      scale=25,
                      width=bins[k],
                      color='black',
                      pivot='mid',
                      alpha=1)

    # ---- LIGHT REGIONS → WHITE ARROWS ----
    for k in range(n_bins - 1):
        mask_width = (desired_widths >= bins[k]) & (desired_widths < bins[k+1])
        mask = mask_width & light_regions.flatten()

        if np.any(mask):
            ax.quiver(X_f[mask], Y_f[mask],
                      U_f[mask], V_f[mask],
                      scale=25,
                      width=bins[k],
                      color='white',
                      pivot='mid',
                      alpha=0.9)

    # Last bin edge
    mask_last = desired_widths >= bins[-1]

    mask = mask_last & dark_regions.flatten()
    if np.any(mask):
        ax.quiver(X_f[mask], Y_f[mask],
                  U_f[mask], V_f[mask],
                  scale=25,
                  width=bins[-1],
                  color='black',
                  pivot='mid',
                  alpha=1)

    mask = mask_last & light_regions.flatten()
    if np.any(mask):
        ax.quiver(X_f[mask], Y_f[mask],
                  U_f[mask], V_f[mask],
                  scale=25,
                  width=bins[-1],
                  color='white',
                  pivot='mid',
                  alpha=1)

    # -----------------------------------
    # AXIS FORMATTING
    # -----------------------------------
    # ax.set_xticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=0)
    # ax.set_yticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=0)
    ax.set_xticks([])
    ax.set_yticks([])

ax1.annotate('(a)', xy=(0.5, 17.5), fontsize=24)
ax2.annotate('(b)', xy=(0.5, 17.5), fontsize=24)

# Colorbar
plt.subplots_adjust(top=0.975, bottom=0.120)
cbar_ax = fig.add_axes([0.2, 0.09, 0.6, 0.02])

cbar = fig.colorbar(images[0],
                    cax=cbar_ax,
                    orientation='horizontal')

cbar.ax.tick_params(labelsize=18)
cbar.set_label(r'$\rho$', fontsize=20)

# plt.savefig('densities_vector.png',1 dpi=300)
plt.savefig('figure_6.eps', dpi=300)
plt.show()