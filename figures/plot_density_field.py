# -*- coding: utf-8 -*-
"""
plot_density_field.py -- side-by-side density + arrow-field panels.

Occupation (grayscale) with an overlaid arrow field whose thickness encodes the
local magnitude. Select which paper figure to build by uncommenting ONE `files`
block below (only one may be active at a time):

    figure_2.eps (paper Fig. 4): A_vel.dat, B_vel.dat   persistence velocity
    figure_4.eps (paper Fig. 5): A_pc.dat,  B_pc.dat    probability current

Data files are read from ./data/ . No LaTeX is required for this script.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# ***** Fig. 4 (file figure_2.eps): persistence velocity
files = ['A_vel.dat', 'B_vel.dat']

# ***** Fig. 5 (file figure_4.eps): probability current
# files = ['A_pc.dat', 'B_pc.dat']

fig = plt.figure(figsize=(8, 4))
gs = GridSpec(2, 2, height_ratios=[15, 1], hspace=0.4, wspace=0.3)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

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
                   cmap='Greys',   # inverted Greys
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
    ax.set_xticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=0)
    ax.set_yticks([0,4,9,14,19], ['1','5','10','15','20'], fontsize=0)

# Labels
# ax1.set_xlabel('Lattice $i$', fontsize=14)
# ax1.set_ylabel('Lattice $j$', fontsize=14)
# ax2.set_xlabel('Lattice $i$', fontsize=14)
# ax2.set_ylabel('Lattice $j$', fontsize=14)

ax1.annotate('(a)', xy=(0.5, 17.5), fontsize=18)
ax2.annotate('(b)', xy=(0.5, 17.5), fontsize=18)

# Colorbar
plt.subplots_adjust(top=0.965, bottom=0.210)
cbar_ax = fig.add_axes([0.3, 0.16, 0.4, 0.03])

cbar = fig.colorbar(images[0],
                    cax=cbar_ax,
                    orientation='horizontal')

# Tick labels (numbers)
cbar.ax.tick_params(labelsize=14)

# Label (ρ)
cbar.set_label(r'$\rho$', fontsize=14)


plt.savefig('densities_vector.png', dpi=300)
plt.savefig('densities_vector.eps', dpi=300)
plt.show()