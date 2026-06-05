# -*- coding: utf-8 -*-
"""
plot_spinodal_surface.py -- paper Figure 3 (file: figure_1.eps).

3D wireframe of the square-lattice spinodal surface
    f(phi, wa) = wa^2 (1-phi)(2 phi - 1)/4 - wa/4 ,
the boundary of the region where the homogeneous mean-field state is locally
unstable. Requires a LaTeX installation with the kpfonts package.
"""
import numpy as np
import matplotlib.pyplot as plt
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
# Define the function
def f(phi, wa):
    return wa**2 * (1 - phi) * (2 * phi - 1) / 4 - wa / 4

# Generate data
wa = np.linspace(0, 100, 400)
phi = np.linspace(0, 1, 400)
PHI, WA = np.meshgrid(phi, wa)
Z = f(PHI, WA)

# Mask values where Z < 0
Z[Z < 0.1] = np.nan



# Create figure and 3D axis
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')

phi1, wa1, wt1 = 0.6, 10, 0
phi2, wa2, wt2 = 0.6, 20, 0
# phi3, wa3, wt3 = 0.75, 60, 0
# phi4, wa4, wt4 = 0.75, 60, 5
# phi5, wa5, wt5 = 0.75, 60, 10



# Plot meshgrid
ax.plot_wireframe(PHI, WA, Z, color='blue', linewidth=0.5)

# Black point (fully visible)
# ax.scatter(phi1, wa1, wt1, color='black', s=30, depthshade=True)
# Red point (partially eclipsed)
# ax.scatter(phi2, wa2, wt2, color='red', s=30, alpha=0.6, depthshade=True)
# ax.scatter(phi3, wa3, wt3, color='gold', s=30, alpha=1, depthshade=True)
# ax.scatter(phi4, wa4, wt4, color='gold', s=30, alpha=1, depthshade=True)
# ax.scatter(phi5, wa5, wt5, color='gold', s=30, alpha=1, depthshade=True)


# Labels and title
ax.set_xlabel(r'$\phi$', fontsize=18,labelpad=10)
ax.set_ylabel(r'$w_a~~[w_r]$', fontsize=18,labelpad=10)
ax.set_zlabel(r'$w_t~~[w_r]$', fontsize=18,labelpad=10)
plt.subplots_adjust(left=0)
# Show plot
plt.savefig('figure_1.eps')
plt.show()
