# Exact mean-field phase diagram for self-avoiding active particles in a lattice

Source code and data accompanying the paper

> **Exact mean-field phase diagram for self-avoiding active particles in a lattice**
> *(authors, journal, year — to be filled in)*
> DOI: `<DOI-PLACEHOLDER>`

The code reproduces the mean-field master-equation (ME) results and the
linear-stability / fluctuation analysis of an active lattice gas of
self-avoiding active particles, including the spinodal surface, the
motility-induced phase-separated (MIPS) steady states, and the "active phonon"
band structure across the six Bravais lattices.

---

## Repository layout

```
active-lattice-mips/
├── simulation/        Mean-field ME solver (Fortran)
│   ├── mean_field_solver.f90   RK4 integration of the ME on an L×L square lattice
│   └── run_seeds.sh            compile & run many random initial conditions in parallel
├── stability/         Linear-stability & fluctuation analysis (Python)
│   ├── phonon_lattices.py      active-phonon dispersion for all 6 Bravais lattices
│   ├── phonons.py              command-line front-end (imports phonon_lattices)
│   ├── phonon_bands.py         square-lattice bands + closed-form curvature/sound checks
│   └── thermo.py               structure factor & entropy production (imports phonon_lattices)
└── figures/           Paper figure scripts (Python) + their input data
    ├── plot_spinodal_surface.py        → figure_1.eps   (paper Fig. 3)
    ├── plot_density_field.py           → figure_2/4.eps (paper Figs. 4 & 5)
    ├── plot_density_field_vertical.py  → figure_6.eps   (paper Fig. 6)
    ├── plot_effect_of_wt.py            → figure_3.eps   (paper "effect of w_t")
    └── data/                           steady-state .dat field files (see below)
```

> **Note:** `stability/phonons.py` and `stability/thermo.py` both
> `import phonon_lattices`, so the three files must stay together in
> `stability/` and `phonon_lattices.py` must keep its name.

---

## Model and conventions

Particles carry a director (one of the `z` nearest-neighbour directions) and
hop on a Bravais lattice with three elementary rates, all expressed in units of
the rotational diffusion rate `w_r`:

| symbol | meaning |
|--------|---------|
| `w_a`  | active (director-biased) hopping rate |
| `w_t`  | translational (unbiased) diffusion rate |
| `w_r`  | rotational diffusion rate (= 1, sets the unit) |
| `phi`  | lattice filling fraction |
| `e`    | exclusion strength (= 1, hard self-avoidance) |

Lengths are in units of the nearest-neighbour distance `a = 1`.

---

## File number ↔ paper figure number

The EPS file names do **not** match the figure numbers in the paper. Mapping:

| produced file   | paper figure | script | input |
|-----------------|--------------|--------|-------|
| `figure_1.eps`  | **Fig. 3** (spinodal surface) | `plot_spinodal_surface.py` | analytic (no data) |
| `figure_2.eps`  | **Fig. 4** (density + persistence velocity) | `plot_density_field.py` | `A_vel.dat`, `B_vel.dat` |
| `figure_3.eps`  | **effect of $w_t$** panel | `plot_effect_of_wt.py` | `F_vel.dat`, `I_vel.dat`, `G_vel.dat` |
| `figure_4.eps`  | **Fig. 5** (density + probability current) | `plot_density_field.py` *(toggled)* | `A_pc.dat`, `B_pc.dat` |
| `figure_6.eps`  | **Fig. 6** (state outside the spinodal) | `plot_density_field_vertical.py` | `C_vel.dat`, `C_pc.dat` |

There is no `figure_5.eps`: the corresponding paper figure (`escape.png`, the
mean escape time from the homogeneous state) is produced by a separate
escape-time code that is **not** part of this repository.

`plot_density_field.py` builds either `figure_2.eps` or `figure_4.eps` depending
on which `files = [...]` block is left uncommented at the top of the script
(only one may be active at a time).

---

## Steady-state data files (`figures/data/`)

These are pre-computed steady states of the mean-field ME used to draw the
figures. They are provided as static inputs; see the provenance note below.

**Header line** (commented with `#`):
`L  w_a  w_t  e  phi  branch  iters`
— lattice size, the three rates, exclusion strength, filling fraction, a
solution-branch label, and an iteration counter (`-1` = converged).

**Body columns:**
- `*_vel.dat` : `i  j  rho  v_x  v_y`  — site indices, occupation, persistence velocity
- `*_pc.dat`  : `i  j  rho  J_x  J_y`  — site indices, occupation, probability current

**Parameters per file** (all `L = 20`, `e = 1`):

| files | `w_a` | `w_t` | `phi` | used for |
|-------|------|------|------|----------|
| `A_*`, `B_*` | 20 | 0 | 0.6 | Figs. 4 & 5 (two MIPS branches) |
| `C_*`        | 20 | 0 | 0.5 | Fig. 6 (state outside spinodal) |
| `F_vel`      | 60 | 0  | 0.75 | effect-of-$w_t$ panel (a) |
| `I_vel`      | 60 | 5  | 0.75 | effect-of-$w_t$ panel (b) |
| `G_vel`      | 60 | 10 | 0.75 | effect-of-$w_t$ panel (c) |

---

## Requirements

**Python** (`stability/` and `figures/`): Python ≥ 3.8 with the packages in
[`requirements.txt`](requirements.txt) — `numpy`, `scipy`, `matplotlib`.

```bash
pip install -r requirements.txt
```

**LaTeX (figures only):** `plot_spinodal_surface.py`, `plot_effect_of_wt.py`
and `plot_density_field_vertical.py` set `text.usetex=True` with the **kpfonts**
package to match the paper typography. A working LaTeX installation including
`kpfonts` and `type1cm` is required to run them. `plot_density_field.py` and all
`stability/` scripts do **not** need LaTeX.

**Fortran** (`simulation/`): a Fortran compiler, e.g. `gfortran`.

---

## Usage

### Mean-field ME solver (Fortran)

`mean_field_solver.f90` integrates the ME on an `L×L` square lattice with RK4,
sweeping over `(phi, w_a)` and writing the Gibbs-entropy time series of each run
to `ENTROPY_wa_*_phi_*_wt_*.txt`. Edit the `PARAMETER` block at the top to change
`L`, the time step, the sweep ranges, `w_t`, and the initial-condition seed.

```bash
cd simulation
gfortran -O2 -o solver mean_field_solver.f90
./solver
```

To average over many random initial conditions in parallel (one subfolder and
seed per run):

```bash
cd simulation
bash run_seeds.sh            # uses all cores
bash run_seeds.sh 4          # cap at 4 parallel jobs
```

### Active-phonon dispersion (Python)

```bash
cd stability

# all six lattices: geometry + closed-form curvature checks + dispersion figures
python3 phonon_lattices.py

# command-line front-end for a single lattice / parameter set
python3 phonons.py --lattice square --wa 40 --wt 0 --phi 0.3
python3 phonons.py --lattice fcc --wa 60 --wt 10 --phi 0.25 --save
python3 phonons.py --lattice linear --wa 40 --phi 0.0     # ideal-gas sound

# square-lattice bands + linear-sound / curvature validation
python3 phonon_bands.py

# structure factor & entropy production (validates detailed balance at w_a=0)
python3 thermo.py
```

Running `phonon_lattices.py` and `thermo.py` prints self-consistency checks that
compare the numerical spectra against the closed-form expressions in the paper
and SI (curvature coefficients, detailed-balance, active-sound relation).

### Figures

```bash
cd figures
python3 plot_spinodal_surface.py         # figure_1.eps  (Fig. 3)   [needs LaTeX]
python3 plot_density_field.py            # figure_2.eps  (Fig. 4)
python3 plot_effect_of_wt.py             # figure_3.eps             [needs LaTeX]
python3 plot_density_field_vertical.py   # figure_6.eps  (Fig. 6)   [needs LaTeX]
```

For `figure_4.eps` (Fig. 5, probability current), open `plot_density_field.py`
and switch the active `files = [...]` block from the `*_vel.dat` pair to the
`*_pc.dat` pair.

---

## Provenance and scope

- `mean_field_solver.f90` integrates the ME and outputs the **entropy time
  series** only; it does **not** by itself write the `*_vel.dat` / `*_pc.dat`
  field files. Those are provided as pre-computed converged steady states.
- The escape-time figure (`escape.png` in the paper) is generated by a separate
  code not included here.

---

## Citation

If you use this code, please cite the paper. See [`CITATION.cff`](CITATION.cff);
the DOI is a placeholder (`<DOI-PLACEHOLDER>`) to be updated on publication.

## License

See [`LICENSE`](LICENSE). **Placeholder — choose and fill in before publishing.**
