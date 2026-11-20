# SPARC Vacuum Model: Continuum Mechanics Approach to Galactic Rotation Curves

Interactive research dashboard and numerical verification suite for modeling galactic rotation curves from the **SPARC** catalog (*Spitzer Photometry and Accurate Rotation Curves*, Lelli et al., 2016) without cold dark matter (CDM) halos.

## Physical Concept
Instead of introducing non-baryonic dark matter particles, this model treats the physical vacuum as a polarizable, shear-elastic continuous medium whose rheological response depends on the local rarefaction of baryonic matter — parameterized via a geometric **Knudsen number ($\text{Kn}$)**:

$$\text{Kn}(R) = \frac{\lambda_0}{\rho_{\text{torus}}(R)}, \quad \rho_{\text{torus}}(R) = \frac{M_{\text{bar}}(R)}{2\pi R \, z_0(R)}$$

where $z_0(R) = z_c \left(1 + (R / R_d)^2\right)$ describes the vertical scale height flaring of the galactic disk.

The total effective acceleration combines the Newtonian baryonic field (gas, stellar disk, central bulge, and central black hole) with the shear-elastic response of the vacuum medium:

$$a_{\text{total}} = \frac{a_{\text{newton}} + \sqrt{a_{\text{newton}}^2 + 4 \, a_{\text{newton}} \, a_{\text{eff}}(\text{Kn})}}{2}$$

* **High-density regime ($\text{Kn} \to 0$, Solar System & galactic cores):** Vacuum shear stress relaxes rapidly; the medium behaves as an ideal superfluid and recovers standard Newtonian/Einsteinian dynamics.
* **Transition regime (galactic outskirts & dwarf disks):** Rarefaction increases $\text{Kn}$, activating shear elasticity and naturally producing asymptotically flat rotation curves.
* **Deep-space limit ($\text{Kn} \gg 1$):** Exponential damping ($e^{-\text{Kn}}$) prevents infinite halo mass divergence at large radii.

## Key Benchmark Results (175 SPARC Galaxies)
Unlike standard dark matter halo fits (NFW, Burkert) that require 2–3 free parameters *per galaxy* (350–500+ parameters across the catalog), this model uses **1-sample cross-validation**: the shear modulus $k_{\text{shear}}$ is calibrated on a **single reference galaxy** (`UGC00891`), frozen, and evaluated across all remaining **174 galaxies** without per-galaxy tuning ($\alpha = 0.18$, $\lambda_0 = 0.05$).

| Kinematic Group ($V_{\max}$) | Sample Size ($n$) | Mean MAPE | Best Group Reference |
| :--- | :---: | :---: | :---: |
| **HSB** ($V_{\max} \ge 150\text{ km/s}$) | 56 | **12.66%** | `NGC2998` |
| **LSB** ($80 \le V_{\max} < 150\text{ km/s}$) | 59 | **20.69%** | `NGC2403` |
| **CDW** ($50 \le V_{\max} < 80\text{ km/s}$) | 37 | **21.48%** | `UGC05750` |
| **LDW** ($V_{\max} < 50\text{ km/s}$) | 22 | **40.13%** | `D564-8` |
| **Full SPARC Sample (excl. calibrator)** | **174** | **21.26%** | **`UGC00891`** |

*(Baryonic velocity profiles include gas, stellar disk with $\Upsilon_{\text{disk}} = 0.5$, and central stellar bulge with $\Upsilon_{\text{bul}} = 0.7$).*

## Repository Structure
* `streamlit_app.py` — Interactive Streamlit & Plotly web interface for real-time parameter exploration and curve visualization.
* `vacuum_physics.py` — Core continuum-mechanics equations (`get_combined_knudsen`, `model_velocity_knudsen`).
* `vacuum_calibration.py` — Single-galaxy L2 optimization routine for $k_{\text{shear}}$.
* `vacuum_stats.py` — Global MAPE evaluation and automated leave-one-in calibrator search across all 175 galaxies.
* `vacuum_groups.py` — Kinematic classification and group-level benchmarking (`HSB`, `LSB`, `CDW`, `LDW`).
* `sparc_data.py` & `sparc_parser.py` — Automated parser for `Rotmod_LTG.zip` and morphological disk scale parameters (`Table1.mrt`).

## Quick Start
```bash
git clone https://github.com
cd sparc-vacuum-model
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Related Theoretical Work & Citation
This computational study accompanies the theoretical framework on flat-spacetime optical-force and continuum-vacuum representations of gravity:

> **Soklakov, A. V.** (2026). *Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations*. Zenodo. https://doi.org

```bibtex
@misc{soklakov2026vacuum,
  author       = {Soklakov, A. V.},
  title        = {{Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations}},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v1},
  doi          = {10.5281/zenodo.23952161},
  url          = {https://doi.org}
}
```
## Author & Academic Profiles
* **Alexander V. Soklakov** — Independent Researcher (Continuum Mechanics, CAE, Unified Field Theory, Emergent Gravity)
* **ORCID:** [0009-0002-4586-1800](https://orcid.org)
* **Google Scholar:** [A. V. Soklakov](https://google.com)
* **Scopus Author ID:**(https://scopus.com)
* **Zenodo Preprint:** [doi:10.5281/zenodo.23952161](https://doi.org)

