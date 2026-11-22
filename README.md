# SPARC Vacuum Model: Continuum Mechanics Approach to Galactic Rotation Curves

Interactive research dashboard and numerical verification suite for modeling galactic rotation curves from the **SPARC** catalog (*Spitzer Photometry and Accurate Rotation Curves*, Lelli et al., 2016) without cold dark matter (CDM) halos.

## Physical Concept
Instead of postulating non-existent individual dark matter halos around galaxies, this model treats empty geometric space as passive and describes the universe as filled with a continuous elastoplastic non-fermionic energy medium. Its rheological shear response depends on the local rarefaction of visible (fermionic) matter — parameterized via an effective Knudsen number (Kn):

$$\text{Kn}(R) = \frac{\lambda_0}{\Sigma_{\text{torus}}(R)}, \quad \Sigma_{\text{torus}}(R) = \frac{M_{\text{bar}}(R)}{2\pi R \, z_0(R)} \quad \left[\frac{10^9\,M_\odot}{\text{kpc}^2}\right]$$


where $z_0(R) = z_c \left(1 + (R / R_d)^2\right)$ describes the vertical scale height flaring of the galactic disk.

The total effective acceleration combines the Newtonian baryonic field (gas, stellar disk, central bulge, and central black hole) with the shear-elastic response of the non-fermionic medium:

$$a_{\text{total}} = \frac{a_{\text{newton}} + \sqrt{a_{\text{newton}}^2 + 4 \, a_{\text{newton}} \, a_{\text{eff}}(\text{Kn})}}{2}$$

where the effective shear-elastic acceleration is derived from the trace-free stress tensor $\pi_{\mu\nu}$ (Preprint I) via a Maxwell-type viscoelastic relaxation law:
$$a_{\text{eff}}(\text{Kn}) = \frac{a_{\text{shear}}}{(1 + \text{Kn})\,\text{Kn}^{\alpha}} e^{-\text{Kn}}$$
- **Viscoelastic transition $(1 + \text{Kn})^{-1}$:** Smoothly suppresses shear stress in dense Newtonian regimes ($\text{Kn} \ll 1$).
- **Power-law elasticity $\text{Kn}^{-\alpha}$:** Captures scale-free shear stiffening of the rarefied non-fermionic substrate ($\alpha = 0.18$).
- **Cosmological shear scale $a_{\text{shear}} \propto c H_0$:** Sets the characteristic background tension of the non-fermionic medium without per-galaxy dark matter parameters.

* **High-density regime ($\text{Kn} \to 0$, Solar System & galactic cores):** Non-fermionic shear stress relaxes rapidly; the substrate behaves as an ideal superfluid and recovers standard Newtonian/Einsteinian dynamics.
* **Transition regime (galactic outskirts & dwarf disks):** Rarefaction increases $\text{Kn}$, activating shear elasticity and naturally producing asymptotically flat rotation curves.
* **Deep-space limit ($\text{Kn} \gg 1$):** Exponential damping ($e^{-\text{Kn}}$) prevents infinite halo mass divergence at large radii.

Substituting $a_{\text{newton}} = V_{\text{bar}}^2(R)/R$ into $V_{\text{total}}(R) = \sqrt{R \cdot a_{\text{total}}}$ yields the closed-form galactic rotation curve:
$$V_{\text{total}}(R) = \sqrt{ \frac{V_{\text{bar}}^2(R)}{2} + \sqrt{ \frac{V_{\text{bar}}^4(R)}{4} + R \, V_{\text{bar}}^2(R) \frac{a_{\text{shear}}}{(1 + \text{Kn})\,\text{Kn}^{\alpha}} e^{-\text{Kn}} } }$$

### Analytical Cosmological Estimate (Zero-Parameter Limit)
The characteristic shear scale can be derived directly from the cosmological horizon acceleration $a_H = c H_0$ without any galactic calibration:
$$a_{\text{shear}}^{\text{theor}} = \frac{c H_0}{2\pi} \lambda_0^{\alpha}$$
For $c = 2.998 \times 10^8\text{ m/s}$, $H_0 \approx 70\text{--}73\text{ km/(s}\cdot\text{Mpc)}$, $\lambda_0 = 5 \times 10^{-5}$, and $\alpha = 0.18$, this yields $a_{\text{shear}}^{\text{theor}} \approx (1.82\text{--}1.90) \times 10^{-11}\text{ m/s}^2$, which corresponds to $k_{\text{shear}}^{\text{theor}} \approx 0.156\text{--}0.163$ in code units ($3600\,k_{\text{shear}}$ in $[\text{(km/s)}^2/\text{kpc}]$).
This first-principles estimate agrees within ~15% with the single-galaxy calibrated value ($k_{\text{shear}} = 0.1889$ for UGC06818), confirming that $k_{\text{shear}}$ reflects the background cosmological tension of the non-fermionic medium rather than a free fitting parameter.

## Key Benchmark Results (175 SPARC Galaxies)
Unlike standard dark matter halo fits (NFW, Burkert) that require 2–3 free parameters *per galaxy* (350–500+ parameters across the catalog), this model uses **1-sample cross-validation**: the shear modulus $k_{\text{shear}}$ is calibrated on a **single reference galaxy** (`UGC06818`), frozen, and evaluated blindly across the remaining catalog without per-galaxy tuning ($\alpha = 0.18$, $\lambda_0 = 5e-5$), while vacuum_stats.py also evaluates the full leave-one-in distribution across all 175 galaxies.


| Group | Kinematic Regime | All ($N$) | All Mean MAPE | $Q=1$ ($N$) | $Q=1$ Mean MAPE |
| :--- | :--- | :---: | :---: | :---: | :---: |
| `HSB` | High-Speed Spirals ($V_{\max} \ge 150\text{ km/s}$) | 56 | 12.66% | 39 | 12.91% |
| `LSB` | Intermediate Disks ($80 \le V_{\max} < 150\text{ km/s}$) | 60 | 21.76% | 41 | 19.67% |
| `CDW` | Compact Dwarfs ($50 \le V_{\max} < 80\text{ km/s}$) | 36 | 22.29% | 18 | 19.12% |
| `LDW` | Low-Mass Dwarfs ($V_{\max} < 50\text{ km/s}$) | 22 | 40.16% | 1 | 9.40% |
| **Total** | **Out-of-Sample Catalog** | **174** | **21.27% (Med 15.00%)** | **99** | **16.80% (Med 14.30%)** |

Baryonic velocity profiles include gas, stellar disk with $\Upsilon_{\text{disk}} = 0.5$, and central stellar bulge with $\Upsilon_{\text{bul}} = 0.7$.

### All-175 Leave-One-In Cross-Validation (Calibrator Sensitivity)
To verify that the model does not rely on a lucky choice of the reference galaxy, $k_{\text{shear}}$ was calibrated sequentially on each of the 175 SPARC galaxies and tested against the remaining 174:
- **Median MAPE (25th percentile):** 15.64%
- **Median MAPE (50th percentile / center):** 18.25%
- **Median MAPE (75th percentile):** 23.60%
- **Mean MAPE (50th percentile):** 24.12%


## Repository Structure
* `streamlit_app.py` — Interactive Streamlit & Plotly web interface for real-time parameter exploration and curve visualization.
* `vacuum_physics.py` — Core continuum-mechanics equations (`get_combined_knudsen`, `model_velocity_knudsen`).
* `vacuum_calibration.py` — Single-galaxy L2 optimization routine for $k_{\text{shear}}$.
* `vacuum_stats.py` — Global MAPE evaluation and full 175-galaxy leave-one-in cross-validation distribution.
* `vacuum_groups.py` — Kinematic classification and group-level benchmarking (`HSB`, `LSB`, `CDW`, `LDW`).
* `sparc_data.py` & `sparc_parser.py` — Automated parser for `Rotmod_LTG.zip` and morphological disk scale parameters (`Table1.mrt`).

## Quick Start
```bash
git clone https://github.com/soklakov-buhtig/sparc-vacuum-model.git
cd sparc-vacuum-model
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Related Theoretical Work & Citation
This computational study accompanies the theoretical framework on flat-spacetime optical-force and continuum-vacuum representations of gravity:

> **Soklakov, A. V.** (2026). *Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations*. Zenodo. https://doi.org/10.5281/zenodo.23952161

```bibtex
@misc{soklakov2026vacuum,
  author       = {Soklakov, A. V.},
  title        = {{Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations}},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v1},
  doi          = {10.5281/zenodo.23952161},
  url          = {https://doi.org/10.5281/zenodo.23952161}
}
```
## Author & Academic Profiles
* **Alexander V. Soklakov** — Independent Researcher (Continuum Mechanics, CAE, Unified Field Theory, Emergent Gravity)
* **ORCID:** [0009-0002-4586-1800](https://orcid.org/0009-0002-4586-1800)
* **Google Scholar:** [A. V. Soklakov](https://scholar.google.com/citations?user=n4iAGDwAAAAJ)
* **Scopus Author ID:** [57209293868](https://scopus.com/authid/detail.uri?authorId=57209293868)

### Publications

* **Preprint II (Galactic Scale):** A. V. Soklakov, *Shear-Elastic Rheology of the Non-Fermionic Medium: A Zero-to-One Parameter Explanation of Galactic Rotation Curves Without Dark Matter Halos*, Zenodo (2026), [doi:10.5281/zenodo.24446953](https://doi.org)
* **Preprint I (Relativistic Scale):** A. V. Soklakov, *Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations*, Zenodo (2026), [doi:10.5281/zenodo.23952161](https://doi.org)



* **Preprint II (Galactic Scale):** A. V. Soklakov, *Shear-Elastic Rheology of the Non-Fermionic Medium: A Zero-to-One Parameter Explanation of Galactic Rotation Curves Without Dark Matter Halos*, Zenodo (2026), [doi:10.5281/zenodo.24446953]([https://doi.org](https://zenodo.org/records/24446954))
* **Preprint I (Relativistic Scale):** A. V. Soklakov, *Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations*, Zenodo (2026), [doi:10.5281/zenodo.23952161]([https://doi.org](https://doi.org/10.5281/zenodo.23952161))



### Publications


* **Preprint II (Galactic Scale):** A. V. Soklakov, *Shear-Elastic Rheology of the Non-Fermionic Medium: A Zero-to-One Parameter Explanation of Galactic Rotation Curves Without Dark Matter Halos*, Zenodo (2026), [doi:10.5281/zenodo.24446953]
* **Preprint I (Relativistic Scale):** A. V. Soklakov, *Gauge-Equivalent Optical-Force Representation of the Schwarzschild Metric via Isotropic Vacuum Fluctuations*, Zenodo (2026), [doi:10.5281/zenodo.23952161]
