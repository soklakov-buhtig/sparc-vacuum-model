import numpy as np


def get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh=0.0):
    rho_torus = M_bar / (2 * np.pi * np.maximum(R, 1e-3) * np.maximum(z0_profile, 1e-3))
    rho_torus = np.maximum(rho_torus, 1e-25)
    return lambda_0 / rho_torus


def model_velocity_knudsen(params, R, Vbar, M_bar, R_d, z_c, M_bh, alpha):
    k_shear, lambda_0, *_ = params
    a0_base = 3600.0 * k_shear
    G_CONST = 4.30091e-6

    a_bh = (G_CONST * (M_bh * 1e9)) / np.maximum(R**2, 1e-3)
    a_newton = (Vbar**2 / R) + a_bh

    z0_profile = z_c * (1.0 + (R / R_d)**2)
    Kn = get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh)

    a_eff = (a0_base / (1.0 + Kn)) * (1.0 / np.maximum(Kn, 1e-6)**alpha) * np.exp(-Kn)
    a_total = (a_newton + np.sqrt(a_newton**2 + 4 * a_newton * a_eff)) / 2
    V_mod = np.sqrt(a_total * R)

    return V_mod
