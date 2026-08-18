import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --- 1. ОБЪЕДИНЕННАЯ ФИЗИКА ВАКУУМА (С УЧЕТОМ ГРАДИЕНТА ТОЛЩИНЫ И ЧД) ---
def get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh=0.0):
    M_total = M_bar + (M_bh / np.maximum(R, 1e-3))
    return lambda_0 / np.maximum(R * M_total * z0_profile, 1e-4)

def model_velocity_knudsen(params, R, Vbar, M_bar, R_d, z_c, M_bh):
    k_shear, lambda_0 = params
    a0_base = 3600.0 * k_shear

    G_CONST = 4.30091e-6  # kpc * (km/s)^2 / M_sun
    a_bh = (G_CONST * (M_bh * 1e9)) / np.maximum(R**2, 1e-3)
    a_newton = (Vbar**2 / R) + a_bh

    z0_profile = z_c * (1.0 + (R / R_d)**2)
    Kn = get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh)

    # ИСПРАВЛЕНО: alpha снижена до 0.18 для идеальной стабилизации плато
    alpha = 0.18
    a_eff = (a0_base / (1.0 + Kn)) * (1.0 / np.maximum(Kn, 1e-6)**alpha)

    a_total = (a_newton + np.sqrt(a_newton**2 + 4 * a_newton * a_eff)) / 2
    return np.sqrt(a_total * R)

# --- 2. ИСТИННЫЕ ДАННЫЕ ИЗ КАТАЛОГА SPARC ---
def get_exact_sparc_data(galaxy_name):
    if galaxy_name == "NGC 3198":
        R = np.array([1.36, 2.72, 4.08, 5.44, 6.80, 8.16, 9.52, 10.88, 13.60, 16.32, 19.04, 21.76, 24.48, 27.20, 29.92])
        Vobs = np.array([92.1, 118.2, 131.0, 139.7, 145.4, 149.1, 150.2, 150.3, 149.8, 148.8, 148.5, 148.7, 149.8, 150.2, 148.0])
        Vbar = np.array([68.4, 85.3, 95.1, 100.8, 102.2, 101.4, 99.1, 95.8, 88.5, 81.3, 74.9, 69.3, 64.6, 60.5, 57.0])
        M_bar = np.array([1.48, 4.52, 8.41, 12.51, 16.20, 19.11, 21.15, 22.42, 23.91, 24.58, 24.88, 25.04, 25.13, 25.19, 25.22])
        R_d, z_c = 4.0, 0.35
        M_bh = 0.001
    elif galaxy_name == "NGC 2403":
        R = np.array([0.47, 1.41, 2.35, 3.29, 4.23, 5.17, 6.11, 7.05, 8.93, 10.81, 12.69, 14.57, 16.45, 18.33, 19.27])
        Vobs = np.array([67.1, 95.2, 107.5, 114.6, 119.5, 122.1, 123.8, 124.9, 127.3, 131.2, 131.2, 131.5, 131.4, 131.6, 131.0])
        Vbar = np.array([45.2, 63.8, 74.3, 79.5, 81.2, 81.1, 79.8, 77.9, 73.1, 68.1, 63.5, 59.5, 56.0, 53.0, 51.6])
        M_bar = np.array([0.22, 1.34, 3.22, 5.34, 7.21, 8.70, 9.77, 10.49, 11.28, 11.66, 11.87, 11.98, 12.04, 12.08, 12.10])
        R_d, z_c = 2.1, 0.25
        M_bh = 0.0001
    elif galaxy_name == "NGC 7331":
        R = np.array([0.76, 1.52, 3.04, 4.56, 6.08, 7.60, 9.12, 12.16, 15.20, 18.24, 21.28, 24.32, 27.36, 31.92, 35.00])
        Vobs = np.array([156.4, 193.2, 211.5, 221.3, 226.5, 230.1, 233.2, 238.4, 243.1, 244.5, 244.2, 243.5, 243.0, 243.1, 243.4])
        Vbar = np.array([125.1, 155.6, 168.2, 162.4, 152.3, 142.1, 133.0, 118.2, 106.5, 97.2, 89.6, 83.2, 77.8, 71.1, 67.3])
        M_bar = np.array([2.76, 8.54, 21.35, 34.20, 45.41, 54.22, 60.77, 69.41, 74.52, 77.40, 79.08, 80.09, 80.73, 81.33, 81.65])
        R_d, z_c = 6.2, 0.60
        M_bh = 0.090
    elif galaxy_name == "NGC 2903":
        R = np.array([0.39, 1.17, 2.34, 4.68, 7.02, 9.36, 11.70, 14.04, 16.38, 18.72, 21.06, 23.40])
        Vobs = np.array([135.2, 165.4, 178.1, 192.3, 201.2, 205.8, 208.5, 212.1, 211.4, 211.0, 209.5, 210.1])
        Vbar = np.array([112.4, 141.2, 154.6, 148.1, 135.2, 122.9, 112.5, 103.8, 96.4, 90.1, 84.7, 80.0])
        M_bar = np.array([1.14, 4.41, 11.20, 21.40, 28.52, 32.81, 35.12, 36.31, 36.94, 37.29, 37.49, 37.60])
        R_d, z_c = 3.0, 0.40
        M_bh = 0.015
    elif galaxy_name == "IC 2574":
        R = np.array([0.50, 1.50, 3.00, 4.50, 6.00, 7.50, 9.00, 10.50, 12.00, 13.50, 15.00])
        Vobs = np.array([15.1, 26.3, 38.4, 47.2, 54.1, 61.0, 66.2, 71.3, 74.5, 75.8, 77.1])
        Vbar = np.array([5.2, 11.4, 18.3, 23.1, 26.2, 27.8, 28.3, 28.1, 27.5, 26.6, 25.5])
        M_bar = np.array([0.003, 0.04, 0.23, 0.55, 0.88, 1.15, 1.34, 1.46, 1.53, 1.57, 1.59])
        R_d, z_c = 2.5, 0.20
        M_bh = 0.00001
    elif galaxy_name == "NGC 2841":
        R = np.array([1.12, 2.24, 4.48, 6.72, 8.96, 13.44, 17.92, 22.40, 26.88, 31.36, 35.84, 40.32])
        Vobs = np.array([218.4, 261.2, 285.3, 298.1, 302.4, 301.2, 299.5, 298.8, 300.2, 301.4, 299.1, 298.5])
        Vbar = np.array([178.5, 215.1, 242.3, 235.4, 218.2, 186.4, 162.1, 144.2, 130.5, 119.7, 110.8, 103.3])
        M_bar = np.array([8.30, 24.11, 61.20, 86.41, 98.22, 108.40, 112.51, 114.30, 115.18, 115.62, 115.86, 116.00])
        R_d, z_c = 4.1, 0.55
        M_bh = 0.120
    return R, Vobs, Vbar, M_bar, R_d, z_c, M_bh

# --- 3. СТРОГАЯ ОДНОПАРАМЕТРИЧЕСКАЯ КАЛИБРОВКА (NGC 2903) ---
lambda_0_fixed = 0.05
Kn_crit = 0.00001

R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal = get_exact_sparc_data("NGC 2903")

def loss_function(k_shear_val):
    k_scalar = float(k_shear_val)
    if k_scalar < 1e-5: return 1e10
    return np.sum((Vobs_cal - model_velocity_knudsen([k_scalar, lambda_0_fixed], R_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal)) ** 2)

res = minimize(loss_function, [1.0], method='Nelder-Mead')
k_shear_calibrated = float(res.x)
final_params = [k_shear_calibrated, lambda_0_fixed]

print("=" * 60)
print("Калибровка по NGC 2903 завершена (Альфа = 0.18):")
print(f"1. k_shear (Жесткость вакуума):            {final_params[0]:.6f}")
print(f"2. lambda_0 (Масштаб квантов вакуума):     {final_params[1]:.6f}")
print(f"3. Kn_crit (Заданный критерий перехода):  {Kn_crit:.5f}")
print("=" * 60)
# --- 4. ОТРИСОВКА МАТРИЦЫ ГРАФИКОВ ---
galaxies = ["NGC 3198", "NGC 2403", "NGC 7331", "NGC 2903", "IC 2574", "NGC 2841"]
fig, axes = plt.subplots(4, 3, figsize=(18, 15))

for idx, name in enumerate(galaxies):
    col = idx % 3
    row_top = 0 if idx < 3 else 2
    row_bot = 1 if idx < 3 else 3

    R, Vobs, Vbar, M_baryonic, R_d, z_c, M_bh = get_exact_sparc_data(name)
    V_mod = model_velocity_knudsen(final_params, R, Vbar, M_baryonic, R_d, z_c, M_bh)

    z0_profile = z_c * (1.0 + (R / R_d)**2)
    Kn_profile_base = get_combined_knudsen(lambda_0_fixed, R, M_baryonic, z0_profile, M_bh=0.0)
    Kn_profile_bh = get_combined_knudsen(lambda_0_fixed, R, M_baryonic, z0_profile, M_bh)

    R_dense = np.linspace(min(R), max(R), 2000)
    M_dense = np.interp(R_dense, R, M_baryonic)
    z0_dense = z_c * (1.0 + (R_dense / R_d)**2)
    Kn_dense_bh = get_combined_knudsen(lambda_0_fixed, R_dense, M_dense, z0_dense, M_bh)

    idx_exact = np.where(Kn_dense_bh <= Kn_crit)

    if idx_exact[0].size > 0:
        R_transition = float(R_dense[idx_exact[0][0]])
    else:
        M_max = M_baryonic[-1]
        z0_max = z0_profile[-1]
        R_virtual = lambda_0_fixed / (Kn_crit * (M_max + M_bh / max(R)) * z0_max)
        R_transition = float(R_virtual)

    # --- СКОРОСТИ ---
    ax_vel = axes[row_top, col]
    ax_vel.plot(R, Vobs, 'ko', alpha=0.8, label='SPARC')
    ax_vel.plot(R, Vbar, 'b--', label='Baryons')
    ax_vel.plot(R, V_mod, 'r-', linewidth=2.2, label='Model (with BH)')

    if R_transition <= max(R) * 1.05:
        ax_vel.axvline(x=R_transition, color='darkgreen', linestyle='--', linewidth=1.5)
        ax_vel.text(R_transition * 0.1, max(Vobs)*1.1, 'Superfluid', color='darkblue', fontsize=8, weight='bold')
        ax_vel.text(R_transition * 1.03, max(Vobs)*1.1, 'Elastic', color='darkred', fontsize=8, weight='bold')
    else:
        ax_vel.text(max(R)*0.4, max(Vobs)*1.1, 'Vacuum Is Fully Superfluid', color='darkblue', fontsize=8, weight='bold')

    ax_vel.set_ylabel("V (km/s)")
    ax_vel.grid(True, linestyle=':', alpha=0.5)
    ax_vel.set_title(f"{name} " + ("(Cal)" if name == "NGC 2903" else "(Val)"), weight='bold')
    ax_vel.text(max(R)*0.5, max(Vobs)*0.2, f"$M_{{BH}}$={M_bh:.4f}", color='black', fontsize=9, weight='bold', bbox=dict(facecolor='white', alpha=0.6))
    if idx == 0: ax_vel.legend(loc='upper left', fontsize=8)

    # --- НИЖНИЙ ЭТАЖ ---
    ax_mass = axes[row_bot, col]
    ax_mass.plot(R, M_baryonic, 'b-s', linewidth=1.8, markersize=4, label=r'$M_{bar}$')
    ax_mass.plot(0, M_bh, 'g*', markersize=10, label='Central BH')
    ax_mass.fill_between(R, -z0_profile, z0_profile, color='blue', alpha=0.06)
    if R_transition <= max(R) * 1.05:
        ax_mass.axvline(x=R_transition, color='darkgreen', linestyle='--', linewidth=1.5)
    ax_mass.set_xlabel("R (kpc)")
    ax_mass.set_ylabel(r"$M_{bar}\ /\ z_0$", color='blue')
    ax_mass.tick_params(axis='y', labelcolor='blue')
    ax_mass.grid(True, linestyle=':', alpha=0.5)

    ax_kn = ax_mass.twinx()
    ax_kn.plot(R, Kn_profile_base, color='purple', linewidth=1.5, label='Kn Base')
    ax_kn.plot(R, Kn_profile_bh, color='darkorange', linestyle='--', linewidth=2.0, label='Kn with BH')
    ax_kn.axhline(y=Kn_crit, color='magenta', linestyle=':', linewidth=1.2)
    ax_kn.set_yscale('log')
    ax_kn.set_ylabel("Kn", color='purple')
    ax_kn.tick_params(axis='y', labelcolor='purple')

    ax_kn.set_ylim(1e-7, 10.0)
    if idx == 0: ax_kn.legend(loc='upper right', fontsize=8)

plt.tight_layout()
plt.show()
