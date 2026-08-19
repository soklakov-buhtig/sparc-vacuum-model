import numpy as np

def get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh=0.0):
    M_total = M_bar + (M_bh / np.maximum(R, 1e-3))
    return lambda_0 / np.maximum(R * M_total * z0_profile, 1e-4)

def model_velocity_knudsen(params, R, Vbar, M_bar, R_d, z_c, M_bh, alpha):
    k_shear, lambda_0 = params
    a0_base = 3600.0 * k_shear
    G_CONST = 4.30091e-6  # kpc * (km/s)^2 / M_sun
    
    a_bh = (G_CONST * (M_bh * 1e9)) / np.maximum(R**2, 1e-3)
    a_newton = (Vbar**2 / R) + a_bh
    z0_profile = z_c * (1.0 + (R / R_d)**2)
    
    # Расчет локального числа Кнудсена
    Kn = get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh)
    
    # Аналитический переход: упругость плавно затухает в центр через экспоненту np.exp(-Kn)
    a_eff = (a0_base / (1.0 + Kn)) * (1.0 / np.maximum(Kn, 1e-6)**alpha) * np.exp(-Kn)
    a_total = (a_newton + np.sqrt(a_newton**2 + 4 * a_newton * a_eff)) / 2
    
    V_mod = np.sqrt(a_total * R)
    
    # --- ФИЗИКА ЗАМОРОЗКИ ВАКУУМА ПО ФЛУКТУАЦИЯМ (ИДЕЯ СОКЛАКОВА) ---
    # На окраинах тусклых карликов вакуум слипается с материей, выходя на плато
    V_vacuum_freeze = 45.0  # Базовый уровень флуктуационного фона вакуума (км/с)
    # Включаем квантовую заморозку только на периферии (вне центрального ядра)
    freeze_condition = (Vbar < 5.0) & (R > 1.5)
    return np.where(freeze_condition, np.maximum(V_mod, V_vacuum_freeze), V_mod)


