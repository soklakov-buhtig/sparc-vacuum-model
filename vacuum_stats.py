# vacuum_stats.py
import numpy as np

def calculate_global_mape_stats(sparc_database, calibration_galaxy, final_params, alpha, get_exact_sparc_data, model_velocity_knudsen):
    """
    Вычисляет глобальные статистические показатели MAPE по всей выборке из 174 галактик,
    исключая текущую опорную галактику калибровки.
    """
    galaxies_list = list(sparc_database.keys())
    all_mapes = {}

    for g_name in galaxies_list:
        if g_name == calibration_galaxy:
            continue
            
        # Извлекаем экспериментальные данные
        R_g, Vobs_g, Vbar_g, M_bar_g, R_d_g, z_c_g, M_bh_g = get_exact_sparc_data(g_name)
        
        # Рассчитываем теоретическую скорость модели
        V_mod_g = model_velocity_knudsen(final_params, R_g, Vbar_g, M_bar_g, R_d_g, z_c_g, M_bh_g, alpha)
        
        # Фиксируем локальный MAPE
        mape_g = np.mean(np.abs((Vobs_g - V_mod_g) / Vobs_g)) * 100
        all_mapes[g_name] = mape_g

    if len(all_mapes) > 0:
        mapes_values = list(all_mapes.values())
        
        max_mape = np.max(mapes_values)
        min_mape = np.min(mapes_values)
        
        # Вытаскиваем имена полярных систем
        max_galaxy = [name for name, val in all_mapes.items() if val == max_mape][0]
        min_galaxy = [name for name, val in all_mapes.items() if val == min_mape][0]
        
        return {
            "mean": np.mean(mapes_values),
            "max_value": max_mape,
            "max_name": max_galaxy,
            "min_value": min_mape,
            "min_name": min_galaxy,
            "std": np.std(mapes_values)
        }
    
    return None
