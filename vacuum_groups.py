# vacuum_groups.py
import numpy as np
import streamlit as st

def calculate_and_display_group_stats(sparc_database, calibration_galaxy, final_params, alpha, get_exact_sparc_data, model_velocity_knudsen, calibrate_vacuum_parameters):
    """
    Разделяет галактики на 4 группы, находит для каждой группы лучшего
    локального калибратора и выводит индивидуальную статистику MAPE.
    """
    galaxies_list = list(sparc_database.keys())
    
    # 1. Шаг предварительной классификации всей базы
    group_assignments = {}
    for g_name in galaxies_list:
        _, Vobs, _, _, _, _, _ = get_exact_sparc_data(g_name)
        v_max = np.max(Vobs)
        
        if v_max >= 150: group_key = "HSB"
        elif 80 <= v_max < 150: group_key = "LSB"
        elif 50 <= v_max < 80: group_key = "CDW"
        else: group_key = "LDW"
        
        group_assignments[g_name] = group_key

    # Словарь для хранения финальных массивов MAPE текущей калибровки
    current_mapes_by_group = {"HSB": [], "LSB": [], "CDW": [], "LDW": []}
    
    # Сначала собираем стандартные MAPE по группам для текущей выбранной в UI галактики
    for g_name, group_key in group_assignments.items():
        if g_name == calibration_galaxy:
            continue
        R, Vobs, Vbar, Md, Rd, zc, Mbh = get_exact_sparc_data(g_name)
        Vmod = model_velocity_knudsen(final_params, R, Vbar, Md, Rd, zc, Mbh, alpha)
        mape = np.mean(np.abs((Vobs - Vmod) / Vobs)) * 100
        current_mapes_by_group[group_key].append(mape)

    # 2. ШАГ ЧЕТЫРЕХКРАТНОГО ПЕРЕБОРА ПОБЕДИТЕЛЕЙ ГРУПП
    best_galaxies_by_group = {}
    
    for current_group in ["HSB", "LSB", "CDW", "LDW"]:
        # Выделяем только галактики, принадлежащие текущей группе
        group_galaxies = [name for name, grp in group_assignments.items() if grp == current_group]
        
        best_group_gal = "Нет данных"
        lowest_group_mean = float('inf')
        
        # Запускаем перебор калибраторов строго внутри этой группы
        for test_cal_gal in group_galaxies:
            R_c, Vobs_c, Vbar_c, M_c, Rd_c, zc_c, Mbh_c = get_exact_sparc_data(test_cal_gal)
            k_shear_test = calibrate_vacuum_parameters(R_c, Vobs_c, Vbar_c, M_c, Rd_c, zc_c, Mbh_c, alpha, final_params[1])
            test_params = [k_shear_test, final_params[1]]
            
            # Считаем среднее MAPE по галактикам этой же группы (исключая тестовую калибровочную)
            test_group_mapes = []
            for g_name in group_galaxies:
                if g_name == test_cal_gal:
                    continue
                R, Vobs, Vbar, Md, Rd, zc, Mbh = get_exact_sparc_data(g_name)
                Vmod = model_velocity_knudsen(test_params, R, Vbar, Md, Rd, zc, Mbh, alpha)
                mape = np.mean(np.abs((Vobs - Vmod) / Vobs)) * 100
                test_group_mapes.append(mape)
                
            if test_group_mapes:
                group_mean = np.mean(test_group_mapes)
                if group_mean < lowest_global_mean if 'lowest_global_mean' in locals() else group_mean < lowest_group_mean:
                    lowest_group_mean = group_mean
                    best_group_gal = test_cal_gal
                    
        best_galaxies_by_group[current_group] = best_group_gal

    # 3. Вывод результатов на боковую панель Streamlit
    st.sidebar.markdown("---")
    st.sidebar.subheader("📊 Анализ по группам:")
    
    for key in ["HSB", "LSB", "CDW", "LDW"]:
        values = current_mapes_by_group[key]
        n_count = len(values)
        if n_count > 0:
            mean_val = np.mean(values)
            st.sidebar.markdown(f"**{key}** (n={n_count}): **{mean_val:.2f}%**")
            st.sidebar.markdown(f"<p style='font-size:11px; color:#4F8BF9; margin-top:-12px; margin-bottom:8px;'>🏆 Лучшая в группе: <b>{best_galaxies_by_group[key]}</b></p>", unsafe_allow_html=True)
