# vacuum_groups.py
import numpy as np
import streamlit as st

def calculate_and_display_group_stats(sparc_database, calibration_galaxy, final_params, alpha, get_exact_sparc_data, model_velocity_knudsen):
    """
    Классифицирует галактики по физическим группам и выводит детальную статистику MAPE.
    """
    groups = {"HSB": [], "LSB": [], "CDW": [], "LDW": []}
    
    for g_name in sparc_database.keys():
        if g_name == calibration_galaxy: continue
        
        R, Vobs, Vbar, Md, Rd, zc, Mbh = get_exact_sparc_data(g_name)
        Vmod = model_velocity_knudsen(final_params, R, Vbar, Md, Rd, zc, Mbh, alpha)
        
        mape = np.mean(np.abs((Vobs - Vmod) / Vobs)) * 100
        v_max = np.max(Vobs)
        
        # Инженерная классификация
        if v_max >= 150: groups["HSB"].append(mape)
        elif 80 <= v_max < 150: groups["LSB"].append(mape)
        elif 50 <= v_max < 80: groups["CDW"].append(mape)
        else: groups["LDW"].append(mape)

    st.sidebar.markdown("---")
    st.sidebar.subheader("📊 Анализ по группам:")
    for key, values in groups.items():
        if values:
            st.sidebar.markdown(f"**{key}** (n={len(values)}): {np.mean(values):.2f}%")
