import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import minimize
from sparc_data import SPARC_DATABASE
from vacuum_physics import model_velocity_knudsen, get_combined_knudsen
from vacuum_calibration import calibrate_vacuum_parameters
# Добавить к вашим импортам в начале файла
from vacuum_stats import calculate_global_mape_stats, find_best_calibration_galaxy
from vacuum_groups import calculate_and_display_group_stats

# Настройка страницы веб-интерфейса
st.set_page_config(page_title="SPARC Vacuum Model", layout="wide")
st.title("🛸 Интерактивная модель упругого вакуума Соклакова")
st.markdown("### Верификация кривых вращения каталога SPARC на базе механики сплошных сред")


# --- 2. ИСТИННЫЕ ДАННЫЕ ИЗ КАТАЛОГА SPARC ---
def get_exact_sparc_data(galaxy_name):
    g = SPARC_DATABASE[galaxy_name]
    
    # Считаем суммарный барионный профиль скорости по Ньютону из чистых компонентов
    vbul = g.get("Vbul", 0.0)
    vbar = np.sqrt(g["Vgas"]**2 + 0.5 * g["Vdisk"]**2 + 0.7 * vbul**2)
#   vbar = np.sqrt(g["Vgas"]**2 + 0.5 * g["Vdisk"]**2) #старая версия
    
    
    # Считаем эталонный профиль барионной массы
    G_const = 4.30091e-6
#   m_bar = (vbar**2) * (g["R"] * 1000) / G_const / 1e9
    m_bar = (vbar**2) * (g["R"]) / G_const / 1e9
#   m_bar = np.maximum(m_bar, 1e-3) # Защита центра от нуля
    m_bar = np.maximum(m_bar, 1e-6) # Защита центра от нуля
    
    return g["R"], g["Vobs"], vbar, m_bar, g["R_d"], g["z_c"], g["M_bh"]


# --- 3. ИНТЕРАКТИВНАЯ БОКОВАЯ ПАНЕЛЬ УПРАВЛЕНИЯ ---
st.sidebar.header("🎛️ Глобальные параметры среды")

alpha = st.sidebar.slider(
    "Степень стабилизации плато (alpha):", 
    min_value=0.05, max_value=0.50, value=0.18, step=0.01
)

lambda_0_fixed = st.sidebar.slider(
    "Масштаб квантов вакуума (lambda_0):", 
    min_value=1e-5, max_value=2e-4, value=5e-5, step=1e-5, format="%.5f"
#   min_value=0.01, max_value=0.20, value=0.05, step=0.01
)

# --- 4. ДИНАМИЧЕСКИЙ ВЫБОР ГАЛАКТИКИ ДЛЯ КАЛИБРОВКИ ---
st.sidebar.markdown("---")
st.sidebar.subheader("🎯 Настройка калибровки")

# Перебираем базу и находим имя лучшей галактики для текущих ползунков alpha и lambda_0
dummy_params = [0.0, lambda_0_fixed]  # k_shear откалибруется внутри
best_galaxy_found = find_best_calibration_galaxy(
    SPARC_DATABASE, dummy_params, alpha, 
    get_exact_sparc_data, model_velocity_knudsen, calibrate_vacuum_parameters
)

# Передаем найденное имя как index по умолчанию в ваш селектбокс
galaxies_keys = list(SPARC_DATABASE.keys())
default_index = galaxies_keys.index(best_galaxy_found) if best_galaxy_found in galaxies_keys else 0

calibration_galaxy = st.sidebar.selectbox(
    "Опорная галактика для калибровки:",
    galaxies_keys,
    index=default_index
)



R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal = get_exact_sparc_data(calibration_galaxy)

k_shear_calibrated = calibrate_vacuum_parameters(R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal, alpha, lambda_0_fixed)
final_params = [k_shear_calibrated, lambda_0_fixed]



st.sidebar.subheader("📈 Результаты калибровки:")
st.sidebar.code(f"k_shear: {k_shear_calibrated:.6f}\nlambda_0: {lambda_0_fixed:.6f}")

# --- ЕДИНАЯ ФУНКЦИЯ ОТРИСОВКИ БЛОКА СТАТИСТИКИ И ГРУПП ---
def render_stats_block(title_prefix, quality_filter=None, show_best_cal=False):
    stats = calculate_global_mape_stats(
        SPARC_DATABASE, calibration_galaxy, final_params, alpha,
        get_exact_sparc_data, model_velocity_knudsen, quality_filter=quality_filter
    )
    if stats:
        st.sidebar.markdown("---")
        st.sidebar.subheader(f"📊 {title_prefix} ({stats['count']} г.):")
        st.sidebar.metric(label="Среднее значение MAPE", value=f"{stats['mean']:.2f}%")
        st.sidebar.metric(label="Медианное значение MAPE", value=f"{stats['median']:.2f}%")
        if show_best_cal:
            st.sidebar.markdown(f"<p style='font-size:12px; color:#4F8BF9; margin-top:-10px; margin-bottom:10px;'>🏆 Лучшая калибровка: <b>{best_galaxy_found}</b></p>", unsafe_allow_html=True)
        st.sidebar.metric(label="Максимальное MAPE", value=f"{stats['max_value']:.2f}%")
        st.sidebar.markdown(f"<p style='font-size:11px; color:gray; margin-top:-15px; margin-bottom:10px;'>Галактика: <b>{stats['max_name']}</b></p>", unsafe_allow_html=True)
        st.sidebar.metric(label="Минимальное MAPE", value=f"{stats['min_value']:.2f}%")
        st.sidebar.markdown(f"<p style='font-size:11px; color:gray; margin-top:-15px; margin-bottom:10px;'>Галактика: <b>{stats['min_name']}</b></p>", unsafe_allow_html=True)
        st.sidebar.metric(label="Стандартное отклонение", value=f"{stats['std']:.2f}%")

        calculate_and_display_group_stats(
            SPARC_DATABASE, calibration_galaxy, final_params, alpha,
            get_exact_sparc_data, model_velocity_knudsen, calibrate_vacuum_parameters,
            quality_filter=quality_filter
        )

# 1. Полная выборка (все 174 галактики)
render_stats_block("Статистика по выборке", quality_filter=None, show_best_cal=True)

# 2. Эталонная выборка высокого качества (Q = 1)
render_stats_block("Выборка высокого качества Q=1", quality_filter=1, show_best_cal=False)



# --- 5. МАТРИЧНЫЙ ВЫВОД ВСЕХ 6 ГАЛАКТИК ---
#galaxies_list = list(SPARC_DATABASE.keys())
galaxies_list = list(SPARC_DATABASE.keys())[:4]  # Временный лимит на 4 графика


for row_idx in range(int(np.ceil(len(galaxies_list) / 2))):
    col1, col2 = st.columns(2)
    for col_idx, current_col in enumerate([col1, col2]):
        if (row_idx * 2 + col_idx) >= len(galaxies_list):
            break
        g_name = galaxies_list[row_idx * 2 + col_idx]
        
        R, Vobs, Vbar, M_bar, R_d, z_c, M_bh = get_exact_sparc_data(g_name)
        # Рассчитываем скорость
        V_mod = model_velocity_knudsen(final_params, R, Vbar, M_bar, R_d, z_c, M_bh, alpha)
        
        z0_profile = z_c * (1.0 + (R / R_d)**2)
        # Рассчитываем кумулятивную (накопленную) барионную массу в тоннах/массах Солнца по радиусу
        G_const_solar = 4.30091e-6 
        accumulated_mass_g = (Vbar ** 2) * (R * 1000) / G_const_solar / 1e9
        Kn_profile_bh = get_combined_knudsen(lambda_0_fixed, R, M_bar, z0_profile, M_bh)
        

        mape = np.mean(np.abs((Vobs - V_mod) / Vobs)) * 100
        
        with current_col:
            if g_name == calibration_galaxy:
                st.markdown(f"### 🌌 {g_name} 👑 <span style='color:#FFD700; font-size:14px; font-weight:bold;'>[ ОПОРНАЯ ГАЛАКТИКА КАЛИБРОВКИ ]</span>", unsafe_allow_html=True)
            else:
                st.markdown(f"### 🌌 {g_name}")
            
            st.info(f"MAPE: {mape:.2f}%")
            
            fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.65, 0.35], vertical_spacing=0.07, specs=[[{"secondary_y": False}], [{"secondary_y": True}]])

            fig.add_trace(go.Scatter(x=R, y=Vobs, mode='markers', name='SPARC', marker=dict(color='yellow', size=6)), row=1, col=1)
            fig.add_trace(go.Scatter(x=R, y=Vbar, mode='lines', name='Baryons', line=dict(color='blue', dash='dash')), row=1, col=1)
            fig.add_trace(go.Scatter(x=R, y=V_mod, mode='lines', name='Model', line=dict(color='red', width=2.5)), row=1, col=1)
            
            fig.add_trace(go.Scatter(x=R, y=Kn_profile_bh, mode='lines', name='Kn', line=dict(color='orange', width=1.5)), row=2, col=1, secondary_y=False)

            fig.update_layout(height=450, showlegend=False, margin=dict(l=10, r=10, t=10, b=10))
            fig.update_yaxes(title_text="Kn / M_bar", type="log", row=2, col=1, secondary_y=False)

            fig.update_yaxes(title_text="Kn", type="log", row=2, col=1)

            # --- НАША НОВАЯ СТРОКА НАСТРОЙКИ ПРАВОЙ ОСИ ДЛЯ ТОЛЩИНЫ ДИСКА ---
# Локальная привязка и калибровка правой линейной оси y3 для толщины текущей галактики

        # Направляем фиолетовый пунктир строго на правую независимую ось нижнего подграфика
            fig.add_trace(go.Scatter(x=R, y=z0_profile, mode='lines', name='z0 Profile', line=dict(color='violet', width=1.7, dash='dot')), row=2, col=1, secondary_y=True)
        # Явно переименовываем текстовый ярлык правой вертикальной оси на нижнем этаже карточки
            fig.update_yaxes(title_text="z0 (kpc)", secondary_y=True, row=2, col=1)

            # Яркая голубая линия накопленной барионной массы на нижнем этаже
            fig.add_trace(go.Scatter(x=R, y=accumulated_mass_g, mode='lines', name='M_bar(R)', line=dict(color='cyan', width=2)), row=2, col=1)

            fig.update_xaxes(title_text="R (kpc)", row=2, col=1)
            
            st.plotly_chart(fig, use_container_width=True)
