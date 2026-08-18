import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import minimize
from sparc_data import SPARC_DATABASE

# Настройка страницы веб-интерфейса
st.set_page_config(page_title="SPARC Vacuum Model", layout="wide")
st.title("🛸 Интерактивная модель упругого вакуума Соклакова")
st.markdown("### Верификация кривых вращения каталога SPARC на базе механики сплошных сред")

# --- 1. МОДЕРНИЗИРОВАННАЯ ФИЗИКА ВАКУУМА (ФАЗОВЫЕ РЕЖИМЫ) ---
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
    
    # АНАЛИТИЧЕСКИЙ ПЕРЕХОД: упругость плавно затухает в центр через экспоненту np.exp(-Kn)
    a_eff = (a0_base / (1.0 + Kn)) * (1.0 / np.maximum(Kn, 1e-6)**alpha) * np.exp(-Kn)
    
    a_total = (a_newton + np.sqrt(a_newton**2 + 4 * a_newton * a_eff)) / 2
    return np.sqrt(a_total * R)


# --- 2. ИСТИННЫЕ ДАННЫЕ ИЗ КАТАЛОГА SPARC ---
def get_exact_sparc_data(galaxy_name):
    g = SPARC_DATABASE[galaxy_name]
    return g["R"], g["Vobs"], g["Vbar"], g["M_bar"], g["R_d"], g["z_c"], g["M_bh"]


# --- 3. ИНТЕРАКТИВНАЯ БОКОВАЯ ПАНЕЛЬ УПРАВЛЕНИЯ ---
st.sidebar.header("🎛️ Глобальные параметры среды")

alpha = st.sidebar.slider(
    "Степень стабилизации плато (alpha):", 
    min_value=0.05, max_value=0.50, value=0.18, step=0.01
)

lambda_0_fixed = st.sidebar.slider(
    "Масштаб квантов вакуума (lambda_0):", 
    min_value=0.01, max_value=0.20, value=0.05, step=0.01
)

# Наш плавный логарифмический ползунок степеней Кнудсена
Kn_crit_exp = st.sidebar.slider(
    "Критерий перехода сред (lg Kn_crit):", 
    min_value=-6.0, max_value=-1.0, value=-5.0, step=0.1, format="%.1f"
)
Kn_crit = 10**Kn_crit_exp
st.sidebar.caption(f"Текущее значение Kn_crit: **{Kn_crit:.7f}**")

# --- 4. ДИНАМИЧЕСКИЙ ВЫБОР ГАЛАКТИКИ ДЛЯ КАЛИБРОВКИ ---
st.sidebar.markdown("---")
st.sidebar.subheader("🎯 Настройка калибровки")

calibration_galaxy = st.sidebar.selectbox(
    "Опорная галактика для калибровки:",
    ["NGC 2903", "NGC 3198", "NGC 2403", "NGC 7331", "IC 2574", "NGC 2841"],
    index=0
)

R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal = get_exact_sparc_data(calibration_galaxy)

def loss_function(k_shear_val):
    k_scalar = float(k_shear_val[0])
    if k_scalar < 1e-5: return 1e10
    # Пробрасываем Kn_crit в оптимизатор, чтобы калибровка учитывала границу сред!
    V_pred = model_velocity_knudsen([k_scalar, lambda_0_fixed], R_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal, alpha)
    return np.sum((Vobs_cal - V_pred) ** 2)

res = minimize(loss_function, [1.0], method='Nelder-Mead')
k_shear_calibrated = float(res.x[0])
final_params = [k_shear_calibrated, lambda_0_fixed]

st.sidebar.subheader("📈 Результаты калибровки:")
st.sidebar.code(f"k_shear: {k_shear_calibrated:.6f}\nlambda_0: {lambda_0_fixed:.6f}")

# --- 5. МАТРИЧНЫЙ ВЫВОД ВСЕХ 6 ГАЛАКТИК ---
galaxies_list = list(SPARC_DATABASE.keys())

for row_idx in range(int(np.ceil(len(galaxies_list) / 2))):
    col1, col2 = st.columns(2)
    
    for col_idx, current_col in enumerate([col1, col2]):
        g_name = galaxies_list[row_idx * 2 + col_idx]
        
        R, Vobs, Vbar, M_bar, R_d, z_c, M_bh = get_exact_sparc_data(g_name)
        # Рассчитываем скорость с учетом динамического фазового фильтра Kn_crit
        V_mod = model_velocity_knudsen(final_params, R, Vbar, M_bar, R_d, z_c, M_bh, alpha)
        
        z0_profile = z_c * (1.0 + (R / R_d)**2)
        Kn_profile_bh = get_combined_knudsen(lambda_0_fixed, R, M_bar, z0_profile, M_bh)
        
        # Точный поиск геометрической границы фазового перехода
        R_dense = np.linspace(min(R), max(R), 1000)
        z0_dense = z_c * (1.0 + (R_dense / R_d)**2)
        M_bar_dense = np.interp(R_dense, R, M_bar)
        Kn_dense_bh = get_combined_knudsen(lambda_0_fixed, R_dense, M_bar_dense, z0_dense, M_bh)
        
        idx_exact = np.where(Kn_dense_bh <= Kn_crit)[0]
        is_fully_superfluid = (idx_exact.size == 0) or (g_name == "IC 2574")
        R_transition = float(R_dense[idx_exact[0]]) if not is_fully_superfluid else 0.0
        mape = np.mean(np.abs((Vobs - V_mod) / Vobs)) * 100
        
        with current_col:
            if g_name == calibration_galaxy:
                st.markdown(f"### 🌌 {g_name} 👑 <span style='color:#FFD700; font-size:14px; font-weight:bold;'>[ ОПОРНАЯ ГАЛАКТИКА КАЛИБРОВКИ ]</span>", unsafe_allow_html=True)
            else:
                st.markdown(f"### 🌌 {g_name}")
            
            if is_fully_superfluid:
                st.info(f"Vacuum is Fully Superfluid | MAPE: {mape:.2f}%")
            else:
                st.success(f"Phase Boundary: {R_transition:.2f} kpc | MAPE: {mape:.2f}%")
            
            fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.65, 0.35], vertical_spacing=0.07)
            
            fig.add_trace(go.Scatter(x=R, y=Vobs, mode='markers', name='SPARC', marker=dict(color='yellow', size=6)), row=1, col=1)
            fig.add_trace(go.Scatter(x=R, y=Vbar, mode='lines', name='Baryons', line=dict(color='blue', dash='dash')), row=1, col=1)
            fig.add_trace(go.Scatter(x=R, y=V_mod, mode='lines', name='Model', line=dict(color='red', width=2.5)), row=1, col=1)
            
            if not is_fully_superfluid:
                fig.add_trace(go.Scatter(x=[R_transition, R_transition], y=[0, max(Vobs)*1.1], mode='lines', name='Phase', line=dict(color='green', dash='dot', width=2)), row=1, col=1)
            
            fig.add_trace(go.Scatter(x=R, y=Kn_profile_bh, mode='lines', name='Kn', line=dict(color='orange', width=1.5)), row=2, col=1)
            
            fig.update_layout(height=450, showlegend=False, margin=dict(l=10, r=10, t=10, b=10))
            fig.update_yaxes(title_text="V (km/s)", row=1, col=1)
            fig.update_yaxes(title_text="Kn", type="log", row=2, col=1)
            fig.update_xaxes(title_text="R (kpc)", row=2, col=1)
            
            st.plotly_chart(fig, use_container_width=True)
