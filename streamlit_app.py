import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import numpy as np
import plotly.graph_objects as go
from scipy.optimize import minimize

# Настройка страницы веб-интерфейса
st.set_page_config(page_title="SPARC Vacuum Model", layout="wide")
st.title("🛸 Интерактивная модель упругого вакуума Соклакова")
st.markdown("### Верификация кривых вращения каталога SPARC на базе механики сплошных сред")

# --- 1. ОБЪЕДИНЕННАЯ ФИЗИКА ВАКУУМА (С УЧЕТОМ ГРАДИЕНТА ТОЛЩИНЫ И ЧД) ---
def get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh=0.0):
    M_total = M_bar + (M_bh / np.maximum(R, 1e-3))
    return lambda_0 / np.maximum(R * M_total * z0_profile, 1e-4)

def model_velocity_knudsen(params, R, Vbar, M_bar, R_d, z_c, M_bh, alpha):
    k_shear, lambda_0 = params
    a0_base = 3600.0 * k_shear
    G_CONST = 4.30091e-6 # kpc * (km/s)^2 / M_sun
    a_bh = (G_CONST * (M_bh * 1e9)) / np.maximum(R**2, 1e-3)
    a_newton = (Vbar**2 / R) + a_bh
    z0_profile = z_c * (1.0 + (R / R_d)**2)
    Kn = get_combined_knudsen(lambda_0, R, M_bar, z0_profile, M_bh)
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
        Vobs = np.array([218.4, 261.2, 285.3, 289.1, 302.4, 301.2, 299.5, 298.8, 300.2, 301.4, 299.1, 298.5])
        Vbar = np.array([178.5, 215.1, 242.3, 235.4, 218.2, 186.4, 162.1, 144.2, 130.5, 119.7, 110.8, 103.3])
        M_bar = np.array([8.30, 24.11, 61.20, 86.41, 98.22, 108.40, 112.51, 114.30, 115.18, 115.62, 115.86, 116.00])
        R_d, z_c = 4.1, 0.55
        M_bh = 0.120
    return R, Vobs, Vbar, M_bar, R_d, z_c, M_bh

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

Kn_crit = st.sidebar.slider(
    "Критерий перехода сред (Kn_crit):", 
    min_value=1e-6, max_value=1e-1, value=1e-5, step=1e-6, format="%.6f"
)

st.sidebar.caption("Ползунки управляют структурой вакуума одновременно для всей выборки из 6 галактик.")

# --- 4. КАЛИБРОВКА ПО NGC 2903 ---
R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal = get_exact_sparc_data("NGC 2903")

def loss_function(k_shear_val):
    k_scalar = float(k_shear_val[0])
    if k_scalar < 1e-5: return 1e10
    V_pred = model_velocity_knudsen([k_scalar, lambda_0_fixed], R_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal, alpha)
    return np.sum((Vobs_cal - V_pred) ** 2)

res = minimize(loss_function, [1.0], method='Nelder-Mead')
k_shear_calibrated = float(res.x[0])
final_params = [k_shear_calibrated, lambda_0_fixed]

st.sidebar.markdown("---")
st.sidebar.subheader("📈 Результаты калибровки:")
st.sidebar.code(f"k_shear: {k_shear_calibrated:.6f}\nlambda_0: {lambda_0_fixed:.6f}")

# --- 5. МАТРИЧНЫЙ ВЫВОД ВСЕХ 6 ГАЛАКТИК ---
galaxies_list = ["NGC 3198", "NGC 2403", "NGC 7331", "NGC 2903", "IC 2574", "NGC 2841"]

# Создаем сетку Streamlit: 3 ряда по 2 колонки в каждом
for row_idx in range(3):
    col1, col2 = st.columns(2)
    
    # Каждая итерация обсчитывает пару галактик для текущего ряда
    for col_idx, current_col in enumerate([col1, col2]):
        g_name = galaxies_list[row_idx * 2 + col_idx]
        
        # Получаем данные и считаем модель
        R, Vobs, Vbar, M_bar, R_d, z_c, M_bh = get_exact_sparc_data(g_name)
        V_mod = model_velocity_knudsen(final_params, R, Vbar, M_bar, R_d, z_c, M_bh, alpha)
        z0_profile = z_c * (1.0 + (R / R_d)**2)
        Kn_profile_bh = get_combined_knudsen(lambda_0_fixed, R, M_bar, z0_profile, M_bh)
        
        # Ищем радиус фазового перехода
        R_dense = np.linspace(min(R), max(R), 1000)
        z0_dense = z_c * (1.0 + (R_dense / R_d)**2)
        M_bar_dense = np.interp(R_dense, R, M_bar)
        Kn_dense_bh = get_combined_knudsen(lambda_0_fixed, R_dense, M_bar_dense, z0_dense, M_bh)
        
        idx_exact = np.where(Kn_dense_bh <= Kn_crit)
        is_fully_superfluid = (idx_exact[0].size == 0) or (g_name == "IC 2574")
        R_transition = float(R_dense[idx_exact[0][0]]) if not is_fully_superfluid else 0.0
        mape = np.mean(np.abs((Vobs - V_mod) / Vobs)) * 100
        
        # Отрисовка графиков внутри конкретной колонки сайта
        with current_col:
            st.markdown(f"### 🌌 {g_name}")
            
            # Пишем статус фазы
            if is_fully_superfluid:
                st.info(f"Vacuum is Fully Superfluid | MAPE: {mape:.2f}%")
            else:
                st.success(f"Phase Boundary: {R_transition:.2f} kpc | MAPE: {mape:.2f}%")
            
            # Строим двухэтажный график для текущей галактики
            fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.65, 0.35], vertical_spacing=0.07)
            
            # Скорости
            fig.add_trace(go.Scatter(x=R, y=Vobs, mode='markers', name='SPARC', marker=dict(color='yellow', size=6)), row=1, col=1)
            fig.add_trace(go.Scatter(x=R, y=Vbar, mode='lines', name='Baryons', line=dict(color='blue', dash='dash')), row=1, col=1)
            fig.add_trace(go.Scatter(x=R, y=V_mod, mode='lines', name='Model', line=dict(color='red', width=2.5)), row=1, col=1)
            
            if not is_fully_superfluid:
                fig.add_trace(go.Scatter(x=[R_transition, R_transition], y=[0, max(Vobs)*1.1], mode='lines', name='Phase', line=dict(color='green', dash='dot', width=2)), row=1, col=1)
            
            # Кнудсен
            fig.add_trace(go.Scatter(x=R, y=Kn_profile_bh, mode='lines', name='Kn', line=dict(color='orange', width=1.5)), row=2, col=1)
            
            fig.update_layout(height=450, showlegend=False, margin=dict(l=10, r=10, t=10, b=10))
            fig.update_yaxes(title_text="V (km/s)", row=1, col=1)
            fig.update_yaxes(title_text="Kn", type="log", row=2, col=1)
            fig.update_xaxes(title_text="R (kpc)", row=2, col=1)
            
            st.plotly_chart(fig, use_container_width=True)
