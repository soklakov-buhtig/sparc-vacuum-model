import streamlit as st
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

# --- 2. ИСТИННЫЕ ДАННЫЕ ИЗ КАТАЛОГА SPARC (Сокращенная база) ---
def get_exact_sparc_data(galaxy_name):
    # Полный набор данных для галактик (NGC 3198, NGC 2403, NGC 7331, NGC 2903, IC 2574, NGC 2841)
    # содержится в исходном коде, здесь представлены заглушки для структуры.
    # Данные включают: R, Vobs, Vbar, M_bar, R_d, z_c, M_bh
    if galaxy_name == "NGC 3198":
        # ... (данные как в оригинале)
        pass 
    # ... (аналогично для других галактик)
    return R, Vobs, Vbar, M_bar, R_d, z_c, M_bh
# --- 3. ИНТЕРАКТИВНАЯ БОКОВАЯ ПАНЕЛЬ УПРАВЛЕНИЯ КОНСТАНТАМИ ---
st.sidebar.header("🎛️ Физические параметры среды")

# Интерактивный выбор исследуемой галактики
galaxy_name = st.sidebar.selectbox(
    "Выберите галактику для анализа:",
    ["NGC 3198", "NGC 2403", "NGC 7331", "NGC 2903", "IC 2574", "NGC 2841"]
)

st.sidebar.markdown("---")

# Ползунки для динамического изменения параметров вакуума в реальном времени
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

st.sidebar.caption("При достижении Kn <= Kn_crit вакуум переходит из ламинарного течения в упругую сетку.")

# --- 4. ДИНАМИЧЕСКАЯ ОДНОПАРАМЕТРИЧЕСКАЯ КАЛИБРОВКА (По NGC 2903) ---
# Подгружаем данные опорной галактики для калибровки масштабного коэффициента жесткости k_shear
R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal = get_exact_sparc_data("NGC 2903")

def loss_function(k_shear_val):
    k_scalar = float(k_shear_val)
    if k_scalar < 1e-5: 
        return 1e10
    # Вычисляем квадратичное отклонение модели от наблюдений с текущим alpha и lambda_0
    V_pred = model_velocity_knudsen([k_scalar, lambda_0_fixed], R_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal, alpha)
    return np.sum((Vobs_cal - V_pred) ** 2)

# Минимизация методом Нелдера-Мида (симплекс-метод)
res = minimize(loss_function, [1.0], method='Nelder-Mead')
k_shear_calibrated = float(res.x)
final_params = [k_shear_calibrated, lambda_0_fixed]

# Отображаем текущие параметры калибровки в интерфейсе
st.sidebar.markdown("---")
st.sidebar.subheader("📈 Результаты калибровки:")
st.sidebar.code(f"k_shear: {final_params[0]:.6f}\nlambda_0: {final_params[1]:.6f}")
# --- 5. ОБСЧЕТ ВЫБРАННОЙ ГАЛАКТИКИ И ПОИСК ГРАНИЦЫ ПЕРЕХОДА ---
# Извлекаем данные и рассчитываем профили Кнудсена (см. полную логику в)
R, Vobs, Vbar, M_baryonic, R_d, z_c, M_bh = get_exact_sparc_data(galaxy_name)
V_mod = model_velocity_knudsen(final_params, R, Vbar, M_baryonic, R_d, z_c, M_bh, alpha)
z0_profile = z_c * (1.0 + (R / R_d)**2)
Kn_profile_bh = get_combined_knudsen(lambda_0_fixed, R, M_baryonic, z0_profile, M_bh)

# Поиск точки перехода R_transition и расчет MAPE
R_dense = np.linspace(min(R), max(R), 2000)
Kn_dense_bh = get_combined_knudsen(lambda_0_fixed, R_dense, np.interp(R_dense, R, M_baryonic), z0_dense, M_bh)
idx_exact = np.where(Kn_dense_bh <= Kn_crit)
R_transition = float(R_dense[idx_exact[0][0]]) if idx_exact[0].size > 0 else float(R_virtual)
is_fully_superfluid = idx_exact[0].size == 0
mape = np.mean(np.abs((Vobs - V_mod) / Vobs)) * 100

# Вывод результатов в Streamlit
st.markdown("---")
col1, col2 = st.columns(2)
with col1:
    if is_fully_superfluid: st.info("Vacuum Is Fully Superfluid")
    else: st.success(f"Superfluid/Elastic на R > {R_transition:.2f} kpc")
with col2: st.metric(label="MAPE", value=f"{mape:.2f} %")

# --- 6. ПОСТРОЕНИЕ ИНТЕРАКТИВНЫХ ГРАФИКОВ (PLOTLY) ---
# Создание двухуровневого графика (кривая вращения + профиль Кнудсена)
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.65, 0.35], specs=[[{"secondary_y": False}], [{"secondary_y": True}]])

# Добавление графиков: SPARC, Модель, Граница, Кнудсен (детальный код в)
fig.add_trace(go.Scatter(x=R, y=Vobs, mode='markers', name='SPARC'), row=1, col=1)
fig.add_trace(go.Scatter(x=R, y=V_mod, mode='lines', name='Model'), row=1, col=1)
if not is_fully_superfluid: fig.add_trace(go.Scatter(x=[R_transition, R_transition], y=[0, max(Vobs)], name='Phase'), row=1, col=1)
fig.add_trace(go.Scatter(x=R, y=Kn_profile_bh, name='Kn with BH'), row=2, col=1, secondary_y=True)

fig.update_layout(height=750, title="Анализ SPARC")
st.plotly_chart(fig, use_container_width=True)
st.markdown("---")
