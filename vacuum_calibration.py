import numpy as np
from scipy.optimize import minimize
from vacuum_physics import model_velocity_knudsen

def calibrate_vacuum_parameters(R_cal, Vobs_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal, alpha, lambda_0_fixed):
    
    def loss_function(k_shear_val):
        k_scalar = float(k_shear_val[0])
        if k_scalar < 1e-5: 
            return 1e10
            
        # Строго соблюдаем порядок аргументов физического движка: params, R, Vbar, M_bar...
        V_pred = model_velocity_knudsen(
            [k_scalar, lambda_0_fixed], 
            R_cal, Vbar_cal, M_cal, Rd_cal, zc_cal, Mbh_cal, alpha
        )
        return np.sum((Vobs_cal - V_pred) ** 2)

    # Проводим оптимизацию Nelder-Mead по реальным точкам Vobs_cal
    res = minimize(loss_function, [1.0], method='Nelder-Mead')
    k_shear_calibrated = float(res.x[0])
    
    return k_shear_calibrated
