import math
import numpy as np
from scipy import constants
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")

def free_param_calc(norm_Te, norm_ne, norm_Z, norm_Ar):
    def lambda_ei(n, T, T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0):
        if T * T_norm < 10.00 * Z_norm ** 2:

            result = 23.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) * Z_norm * (T * T_norm) ** (-3.00/2.00))

        else:

            result = 24.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) / (T * T_norm))   

        return result


    #9.10938E-31            # Electron mass in kg
    el_charge = e #1.602189E-19    # Electron charge
    r_b = bohr_radi #5.29E-11              # Bohr radius in metres

    sigma_0 = math.pi * r_b ** 2  # Cross-section normaliation
    gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)

    T_J = norm_Te * el_charge

    gamma_ei_0 = norm_Z ** 2 * gamma_ee_0

    v_t = math.sqrt(2.0 * T_J / me)                  # Velocity normalization
    k = lambda_ei(1.0, 1.0, T_norm = norm_Te, n_norm = norm_ne, Z_norm = norm_Z)
    print(k)
    time_norm = v_t ** 3 / (gamma_ei_0 * (norm_ne/norm_Z) * lambda_ei(1.0, 1.0, T_norm = norm_Te, n_norm = norm_ne, Z_norm = norm_Z)) # Time normalization
    x_0 = v_t * time_norm               # Space normalization

    e_0 = me * v_t / (el_charge * time_norm) # Electric field normalization
    q_0 = me * norm_ne * (v_t ** 3)

    param_dict = {}

    param_dict['ne'] = norm_ne
    param_dict['Te'] = norm_Te
    param_dict['Z'] = norm_Z

    # convert ne to 10**21 cm**-3
    ni = norm_ne / norm_Z
    param_dict['ni'] = ni
    # print('Thermal Velocity is : ', v_t)
    # print('Collision Time is : ', time_norm)
    # print('Collion mfp is : ', x_0)

    param_dict['vte'] = v_t
    param_dict['tau_ei'] = time_norm
    param_dict['nu_ei'] = 1/time_norm
    param_dict['lambda_mfp'] = x_0
    param_dict['qe'] = q_0

    gamma = 5/3
    param_dict['sound_speed'] = np.sqrt((gamma * norm_Te * e) / (norm_Ar *mp)) #SI
    # print('Sound Speed is :  ', param_dict['sound_speed'])
    return(param_dict)











