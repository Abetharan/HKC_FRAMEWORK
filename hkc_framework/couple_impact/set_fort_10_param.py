def set_fort_10(
wpe_over_nuei,
c_over_vte,
atomic_Z,
atomic_A,
nv,
nx,
ny,
dt,
tmax, 
fo_Cee_on  = ".true.",
f1_dt_on = ".true." ,
fo_dt_on = ".true." ,
Cee0_iter_on = ".false." ,
initial_cond_on = ".false." ,
heating_cooling_on = ".false." ,
ei_coll_fix_on = ".true." ,
displ_J_on = ".true." ,
Cee0_Krook_on = ".false." ,
disable_force_runtime_size_on = ".false." ,
dvfo_centred_diff_on = ".true." ,
dv2f1_centred_diff_on = ".true." ,
dBdt_centred_diff_on = ".true." ,
Ohmic_Fara_consis_on = ".true." ,
Ohmic_all_on = ".false." ,
hydro_on = ".false." ,
hydro_fo_on = ".false." ,
hydro_f1_on = ".false." ,
p_SG = "2.0d0" ,
q_SG = "2.0d0",
p_SG_init ="2.0d0" ,
x_bc_type = "0",
xmin = "0.0d0" ,
xmax = "159424.7662d0",
grid_x_type = "1",
grid_x_ratio = "1.0d0",
grid_x_alpha = "0.0d0",
grid_x_offset = "0.0d0",
y_bc_type = "0",
ymin = "-3.0d3",
ymax = "+3.0d3",
grid_y_type = "0",
grid_y_ratio = "1.0d0",
grid_y_alpha ="0.0d0",
grid_y_offset = "0.0d0",
vmax = "13.0d0",
grid_v_type = "2",
grid_v_ratio = "1.0d0",
grid_v_alpha = "0.0d0",
grid_v_offset = "0.0d0" ,
do_user_prof_sub = ".false.",
prof_Bz_ave = "0.0d0",
do_user_heating_sub = ".true.",
matrix_solver = "'PETSC'",
matrix_solver_tol = "1.0d-15",
nonlin_tol = "1.0d-12",
nonlin_itmax = "25",
CCee0_00_its_delta = "10",
initial_cond_rel_dt = "1.0d-2",
initial_cond_nt = 2,
op_restart_freq = 100):

    kappa = {
    'user_inp' : "$user_inp",
    'fo_Cee_on':fo_Cee_on,
    'f1_dt_on':f1_dt_on,
    'fo_dt_on':fo_dt_on,
    'Cee0_iter_on':Cee0_iter_on,
    'initial_cond_on':initial_cond_on,
    'heating_cooling_on':heating_cooling_on,
    'ei_coll_fix_on':ei_coll_fix_on,
    'displ_J_on':displ_J_on,
    'Cee0_Krook_on':Cee0_Krook_on,
    'disable_force_runtime_size_on':disable_force_runtime_size_on,
    'dvfo_centred_diff_on':dvfo_centred_diff_on,
    'dv2f1_centred_diff_on':dv2f1_centred_diff_on,
    'dBdt_centred_diff_on':dBdt_centred_diff_on,
    'Ohmic_Fara_consis_on':Ohmic_Fara_consis_on,
    'Ohmic_all_on':Ohmic_all_on,
    'hydro_on':hydro_on,
    'hydro_fo_on':hydro_fo_on,
    'hydro_f1_on':hydro_f1_on,
    'wpe_over_nuei':wpe_over_nuei,
    'c_over_vte':c_over_vte,
    'atomic_Z':atomic_Z,
    'atomic_A':atomic_A,
    'p_SG':p_SG,
    'q_SG':q_SG,
    'p_SG_init':p_SG_init,
    'nv':nv,
    'nx':nx,
    'ny':ny,
    'dt':dt,
    'tmax':tmax,
    'x_bc_type':x_bc_type,
    'xmin':xmin,
    'xmax':xmax,
    'grid_x_type':grid_x_type,
    'grid_x_ratio':grid_x_ratio,
    'grid_x_alpha':grid_x_alpha,
    'grid_x_offset':grid_x_offset,
    'y_bc_type':y_bc_type,
    'ymin':ymin,
    'ymax':ymax,
    'grid_y_type':grid_y_type,
    'grid_y_ratio':grid_y_ratio,
    'grid_y_alpha':grid_y_alpha,
    'grid_y_offset':grid_y_offset,
    'vmax':vmax,
    'grid_v_type':grid_v_type,
    'grid_v_ratio':grid_v_ratio,
    'grid_v_alpha':grid_v_alpha,
    'grid_v_offset':grid_v_offset,
    'do_user_prof_sub':do_user_prof_sub,
    'prof_Bz_ave':prof_Bz_ave,
    'do_user_heating_sub':do_user_heating_sub,
    'matrix_solver_tol':matrix_solver_tol,
    'nonlin_tol':nonlin_tol,
    'nonlin_itmax':nonlin_itmax,
    'CCee0_00_its_delta': CCee0_00_its_delta,
    'initial_cond_rel_dt': initial_cond_rel_dt,
    'initial_cond_nt': initial_cond_nt,
    'op_restart_freq':op_restart_freq,
    'end' : "$end"
    }
    return(kappa)
