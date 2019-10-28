import os 
def fort_generator(path, fort12_times):
    if not os.path.exists(path):
        os.makedirs(path)

    file10 = open(path + "/tmpfort.10", "w")
    file12 = open(path + "/fort.12", "w")
    file14 = open(path + "/fort.14", "w")

    fort10 = """ $user_inp

    switch_fo_Cee_on          = $fo_Cee_on         
    switch_f1_dt_on           = $f1_dt_on          
    switch_fo_dt_on           = $fo_dt_on          
    switch_Cee0_iter_on       = $Cee0_iter_on      
    switch_initial_cond_on    = $initial_cond_on   
    switch_heating_cooling_on = $heating_cooling_on
    switch_ei_coll_fix_on     = $ei_coll_fix_on    
    switch_displ_J_on         = $displ_J_on        
    switch_Cee0_Krook_on      = $Cee0_Krook_on     

    switch_disable_force_runtime_size_on = $disable_force_runtime_size_on

    Bz_implicitness_param       = 1.0      
    switch_dvfo_centred_diff_on = $dvfo_centred_diff_on
    switch_dv2f1_centred_diff_on = $dv2f1_centred_diff_on
    switch_dBdt_centred_diff_on = $dBdt_centred_diff_on

    switch_Ohmic_Fara_consis_on = $Ohmic_Fara_consis_on
    switch_Ohmic_all_on         = $Ohmic_all_on        
    switch_hydro_on           = $hydro_on          
    switch_hydro_fo_on        = $hydro_fo_on       
    switch_hydro_f1_on        = $hydro_f1_on       

    wpe_over_nuei = $wpe_over_nuei
    c_over_vte    = $c_over_vte   
    atomic_Z      = $atomic_Z     
    atomic_A      = $atomic_A     

    p_SG      = $p_SG     
    q_SG      = $q_SG     
    p_SG_init = $p_SG_init

    nv = $nv
    nx = $nx
    ny = $ny

    dt   = $dt  
    tmax = $tmax

    x_bc_type     = $x_bc_type    
    xmin          = $xmin         
    xmax          = $xmax         
    grid_x_type   = $grid_x_type  
    grid_x_ratio  = $grid_x_ratio 
    grid_x_alpha  = $grid_x_alpha 
    grid_x_offset = $grid_x_offset

    y_bc_type     = $y_bc_type    
    ymin          = $ymin         
    ymax          = $ymax         
    grid_y_type   = $grid_y_type  
    grid_y_ratio  = $grid_y_ratio 
    grid_y_alpha  = $grid_y_alpha 
    grid_y_offset = $grid_y_offset

    vmax          = $vmax         
    grid_v_type   = $grid_v_type  
    grid_v_ratio  = $grid_v_ratio 
    grid_v_alpha  = $grid_v_alpha 
    grid_v_offset = $grid_v_offset

    do_user_prof_sub = $do_user_prof_sub

    prof_Bz_ave    = $prof_Bz_ave   

    do_user_heating_sub = $do_user_heating_sub

    switch_packed_sparse_on    = .true.
    switch_precomp_mat_cols_on = .false.
    matrix_solver              = 'PETSC'
    matrix_solver_tol          = $matrix_solver_tol         
    nonlin_tol                 = $nonlin_tol                
    nonlin_itmax               = $nonlin_itmax              
    CCee0_00_its_delta         = $CCee0_00_its_delta        

    initial_cond_rel_dt = $initial_cond_rel_dt
    initial_cond_nt     = $initial_cond_nt    

    do_out_data_compress = .false.
    op_time_mon_skip     = 10000
    op_restart_freq      = 100
    $end
    """
    fort12 =  fort12_times
    #------  Create graphics output selection file, "fort.14" (NAMELIST)  -----
    fort14 = """$user_gra_sel
    op_save_on(20,1) = 1
    op_save_on(20,3) = 1
    op_save_on(21,3) = 1
    op_save_on(22,3) = 1
    op_save_on(23,3) = 1
    op_save_on(24,3) = 1
    op_save_on(33,3) = 1
    $end 
    """


    file10.write(fort10)
    file12.write(fort12)
    file14.write(fort14)
