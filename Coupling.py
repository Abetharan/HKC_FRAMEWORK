
import numpy as np 
file = open("fort.10", "w")
fort10txt = """ $user_inp

        switch_fo_Cee_on          = .true.
        switch_f1_dt_on           = .true.
        switch_fo_dt_on           = .true.
        switch_Cee0_iter_on       = .FALSE. 
        switch_initial_cond_on    = .false.
        switch_heating_cooling_on = .TRUE.
        switch_ei_coll_fix_on     = .true.
        switch_displ_J_on         = .false.
        switch_Cee0_Krook_on      = .false.

        switch_disable_force_runtime_size_on = .false.

        Bz_implicitness_param       =  1.0d0
        switch_dvfo_centred_diff_on = .true.  
        switch_dv2f1_centred_diff_on = .true.
        switch_dBdt_centred_diff_on = .true.

        switch_Ohmic_Fara_consis_on = .true.
        switch_Ohmic_all_on         = .false.
        switch_hydro_on           = .false.
        switch_hydro_fo_on        = .false.
        switch_hydro_f1_on        = .false.

        wpe_over_nuei = 11.14576d0
        c_over_vte    = 30.455466d0
        atomic_Z      = 37.664907d0
        atomic_A      = 70.0d0

        p_SG      = 2.0d0
        q_SG      = 2.0d0
        p_SG_init = 2.0d0

        nv = 300
        nx = 6
        ny = 1

        dt   =     0.5d0
        tmax =   750.0d0

        x_bc_type     = 1
        xmin          = 0.0d0
        xmax          = 159424.7662d0
        grid_x_type   = 3
        grid_x_ratio  =  1.0d0
        grid_x_alpha  =  0.0d0
        grid_x_offset =  0.0d0

        y_bc_type     = 0
        ymin          = -3.0d3
        ymax          = +3.0d3
        grid_y_type   = 0
        grid_y_ratio  = 1.0d0
        grid_y_alpha  = 0.0d0
        grid_y_offset = 0.0d0

        vmax          = 30.0d0
        grid_v_type   =  2
        grid_v_ratio  = 30.0d0
        grid_v_alpha  =  0.0d0
        grid_v_offset =  0.0d0

        do_user_prof_sub = .TRUE.

        prof_Bz_ave    = 0.0d0

        do_user_heating_sub = .TRUE.

        switch_packed_sparse_on    = .true.
        switch_precomp_mat_cols_on = .false.
        matrix_solver              = 'PETSC'
        matrix_solver_tol          = 1.0d-15
        nonlin_tol                 = 1.0d-12
        nonlin_itmax               = 25
        CCee0_00_its_delta         = 10

        initial_cond_rel_dt = 1.0d-2
        initial_cond_nt     = 2

        do_out_data_compress = .false.
        op_time_mon_skip     = 10000
        op_restart_freq      = 100
 $end
"""

prof = """
${RUN}_prof.f"
c234567890---------20--------3---------4---------5---------6---------7-2-----8
      subroutine init_prof_user

      implicit none

c ...  Common blocks  ...
      include 'fp2df1_dims.com.f'
      include 'fp2df1_grids.com.f'
      include 'fp2df1_results.com.f'
      include 'fp2df1_user_params.com.f'
      include 'fp2df1_io.com.f'

c ... COMMON BLOCK for sharing heating profile (LOCAL sized arrays) ...
c (09/11/16) RK  --  Adapted from Polar1_*.sh series
      double precision vol_laser_heating_rate(ccg_ydl:ccg_ydu,
     +                                        ccg_xdl:ccg_xdu)
      double precision vol_rad_cooling_rate(ccg_ydl:ccg_ydu,
     +                                      ccg_xdl:ccg_xdu)
      common /heat_prof/vol_laser_heating_rate,vol_rad_cooling_rate

c ...  Local variables  ...
      integer ix,iy, ix_glb
      integer bn_ll
      double precision ne, te, ni, ZZ, Bz, vte
c GLOBAL sized arrays for storage of file data
      double precision
     +     Te_arr(ccg_ydl:ccg_ydu,ccg_xdl_glb:ccg_xdu_glb),
     +     ne_arr(ccg_ydl:ccg_ydu,ccg_xdl_glb:ccg_xdu_glb),
     +     ni_arr(ccg_ydl:ccg_ydu,ccg_xdl_glb:ccg_xdu_glb),
     +      Z_arr(ccg_ydl:ccg_ydu,ccg_xdl_glb:ccg_xdu_glb),
     +     vhr_laser_arr(ccg_ydl:ccg_ydu,ccg_xdl_glb:ccg_xdu_glb),
     +     vhr_rad_arr(ccg_ydl:ccg_ydu,ccg_xdl_glb:ccg_xdu_glb)
     

      character*100 fname

c ...  Function  ...
      integer ix_to_glb
      logical cc_redundant_cell

c (31/07/17) RK
      character*60 ip_bname
      character*61 dir 
      character*61 cycle
      character*120 fname_tmp
      
      parameter( dir = '/home/abetharan/IMPACT/RUNS/A6_INITIAL_DATA/')
      parameter( cycle = 'cycle_1/')
      parameter( ip_bname = 'hydra_impact_testA3' )

      fname_tmp = trim(dir)//trim(cycle)//trim(ip_bname)
c----------------------------------------------------------------------------
       write(os,*)
       write(os,'(1x,80(''*''))')
       write(os,'(1x,20x,a)')' LOADING for HYDRA 1D Gd+He Surrogate '//
     +                      'Hohlraum (initial conditions)'
       write(os,'(1x,80(''*''))')

c ----- Load Data -----

      fname = trim(fname_tmp)//'_tmat.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, Te_arr,
     +       ccg_xdl_glb,ccg_xdu_glb,ccg_ydl,ccg_ydu,
     +       ccg_xdl_glb,ccg_xdu_glb, ccg_yl, ccg_yu)

c ne
      fname = trim(fname_tmp)//'_eden.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, ne_arr,
     +       ccg_xdl_glb,ccg_xdu_glb,ccg_ydl,ccg_ydu,
     +       ccg_xdl_glb,ccg_xdu_glb, ccg_yl, ccg_yu)

c ni
      fname = trim(fname_tmp)//'_ionden.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, ni_arr,
     +       ccg_xdl_glb,ccg_xdu_glb,ccg_ydl,ccg_ydu,
     +       ccg_xdl_glb,ccg_xdu_glb, ccg_yl, ccg_yu)

c Z
      fname = trim(fname_tmp)//'_zstar.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, Z_arr,
     +       ccg_xdl_glb,ccg_xdu_glb,ccg_ydl,ccg_ydu,
     +       ccg_xdl_glb,ccg_xdu_glb, ccg_yl, ccg_yu)

c Laser volumetric heating rate (in IMPACT norms)
      fname = trim(fname_tmp)//'_laserdep.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, vhr_laser_arr,
     +       ccg_xdl_glb,ccg_xdu_glb,ccg_ydl,ccg_ydu,
     +       ccg_xdl_glb,ccg_xdu_glb, ccg_yl, ccg_yu)

c Electron radiative cooling rate (in IMPACT norms)
      fname = trim(fname_tmp)//'_rad_to_electron.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, vhr_rad_arr,
     +       ccg_xdl_glb,ccg_xdu_glb,ccg_ydl,ccg_ydu,
     +       ccg_xdl_glb,ccg_xdu_glb, ccg_yl, ccg_yu)

c ...  Clear  ...
      vol_laser_heating_rate = 0.0d0
      vol_rad_cooling_rate   = 0.0d0

c -----  Set fo to Maxwellian  -----
      do ix=ccg_xl,ccg_xu
	     ix_glb = ix_to_glb(ix)
         do iy=ccg_yl,ccg_yu
            ne = ne_arr(iy,ix_glb)
            te = Te_arr(iy,ix_glb)
            ZZ =  Z_arr(iy,ix_glb)
            ni = ni_arr(iy,ix_glb)
            Bz = prof_Bz_ave
c Copy to LOCAL heating arrays (in common block)
            vol_laser_heating_rate(iy,ix) = vhr_laser_arr(iy,ix_glb)
            vol_rad_cooling_rate(iy,ix)   = vhr_rad_arr(iy,ix_glb)

            vte=sqrt(2.0d0*te)             
c ...  Initialise  ...
            if (.not.cc_redundant_cell(ix,iy)) then
               call init_fo_vc(iy,ix,vte,ne)
            endif
            Res_Bz(iy,ix,tl_new) = Bz
            fixed_prof_ni(iy,ix) = ni
            fixed_prof_Z(iy,ix)  = ZZ
          enddo
       enddo
       return
       end
$end
"""

heating_routine = """ 
c234567890---------20--------3---------4---------5---------6---------7-2-----8
c
c  Adapted from Polar1_*.sh series
c
c  Simple approach:
c  - Total heating rate is sum if laser + cooling.
c  - This is sent to the MW operator.
c
c----------------------------------------------------------------------------+

      subroutine Cee0_00_heating_user

      implicit none

c ...  Common blocks  ...
      include 'fp2df1_dims.com.f'
      include 'fp2df1_grids.com.f'
      include 'fp2df1_Cee0_00.com.f'
      include 'fp2df1_prog_cntr.com.f'
      include 'fp2df1_results.com.f'
      include 'fp2df1_diagno.com.f'

c ... COMMON BLOCK for sharing heating profile  ...
      double precision vol_laser_heating_rate(ccg_ydl:ccg_ydu,
     +                                        ccg_xdl:ccg_xdu)
      double precision vol_rad_cooling_rate(ccg_ydl:ccg_ydu,
     +                                      ccg_xdl:ccg_xdu)
      common /heat_prof/vol_laser_heating_rate,vol_rad_cooling_rate


c ...  Local variables  ...
      logical Maxw_heating
      integer ixc, iyc,ivb
      double precision vo2, en_in, en_in_novol
      double precision tdep, time_heat_start
      double precision fo_at_v0(fo_ydl:fo_ydu,fo_xdl:fo_xdu)
      double precision heat_coef(cb_vdl:cb_vdu)
c
      double precision pi, edr, dTdt, Te_max
c (27/04/09)
      double precision edr_max, dTdt_max, edr_min, dTdt_min
c GLOBAL versions of accounting parameters
      double precision en_in_glb, en_in_novol_glb
      double precision edr_max_glb, dTdt_max_glb
      double precision edr_min_glb, dTdt_min_glb
c ...  Static data  ...
      double precision en_in_tot_novol_glb
      double precision en_in_tot_glb

c ...  Parameters  ...
      parameter( time_heat_start = 0.1d0 ,  Te_max = 20.0d0 )
      parameter( Maxw_heating = .false. )

c ...  Data  ...
      data en_in_tot_glb,en_in_tot_novol_glb/0.0d0,0.0d0/

c----------------------------------------------------------------------+
c        Calculate IB heating   (CORRECTED  29/08/06)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c                   d(Ue)    3    d(Te)
c   en_dep_rate  =   --   =  - ne  --  
c                    dt      2     dt
c
c Relation of d(Te)/dt to IB parameters (IMPACT norms) is:
c
c     d(Te)      4         (v_osc)^2         [fo = f_Maxw]
c     --    =  ----------  -----------
c     dt       9 sqrt{pi}  (2Te)^(3/2)
c
c More generally:
c
c     d(Ue)    2 pi 
c     --    =  ----  (v_osc)^2  fo(v=0)      [fo = arb.]
c     dt        3  
c
c Hence, will calculate necessary  v_osc  to yield supplied EDR
c  given fo(v=0) and ne.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c        Calculate Maxwellian heating 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c                   d(Ue)    3    d(Te)
c   en_dep_rate  =   --   =  - ne  --  
c                    dt      2     dt
c
c Relation of d(Te)/dt to MW parameters (IMPACT norms) is:
c
c     dT/dt  = 2 * D_o / Z(x,y)    (in IMPACT norms)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+


      edr_max  = -123456789.0d+100
      dTdt_max = edr_max
      edr_min  = abs(edr_max)
      dTdt_min = edr_min

      pi = acos(-1.0d0)

      tdep=0.0d0
      if (time.ge.time_heat_start) tdep=1.0d0
        

      en_in       = 0.0d0
      en_in_novol = 0.0d0

      call unpk_fo_xy(tl_old,cc_vl,fo_at_v0)

      do 10 ixc=fo_xl,fo_xu
         do 20 iyc=fo_yl,fo_yu

c ...  Sum laser heating & elec. rad. cooling rates  ...
c Note:
c  vol_laser_heating_rate() is  volumetric heating rate due to laser IB deposition.
c  vol_rad_cooling_rate()   is  net volumetric heating rate of the electrons due 
c         to radiation absorption.  So a negative value is net radiation loss.
c
            edr = vol_laser_heating_rate(iyc,ixc) 
     +                 + vol_rad_cooling_rate(iyc,ixc) 

c ...  Calc.  dT/dt  equiv. to EDR  ...
            dTdt = tdep * edr / ((3.0d0/2)*ne(iyc,ixc))


            if (.not.Maxw_heating) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c  Set IB heating operator
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c ...  Calc.  (v_osc)^2  ...
           vo2=tdep*edr*(3.0d0/2.0d0)/(pi*fo_at_v0(iyc,ixc))

c ...  Compute then apply local IB operator H_IB(v)  ...
            call Cee0_00_heating_IB_red(vo2,heat_coef)
            do ivb=cb_vl,cb_vu
               CCee0_00_diff(ivb,iyc,ixc)=CCee0_00_diff(ivb,iyc,ixc)+
     +              heat_coef(ivb)
            enddo

        else
        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c  Set MW heating operator
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
c ...  Compute then apply local MW operator H_Mw(v)  ...
	    call Cee0_00_heating_maxwell(dTdt,Te_max,ixc,iyc,heat_coef)
	    
            do  ivb=cb_vl,cb_vu
               CCee0_00_diff(ivb,iyc,ixc)=CCee0_00_diff(ivb,iyc,ixc)+
     +              heat_coef(ivb)
            enddo

         endif

c ...  Account energy  ...            
            en_in_novol = en_in_novol + edr*dt*tdep
            en_in       = en_in + (edr*dt*tdep)*dxc(ixc)*dyc(iyc)

c (27/04/09)
            edr_max =  max(edr, edr_max)
            dTdt_max = max(dTdt, dTdt_max)
            edr_min =  min(edr, edr_min)
            dTdt_min = min(dTdt, dTdt_min)

 20       enddo
 10   enddo

c ...  Collection reduction over local statistics  ...
      call para_data_sum( en_in, en_in_glb, 1 )
      call para_data_sum( en_in_novol, en_in_novol_glb,1)
      call para_data_calc(  edr_max,  edr_max_glb, 1, 'MPI_MAX' )
      call para_data_calc( dTdt_max, dTdt_max_glb, 1, 'MPI_MAX' )
      call para_data_calc(  edr_min,  edr_min_glb, 1, 'MPI_MIN' )
      call para_data_calc( dTdt_min, dTdt_min_glb, 1, 'MPI_MIN' )


      en_in_tot_novol_glb = en_in_tot_novol_glb + en_in_novol_glb
      en_in_tot_glb       = en_in_tot_glb + en_in_glb

    
      write(os,*)
      write(os,'(1x,40(''- ''))')
      write(os,*)'  Cee0_heating_user'
      write(os,*)'  -----------------'
      write(os,*)'     Maxellian heating? (else IB) = ',Maxw_heating
      write(os,*)'                             time = ',time
      write(os,*)'        time for heating to start = ',time_heat_start
      write(os,*)'                             tdep = ',tdep 
      write(os,*)'             dT/dt (max, min GLB) = ',dTdt_max_glb,
     +                                                  dTdt_min_glb
      write(os,*)'  Energy dep. rate (max, min GLB) = ',edr_max_glb,
     +                                                  edr_min_glb
      write(os,*)'  Energy (GLB) in this dt (novol) = ',en_in_novol_glb
      write(os,*)'  Energy (GLB) in so far (novol)  = ',
     +                                             en_in_tot_novol_glb
      write(os,*)'  Energy (GLB) in this dt         = ',en_in_glb
      write(os,*)'  Energy (GLB) in so far          = ',en_in_tot_glb
      write(os,'(1x,40(''- ''))')

      return
      end

c234567890---------20--------3---------4---------5---------6---------7-2-----8

      subroutine Cee0_00_cooling_user
      implicit none
      return
      end 
      
$end
"""

user_custom_routine = """
c234567890---------20--------3---------4---------5---------6---------7-2-----8
c
c RK  (08/11/16)   --  Load specific spatial grid 
c
c Notes:
c  - Read xb non-uniform grid from files.
c    (IMPACT deals with creating xc & xcg grids.)
ccc  - Files are standard *.xy format file that have xb & xcg as data!
c
c234567890---------20--------3---------4---------5---------6---------7-2-----8

      subroutine grid_custom(n,xin,xout,char_ind)
      implicit none

c ...  Common blocks  ...
      include 'fp2df1_dims.com.f'
      include 'fp2df1_io.com.f'

c ...  Passed arguments  ...
      integer n
      double precision xin(1:n), xout(1:n)
      character char_ind*2

c ...  Local variables  ...
      character*100 fname
      logical time_mon
      integer ndims, nd(4)
      integer ix,iy
      integer nxb,  dest_proc
      parameter( nxb=nxbm_glb )
      double precision dxc, grid_loaded_xb(1:nycm, 1:nxb )
c
	  integer ng_max
	  parameter( ng_max = max(nxm_glb,nym)+5 )
      double precision grid_arr(ng_max,2), grid_arr_root(ng_max,2)
      double precision sample_time

c (31/07/17) RK
      character*60 ip_bname
      character*61 dir 
      character*61 cycle
      character*120 fname_tmp
c (05/03/19) AA 
      parameter( dir = '/home/abetharan/IMPACT/RUNS/A6_INITIAL_DATA/')
      parameter( cycle = 'cycle_1/')
      parameter( ip_bname = 'hydra_impact_testA3' )

      fname_tmp = trim(dir)//trim(cycle)//trim(ip_bname)
c----------------------------------------------------------------------+
      if (char_ind.ne.'xx') then
         write(os,*)'ERROR: this grid_custom() routine only works'//
     :        ' for the x-grid'
         stop
      endif

      grid_loaded_xb = 0.0d0

c DEBUG  --->
      write(os,*)'Starting grid_custom()'
c <---  DEBUG

c ----- Load Data (every processor!)  -----
      write(os, *)'Tmp name', fname_tmp
      fname = trim(fname_tmp)//'_xf.xy'
      write(os,*)'Loading: ',fname
      call init_prof_user_load_xy(fname, grid_loaded_xb,
     +       cb_xdl_glb,cb_xdu_glb,cc_ydl,cc_ydu,
     +       cb_xdl_glb,cb_xdu_glb, cc_yl, cc_yu)
      nd(1) = cc_yu  - cc_yl  + 1
      nd(2) = cb_xdu_glb - cb_xdl_glb + 1

c ...  Copy into xout  ...
      iy = 1
      write(os,*)'nxb_requested, nxb_file, nxb = ', n, nd(2), nxb
c (05/11/16) RK  --  Deal with IMPACTs new protocol for setting custom grids
	  if (n.eq.nxb) then
c (05/11/16) RK  --  For v4.11.7 and earlier
	    write(os,*)'Assuming v4.11.7-: Ghost cell outer boundaries'//
     +           ' are redundant'
        do ix=1,nxb
    	  xout(ix) = grid_loaded_xb(iy,ix)
        enddo
	  elseif ((n.eq.nxb+2).and.(nd(2).eq.nxb)) then
c (05/11/16) RK  --  For v4.11.8 onwards; must set outer boundary of ghost-cells
c                    (Not in loaded data, so generate.)
	    write(os,*)'Assuming v4.11.8+: Ghost cell outer boundaries'//
     +           ' being generated by extension.'
      	do ix=1,nxb
          xout(ix+1) = grid_loaded_xb(iy,ix)
        enddo
        dxc = xout(3)-xout(2)
        xout(1) = xout(2) - dxc
        dxc = xout(nxb+1)-xout(nxb)
        xout(nxb+2) = xout(nxb+1) + dxc
	  elseif ((n.eq.nxb+2).and.(nd(2).eq.nxb+2)) then
	    write(os,*)'Assuming v4.11.8+: Ghost cell outer boundaries'//
     +           ' specified in file.'
      	do ix=1,nxb+2
          xout(ix) = grid_loaded_xb(iy,ix)
        enddo

	  else
         write(os,*)'ERROR: length of loaded data & return array'//
     +         ' incompatible!'
         stop
	  endif

      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +
 999  write(os,'(1x,a,a)')'The following file not found:  ',fname

      close(UNIT=io_gen_unit)
   	  stop
      end

c234567890---------20--------3---------4---------5---------6---------7-2-----8
      subroutine user_custom_tloop
      return
      end

      subroutine user_custom_nl_it
      return
      end 
$end
"""

fort12 =""" 
    7
   0.1d0
   2.0d0
 150.0d0
 300.0d0
 450.0d0
 600.0d0
 750.0d0
   0.0d0
$end
 """
#------  Create graphics output selection file, "fort.14" (NAMELIST)  -----
fort14 = """ 
 \$user_gra_sel

    op_save_on(20,1) = 1
    op_save_on(20,3) = 1
    op_save_on(21,3) = 1
    op_save_on(22,3) = 1
    op_save_on(23,3) = 1
    op_save_on(24,3) = 1

    op_save_on(33,3) = 1

 \$end 
"""


file.write(fort10txt)








