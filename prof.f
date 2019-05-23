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
c (05/03/19) AA 
      parameter( dir = '/home/abetharan/IMPACT/RUNS/A6_INITIAL_DATA/')
      parameter( cycle = 'cycle_1/')
      parameter( ip_bname = 'hydra_impact_testA3' )
c      parameter( dir = $PATH)
c      parameter( ip_bname = $RUNNAME)
      fname_tmp = trim(dir)//trim(cycle)//trim(ip_bname)
c      fname_tmp = trim(dir)//trim(ip_bname)
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
