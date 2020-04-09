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

 