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
      
c      parameter( dir = '/home/abetharan/IMPACT/RUNS/A6_INITIAL_DATA/')
c      parameter( cycle = 'cycle_1/')
c      parameter( ip_bname = 'hydra_impact_testA3' )
      parameter( dir = $PATH)
      parameter( ip_bname = $RUNNAME)
c      fname_tmp = trim(dir)//trim(cycle)//trim(ip_bname)
      fname_tmp = trim(dir)//trim(ip_bname)
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
