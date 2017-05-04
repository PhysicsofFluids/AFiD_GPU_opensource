!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitVariables
      use param
      use local_arrays
      use stat_arrays
      use decomp_2d
      use AuxiliaryRoutines
      implicit none
#ifdef USE_CUDA
      type(cudaDeviceProp) :: prop
      integer :: dev
#endif      
!-------------------------------------------------
! Arrays for grid making
!-------------------------------------------------

      call AllocateReal1DArray(zc,1,nz)
      call AllocateReal1DArray(zm,1,nz)
      call AllocateReal1DArray(ak1,1,nz)
      call AllocateReal1DArray(ao,1,nz)

      call AllocateReal1DArray(yc,1,ny)
      call AllocateReal1DArray(ym,1,ny)
      call AllocateReal1DArray(ak2,1,ny)
      call AllocateReal1DArray(ap,1,ny)

      call AllocateReal1DArray(xc,1,nx)
      call AllocateReal1DArray(xm,1,nx)
      call AllocateReal1DArray(g3rc,1,nx)
      call AllocateReal1DArray(g3rm,1,nx)

      call AllocateReal1DArray(udx3c,1,nx)
      call AllocateReal1DArray(udx3m,1,nx)

      call AllocateReal1DArray(ap3ck,1,nx)
      call AllocateReal1DArray(ac3ck,1,nx)
      call AllocateReal1DArray(am3ck,1,nx)

      call AllocateReal1DArray(ap3sk,1,nx)
      call AllocateReal1DArray(ac3sk,1,nx)
      call AllocateReal1DArray(am3sk,1,nx)

      call AllocateReal1DArray(ap3ssk,1,nx)
      call AllocateReal1DArray(ac3ssk,1,nx)
      call AllocateReal1DArray(am3ssk,1,nx)

      call AllocateReal1DArray(amphk,1,nx)
      call AllocateReal1DArray(acphk,1,nx)
      call AllocateReal1DArray(apphk,1,nx)
 
      call AllocateInt1dArray(kmc,1,nx)
      call AllocateInt1dArray(kpc,1,nx)
      call AllocateInt1dArray(kmv,1,nx)
      call AllocateInt1dArray(kpv,1,nx)

!-------------------------------------------------
! Arrays for implicit update 
!-------------------------------------------------
!JR   Arrays preallocated here to avoid device synchronization caused by
!     static array allocation

      call AllocateReal1DArray(amkl,1,nx)
      call AllocateReal1DArray(apkl,1,nx)
      call AllocateReal1DArray(ackl,1,nx)
      call AllocateReal1DArray(fkl,1,nx)

!-------------------------------------------------
! Arrays for temperature boundary conditions    
!-------------------------------------------------

      call AllocateReal2DArray(tempbp,1,ny,1,nz)
      call AllocateReal2DArray(temptp,1,ny,1,nz)

!-------------------------------------------------
! Arrays for statistics    
!-------------------------------------------------

      if (statcal) then
       stat_columns = 1
#if defined(USE_CUDA)
       stat_columns = xend(2) - xstart(2) + 1
       istat = cudaGetDevice(dev)
       istat = cudaGetDeviceProperties(prop,dev)
       !if(ismaster) print *,"SMs = ",prop%multiProcessorCount
       do while(stat_columns*(((nxm+63)/64)*64) > 1024*prop%multiProcessorCount)
         stat_columns = stat_columns-1
       enddo
#endif

       call AllocateReal2DArray(vx_m1,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vy_m1,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vz_m1,1,nxm,1,stat_columns)
       call AllocateReal2DArray(tp_m1,1,nxm,1,stat_columns)

       call AllocateReal2DArray(vx_m2,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vy_m2,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vz_m2,1,nxm,1,stat_columns)
       call AllocateReal2DArray(tp_m2,1,nxm,1,stat_columns)

       call AllocateReal2DArray(vx_m3,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vy_m3,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vz_m3,1,nxm,1,stat_columns)
       call AllocateReal2DArray(tp_m3,1,nxm,1,stat_columns)

       call AllocateReal2DArray(vx_m4,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vy_m4,1,nxm,1,stat_columns)
       call AllocateReal2DArray(vz_m4,1,nxm,1,stat_columns)
       call AllocateReal2DArray(tp_m4,1,nxm,1,stat_columns)

       call AllocateReal2DArray(tpvx_m1,1,nxm,1,stat_columns)

       if (disscal) then
        call AllocateReal2DArray(disste,1,nxm,1,stat_columns)
        call AllocateReal2DArray(dissth,1,nxm,1,stat_columns)
       end if

      end if

      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
!JR   Using basic allocation of these variables to preserve pinned attribute
      allocate(vy(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
      allocate(vz(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
      allocate(vx(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
      allocate(pr(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
      allocate(temp(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
      allocate(dphhalo(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))

      !-----------------------------------------------
      ! Arrays without ghost cells
      !-----------------------------------------------
      call AllocateReal3DArray(rhs,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dq,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(qcap,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(rux,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruy,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruz,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(hro,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(rutemp,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      allocate(dph(1:nxm,xstart(2):xend(2),xstart(3):xend(3)))


!JR   Initialize data allocated using basic allocate statement
      vx = 0.0_fp_kind
      vy = 0.0_fp_kind
      vz = 0.0_fp_kind
      pr = 0.0_fp_kind
      temp = 0.0_fp_kind
      dphhalo = 0.0_fp_kind
      dph = 0.0_fp_kind


#ifdef USE_CUDA

!     Allocate the variables on GPU.
!     To work around a compiler bug, instead of calling the equivalent of AllocateRealXDArray
!     we use  a F2003 extension of allocate: 
!           allocate(array_gpu, source=array_cpu)
!     It will allocate array_gpu  with the same bounds of arrays_cpu and will also copy the values.
!     Since these arrays are defined with the device attribute, they will be allocated on the GPU

      if (statcal) then
       allocate(vx_m1_d,source=vx_m1)
       allocate(vy_m1_d,source=vy_m1)
       allocate(vz_m1_d,source=vz_m1)
       allocate(tp_m1_d,source=tp_m1)

       allocate(vx_m2_d,source=vx_m2)
       allocate(vy_m2_d,source=vy_m2)
       allocate(vz_m2_d,source=vz_m2)
       allocate(tp_m2_d,source=tp_m2)

       allocate(vx_m3_d,source=vx_m3)
       allocate(vy_m3_d,source=vy_m3)
       allocate(vz_m3_d,source=vz_m3)
       allocate(tp_m3_d,source=tp_m3)

       allocate(vx_m4_d,source=vx_m4)
       allocate(vy_m4_d,source=vy_m4)
       allocate(vz_m4_d,source=vz_m4)
       allocate(tp_m4_d,source=tp_m4)

       allocate(tpvx_m1_d,source=tpvx_m1)

       if(disscal) then
         allocate(disste_d,source=disste)
         allocate(dissth_d,source=dissth)
       endif

      endif

#ifdef USE_HYBRID
!JR If using hybrid, only allocate partial pencil data on GPU
      allocate(vx_d(1:nx,xstart_gpu(2)-lvlhalo:xend_gpu(2)+lvlhalo,xstart_gpu(3)-lvlhalo:xend_gpu(3)+lvlhalo))
      allocate(vy_d(1:nx,xstart_gpu(2)-lvlhalo:xend_gpu(2)+lvlhalo,xstart_gpu(3)-lvlhalo:xend_gpu(3)+lvlhalo))
      allocate(vz_d(1:nx,xstart_gpu(2)-lvlhalo:xend_gpu(2)+lvlhalo,xstart_gpu(3)-lvlhalo:xend_gpu(3)+lvlhalo))
      allocate(pr_d(1:nx,xstart_gpu(2)-lvlhalo:xend_gpu(2)+lvlhalo,xstart_gpu(3)-lvlhalo:xend_gpu(3)+lvlhalo))
      allocate(temp_d(1:nx,xstart_gpu(2)-lvlhalo:xend_gpu(2)+lvlhalo,xstart_gpu(3)-lvlhalo:xend_gpu(3)+lvlhalo))

      allocate(rux_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))
      allocate(ruy_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))
      allocate(ruz_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))
      allocate(hro_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))
      allocate(rutemp_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))

!JR These arrays are mapped to existing buffers to save memory. See top of TimeMarcher subroutine to see mapping.
      !allocate(rhs_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))
      !allocate(dq_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))
      !allocate(qcap_d(1:nx, xstart_gpu(2):xend_gpu(2), xstart_gpu(3):xend_gpu(3)))

!JR Initialize data on GPU 
      vx_d = 0.0_fp_kind
      vy_d = 0.0_fp_kind
      vz_d = 0.0_fp_kind
      pr_d = 0.0_fp_kind
      temp_d = 0.0_fp_kind

      rux_d = 0.0_fp_kind
      ruy_d = 0.0_fp_kind
      ruz_d = 0.0_fp_kind
      rutemp_d = 0.0_fp_kind
      hro_d = 0.0_fp_kind
#else
      allocate(vx_d,source=vx)
      allocate(vy_d,source=vy)
      allocate(vz_d,source=vz)
      allocate(pr_d,source=pr)
      allocate(temp_d,source=temp)

!JR These arrays are mapped to existing buffers to save memory. See top of TimeMarcher subroutine to see mapping.
      !allocate(rhs_d,source=rhs)
      !allocate(dq_d,source=dq)
      !allocate(qcap_d,source=qcap)

      allocate(rux_d,source=rux)
      allocate(ruy_d,source=ruy)
      allocate(ruz_d,source=ruz)
      allocate(hro_d,source=hro)
      allocate(rutemp_d,source=rutemp)

#endif

!JR This array is mapped to existing buffers to save memory. See top of TimeMarcher subroutine to see mapping.
      !allocate(dphhalo_d,source=dphhalo)
      allocate(dph_d,source=dph)
      
!JR   Arrays preallocated here to avoid device synchronization caused by
!     static array allocation
      allocate(amkl_d, source=amkl)
      allocate(apkl_d, source=apkl)
      allocate(ackl_d, source=ackl)
      allocate(fkl_d, source=fkl)

#endif



      return 
      end   

