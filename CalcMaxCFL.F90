!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcMaxCFL.F90                                 !
!    CONTAINS: subroutine CalcMaxCFL                      !
!                                                         ! 
!    PURPOSE: Compute the maximum value of the local CFL  !
!     stability condition for the explicit terms          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcMaxCFL(cflm)
#if defined(USE_CUDA) && defined(USE_HYBRID)
  use cudafor
  use param, only: fp_kind, nxm, dy_h => dy, dy => dy_d, dz_h => dz, dz => dz_d, udx3m_h => udx3m, udx3m => udx3m_d
  use local_arrays, only: vx_h => vx, vx => vx_d, vy_h => vy, vy => vy_d, vz_h => vz, vz => vz_d
#elif defined (USE_CUDA)
  use cudafor
  use param, only: fp_kind, nxm, dy => dy_d, dz => dz_d, udx3m => udx3m_d
  use local_arrays, only: vx => vx_d, vy => vy_d, vz => vz_d
#else
  use param, only: fp_kind, nxm, dy, dz, udx3m
  use local_arrays, only: vx,vy,vz
#endif

#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, xstart_cpu, xend_cpu
#else
  use decomp_2d, only: xstart, xend
#endif
  use mpih
  use nvtx
  implicit none
  real(fp_kind),intent(out)    :: cflm
  integer :: j,k,jp,kp,i,ip
  real(fp_kind) :: qcf

#ifdef USE_HYBRID
  real(fp_kind) :: cflm_h
  real(fp_kind), device :: cflm_d
  cflm_h=real(0.00000001,fp_kind) 
  cflm_d=real(0.00000001,fp_kind)
#endif
  
  cflm=real(0.00000001,fp_kind)

#ifdef USE_CUDA
  !$cuf kernel do(3) <<<*,*>>>
#else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(xstart,xend,nxm,vz,vy,vx) &
  !$OMP   SHARED(dz,dy,udx3m) &
  !$OMP   PRIVATE(i,j,k,ip,jp,kp,qcf) &
  !$OMP   REDUCTION(max:cflm)
#endif
  do i=xstart(3),xend(3)
     ip=i+1
     do j=xstart(2),xend(2)
        jp=j+1
        do k=1,nxm
           kp=k+1
           qcf=( abs((vz(k,j,i)+vz(k,j,ip))*real(0.5,fp_kind)*dz) &
                +abs((vy(k,j,i)+vy(k,jp,i))*real(0.5,fp_kind)*dy) &
                +abs((vx(k,j,i)+vx(kp,j,i))*real(0.5,fp_kind)*udx3m(k)))
           
#ifdef USE_HYBRID
           ! JR Reduce into device resident cflm_d to avoid implicit sync after CUF kernel
           cflm_d = max(cflm_d, qcf) 
#else
           cflm = max(cflm,qcf)
#endif
        enddo
     enddo
  enddo
#ifndef USE_CUDA
  !$OMP END PARALLEL DO
#endif

#ifdef USE_HYBRID
call nvtxStartRangeAsync("CPU",9)
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(xstart_cpu,xend_cpu,nxm,vz_h,vy_h,vx_h) &
  !$OMP   SHARED(dz_h,dy_h,udx3m_h) &
  !$OMP   PRIVATE(i,j,k,ip,jp,kp,qcf) &
  !$OMP   REDUCTION(max:cflm_h)
  do i=xstart_cpu(3),xend_cpu(3)
     ip=i+1
     do j=xstart_cpu(2),xend_cpu(2)
        jp=j+1
        do k=1,nxm
           kp=k+1
           qcf=( abs((vz_h(k,j,i)+vz_h(k,j,ip))*real(0.5,fp_kind)*dz_h) &
                +abs((vy_h(k,j,i)+vy_h(k,jp,i))*real(0.5,fp_kind)*dy_h) &
                +abs((vx_h(k,j,i)+vx_h(kp,j,i))*real(0.5,fp_kind)*udx3m_h(k)))
           
           cflm_h = max(cflm_h,qcf)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
call nvtxEndRangeAsync

  ! Reduce device and host CFL
  cflm = cflm_d 
  cflm = max(cflm, cflm_h)

#endif


  call MpiAllMaxRealScalar(cflm)

      return  
      end
