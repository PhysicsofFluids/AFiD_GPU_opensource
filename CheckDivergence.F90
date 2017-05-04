!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CheckDivergence.F90                            !
!    CONTAINS: subroutine CheckDivergence                 !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckDivergence(qmax)

#if defined(USE_CUDA) && defined(USE_HYBRID)
  use cudafor
  use param, only: fp_kind, nxm,dz,dy, udx3m_h => udx3m, udx3m => udx3m_d
  use local_arrays, only: vy_h => vy, vy => vy_d, vx_h => vx, vx => vx_d, vz_h => vz, vz => vz_d
#elif defined (USE_CUDA)
  use cudafor
  use param, only: fp_kind, nxm,dz,dy, udx3m => udx3m_d
  use local_arrays, only: vy => vy_d, vx => vx_d, vz => vz_d
#else
  use param, only: fp_kind, nxm,dz,dy, udx3m
  use local_arrays, only: vy,vx,vz
#endif
  use mpih
#ifdef USE_HYBRID
  use decomp_2d, only: xstart => xstart_gpu, xend => xend_gpu, xstart_cpu, xend_cpu, nrank
#else
  use decomp_2d, only: xstart,xend, nrank
#endif
  use nvtx
  implicit none
  real(fp_kind),intent(out) :: qmax
  integer :: jc,kc,kp,jp,ic,ip, istat
  real(fp_kind)    :: dqcap
#ifdef USE_HYBRID
  real(fp_kind) :: qmax_h
  real(fp_kind), device :: qmax_d
  qmax_h =-huge(real(0.0,fp_kind))
  qmax_d =-huge(real(0.0,fp_kind))
#endif

  qmax =-huge(real(0.0,fp_kind))

#ifdef USE_CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy,udx3m) &
  !$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
  !$OMP   PRIVATE(dqcap) &
  !$OMP   REDUCTION(max:qmax)
#endif
  do ic=xstart(3),xend(3)
     ip=ic+1
     do jc=xstart(2),xend(2)
        jp=jc+1
        do kc=1,nxm
           kp=kc+1
           dqcap= (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
                 +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
                 +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
#ifdef USE_HYBRID
           qmax_d = max(abs(dqcap),qmax_d)          
#else
           qmax = max(abs(dqcap),qmax)          
#endif
        enddo
     enddo
  enddo
#ifndef USE_CUDA
  !$OMP END PARALLEL DO
#endif

#ifdef USE_HYBRID
  call nvtxStartRangeAsync("CPU", 4)
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(xstart_cpu,xend_cpu,nxm,vz_h,vy_h,vx_h,dz,dy,udx3m_h) &
  !$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
  !$OMP   PRIVATE(dqcap) &
  !$OMP   REDUCTION(max:qmax_h)
  do ic=xstart_cpu(3),xend_cpu(3)
     ip=ic+1
     do jc=xstart_cpu(2),xend_cpu(2)
        jp=jc+1
        do kc=1,nxm
           kp=kc+1
           dqcap= (vz_h(kc,jc,ip)-vz_h(kc,jc,ic))*dz &
                 +(vy_h(kc,jp,ic)-vy_h(kc,jc,ic))*dy &
                 +(vx_h(kp,jc,ic)-vx_h(kc,jc,ic))*udx3m_h(kc)
           qmax_h = max(abs(dqcap),qmax_h)          
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  call nvtxEndRangeAsync

  ! Reduce device and host divergence
  qmax = qmax_d
  qmax = max(qmax, qmax_h)
#endif

  call MpiMaxRealScalar(qmax)
  
  return     
end subroutine CheckDivergence
