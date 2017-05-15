!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcLocalDivergence.F90                        !
!    CONTAINS: subroutine CalcLocalDivergence             !
!                                                         ! 
!    PURPOSE: Compute the divergence of the intermediate  !
!     velocity at every point for the pressure            !
!     correction step                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(CalcLocalDivergence)(vx,vy,vz,dph,udx3m)
  use param,only:  fp_kind, dt, al, nxm, dz, dy, nx, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
  real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3))  :: dph
  real(fp_kind), dimension(1:nx), intent(IN) :: udx3m

  #ifdef USE_GPU
  attributes(device) :: vx,vy,vz,dph,udx3m
  #endif

  integer :: jc,jp,kc,kp,ic,ip
  real(fp_kind)    :: usdtal,dqcap   

  usdtal = real(1.0,fp_kind)/(dt*al)

  #ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
  #endif
  do ic=istart(3),iend(3)
    ip=ic+1
  #ifndef USE_GPU
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,ic,ip,vz,vy,vx,dz,dy,udx3m,usdtal) &
  !$OMP   SHARED(dph,nxm,iend) &
  !$OMP   PRIVATE(jc,kc,jp,kp) &
  !$OMP   PRIVATE(dqcap)
  #endif
    do jc=istart(2),iend(2)
      jp=jc+1
      do kc=1,nxm
        kp=kc+1
        dqcap= (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
              +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
              +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
        dph(kc,jc,ic)=dqcap*usdtal
      enddo
    enddo
  #ifndef USE_GPU
  !$OMP END PARALLEL DO
  #endif
  enddo

  return

end subroutine
