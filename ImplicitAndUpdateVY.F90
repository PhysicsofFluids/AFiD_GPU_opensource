!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVY.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVY             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ImplicitAndUpdateVY)(vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv)
  use param, only: fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, nxm, dy, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  use solver_interfaces
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vy,pr
  real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, ruy
  real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(IN)  :: dph
  real(fp_kind), dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
  integer, dimension(1:nx), intent(IN) :: kmv,kpv
  #ifdef USE_GPU
  attributes(device) :: vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv
  #endif
  integer :: kc,jmm,jc,ic
  integer :: kpp,kmm
  real(fp_kind)    :: alre,udy,betadx,ackl_b
  real(fp_kind)    :: amm,acc,app
  real(fp_kind)    :: dyp,dxxvy
  real(fp_kind)    :: time1,time2


  alre=al/ren
  betadx=beta*al
  udy=dy*al

#ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
#else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,iend,nxm,vy,pr) &
  !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
  !$OMP   SHARED(dy,al,ga,ro,alre,dt,dph,betadx) &
  !$OMP   SHARED(udy,rhs,ruy) &
  !$OMP   PRIVATE(ic,jc,kc,kmm,kpp,jmm) &
  !$OMP   PRIVATE(amm,acc,app) &
  !$OMP   PRIVATE(dyp,dxxvy,ackl_b)
#endif
  do ic=istart(3),iend(3)
    do jc=istart(2),iend(2)
      jmm=jc-1
      do kc=1,nxm
        kmm=kmv(kc)
        kpp=kpv(kc)
        amm=am3sk(kc)
        acc=ac3sk(kc)
        app=ap3sk(kc)

        ackl_b=real(1.0,fp_kind)/(real(1.0,fp_kind)-ac3sk(kc)*betadx)


        !   Second derivative in x-direction of vy
        !
        !
        dxxvy=vy(kpp,jc,ic)*app &
             +vy(kc,jc,ic)*acc &
             +vy(kmm,jc,ic)*amm

        !   component of grad(pr) along y direction
        !
        dyp=(pr(kc,jc,ic)-pr(kc,jmm,ic))*udy

        !    Calculate right hand side of Eq. 5 (VO96)
        !    Normalize 
        rhs(kc,jc,ic)=(ga*dph(kc,jc,ic)+ro*ruy(kc,jc,ic) &
                      +alre*dxxvy-dyp)*dt*ackl_b

        !    Store the non-linear terms for the calculation of 
        !    the next timestep

        ruy(kc,jc,ic)=dph(kc,jc,ic)
      enddo
    enddo
  enddo
#ifndef USE_GPU
  !$OMP  END PARALLEL DO
#endif

  !  Solve equation and update velocity
  call SolveImpEqnUpdate_YZ(vy,rhs,am3sk,ac3sk,ap3sk)
  return

end subroutine
