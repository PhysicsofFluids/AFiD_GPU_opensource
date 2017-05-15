!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ImplicitAndUpdateVZ)(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
  use param, only : fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, dz, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  use solver_interfaces
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vz,pr
  real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, ruz, dq
  real(fp_kind), dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
  integer, dimension(1:nx), intent(IN) :: kmv,kpv
  #ifdef USE_GPU
  attributes(device) :: vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv
  #endif
  integer :: kc,jc,ic,imm
  integer :: kmm,kpp
  real(fp_kind)    :: alre,amm,acc,app,udz,betadx,ackl_b
  real(fp_kind)    :: dxxvz,dzp
  real(fp_kind)    :: time1,time2

  alre=al/ren
  betadx=beta*al
  udz=dz*al

  #ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
  #else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,iend,nxm,vz,pr) &
  !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
  !$OMP   SHARED(dz,al,ga,ro,alre,dt,dq,betadx) &
  !$OMP   SHARED(rhs,ruz) &
  !$OMP   PRIVATE(ic,jc,kc,imm,kmm,kpp) &
  !$OMP   PRIVATE(amm,acc,app) &
  !$OMP   PRIVATE(dxxvz,dzp,ackl_b)
  #endif
  do ic=istart(3),iend(3)
    imm=ic-1
    do jc=istart(2),iend(2)
      do kc=1,nxm
        kmm=kmv(kc)
        kpp=kpv(kc)
        amm=am3sk(kc)
        acc=ac3sk(kc)
        app=ap3sk(kc)

        ackl_b=real(1.0,fp_kind)/(real(1.0,fp_kind)-ac3sk(kc)*betadx)      

        !   Second derivative in x-direction of vz
        !
        dxxvz=vz(kpp,jc,ic)*app &
             +vz(kc ,jc,ic)*acc &
             +vz(kmm,jc,ic)*amm

        !   component of grad(pr) along z direction
        !
        dzp=(pr(kc,jc,ic)-pr(kc,jc,imm))*dz*al

        !    Calculate right hand side of Eq. 5 (VO96)
        !
        rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ruz(kc,jc,ic) &
                      +alre*dxxvz-dzp)*dt*ackl_b

        !    Store the non-linear terms for the calculation of 
        !    the next timestep

        ruz(kc,jc,ic)=dq(kc,jc,ic)
      enddo
    enddo
  enddo
  #ifndef USE_GPU
  !$OMP END PARALLEL DO
  #endif

  !  Solve equation and update velocity
  call SolveImpEqnUpdate_YZ(vz,rhs,am3sk,ac3sk,ap3sk)

  return

end subroutine
