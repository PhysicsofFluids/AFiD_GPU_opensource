!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ImplicitAndUpdateVX)(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
  use param, only: fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  use solver_interfaces
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx,pr
  real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3)) :: rhs, rux, qcap
  real(fp_kind), dimension(1:nx), intent(IN) :: am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk
  #ifdef USE_GPU
  attributes(device) :: vx,pr,rhs,rux,qcap,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk
  #endif

  integer :: jc,kc
  integer :: km,kp,ic
  real(fp_kind)    :: alre,udx3,ackl_b,betadx
  real(fp_kind)    :: amm,acc,app,dxp,dxxvx
  real(fp_kind)    :: time1,time2

  alre=al/ren
  betadx=beta*al

#ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
#else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,iend,nxm,vx,pr) &
  !$OMP   SHARED(am3ck,ac3ck,ap3ck) &
  !$OMP   SHARED(al,ga,ro,alre,dt,qcap,betadx) &
  !$OMP   SHARED(udx3c,rhs,rux) &
  !$OMP   PRIVATE(ic,jc,kc,km,kp) &
  !$OMP   PRIVATE(amm,acc,app,udx3) &
  !$OMP   PRIVATE(dxxvx,dxp,ackl_b)
#endif
  do ic=istart(3),iend(3)
    do jc=istart(2),iend(2)
      do kc=2,nxm
        km=kc-1
        kp=kc+1
        udx3 = al*udx3c(kc)
        amm=am3ck(kc)
        acc=ac3ck(kc)
        app=ap3ck(kc)


        ackl_b=real(1.0,fp_kind)/(real(1.0,fp_kind)-ac3ck(kc)*betadx)

        !   Second derivative in x-direction of vx
        !
        dxxvx=vx(kp,jc,ic)*app &
            +vx(kc,jc,ic)*acc &
            +vx(km,jc,ic)*amm

        !  component of grad(pr) along x direction
        !
        dxp=(pr(kc,jc,ic)-pr(km,jc,ic))*udx3

        !    Calculate right hand side of Eq. 5 (VO96)
        !    MF Scaling moved here
        !
        rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*rux(kc,jc,ic) &
                      +alre*dxxvx-dxp)*dt*ackl_b

        !    Store the non-linear terms for the calculation of 
        !    the next timestep

        rux(kc,jc,ic)=qcap(kc,jc,ic)
      enddo
    enddo
  enddo
#ifndef USE_GPU
  !$OMP  END PARALLEL DO
#endif

#ifdef USE_GPU
  !$cuf kernel do(2) <<<*,*>>>
#else
  !$OMP PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,iend,rhs,fp_kind, nx) &
  !$OMP   PRIVATE(ic, jc)
#endif
  do ic=istart(3),iend(3)
    do jc=istart(2),iend(2)
      rhs( 1,jc,ic)=real(0.0,fp_kind)
      rhs(nx,jc,ic)=real(0.0,fp_kind)
    end do
  end do
#ifndef USE_GPU
  !$OMP END PARALLEL DO
#endif


  !  Solve equation and update velocity
  call SolveImpEqnUpdate_X(vx,rhs,am3ssk,ac3ck,ap3ssk)

  !  Set boundary conditions on the vertical velocity at top
  !  and bottom plates. This seems necessary.

  !MF  Using do loops instead of implicit array notation to use CUF kernel on GPU
  !JR  Commenting out reset of BCs, doesn't seem necessary anymore. Leaving just in case.

  !#ifdef ON_GPU
  !!$cuf kernel do(2) <<<*,*>>>
  !#endif
  !do ic=lbound(vx,3),ubound(vx,3)
  !  do jc=lbound(vx,2),ubound(vx,2)
  !    vx(1 ,jc,ic)=real(0.0,fp_kind)
  !    vx(nx,jc,ic)=real(0.0,fp_kind)
  !  end do
  !end do

  return

end subroutine
