!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_X.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_X             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(SolveImpEqnUpdate_X)(vx,rhs,am3ssk,ac3ck,ap3ssk)
#ifdef USE_GPU
  use param,only: fp_kind, nx, nxm, beta, al, lvlhalo, amkl=>amkl_d, apkl=>apkl_d, ackl=>ackl_d
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
  use ep_solve
#else
  use param,only: fp_kind, nx, nxm, beta, al, lvlhalo, amkl, apkl, ackl
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  use nvtx
  implicit none
 
  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx
  real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs
  real(fp_kind), dimension(1:nx), intent(IN) :: am3ssk,ac3ck,ap3ssk
#ifdef USE_GPU
  attributes(device) :: vx,rhs,am3ssk,ac3ck,ap3ssk
#endif

#ifndef USE_GPU
  real(fp_kind) :: amkT(nx-1),apkT(nx-1),appk(nx-2),ackT(nx)
#endif

  integer :: jc,kc,info,ic,nrhs,istat
  integer :: ipkv(nx)
  real(fp_kind) :: betadx,ackl_b

  betadx=beta*al
  nrhs=(iend(3)-istart(3)+1)*(iend(2)-istart(2)+1)


#ifdef USE_GPU
  !$cuf kernel do(1) <<<*,*>>>
#endif
  do kc=1,nx
    ackl_b=real(1.0,fp_kind)/(real(1.0,fp_kind)-ac3ck(kc)*betadx)
    amkl(kc)=-am3ssk(kc)*betadx*ackl_b
    ackl(kc)=real(1.0,fp_kind)
    apkl(kc)=-ap3ssk(kc)*betadx*ackl_b

    if(kc==1 .or. kc==nx) then
      amkl(kc)=real(0.,fp_kind)
      apkl(kc)=real(0.,fp_kind)
      ackl(kc)=real(1.,fp_kind)
    end if
  enddo

#ifdef USE_GPU
  call tepDgtsv_nopivot( nx, nrhs, amkl, ackl, apkl, rhs, nx)      
#else
  call nvtxStartRangeAsync("Solve CPU")
  amkT=amkl(2:nx)
  apkT=apkl(1:(nx-1))
  ackT=ackl(1:nx)

  call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
  call dgttrs('N',nx,nrhs,amkT,ackT,apkT,appk,ipkv,rhs(1, istart(2), istart(3)),nx,info)
  call nvtxEndRangeAsync
#endif

#ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
#else
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(none) &
  !$OMP SHARED(rhs, vx, istart, iend, nxm) &
  !$OMP PRIVATE(ic, jc, kc)
#endif
  do ic=istart(3),iend(3)
    do jc=istart(2),iend(2)
       do kc=2,nxm
        vx(kc,jc,ic)=vx(kc,jc,ic) + rhs(kc,jc,ic)
       end do
    end do
  end do
#ifndef USE_GPU
  !$OMP END PARALLEL DO
#endif

#ifdef DEBUG
  call compare(rhs,"rhsx")
  call compare(vx,"vx")
#endif

  return

end subroutine
