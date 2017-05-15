!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_YZ.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_YZ            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any of the horizontal directions, and updates    !
!     it to time t+dt                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(SolveImpEqnUpdate_YZ)(q,rhs,am3sk,ac3sk,ap3sk)
#ifdef USE_GPU
  use param, only: fp_kind, al, beta, nxm, nx, lvlhalo, amkl=>amkl_d, apkl=>apkl_d, ackl=>ackl_d, fkl=>fkl_d
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
  use ep_solve
#else
  use param, only: fp_kind, al, beta, nxm, nx, lvlhalo, amkl, apkl, ackl, fkl
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  use nvtx
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: q
  real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)) :: rhs
  real(fp_kind), dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
  #ifdef USE_GPU
  attributes(device) :: am3sk,ac3sk,ap3sk,q,rhs
  #endif



  integer :: jc,kc,ic,nrhs,info,ipkv(nxm)
  real(fp_kind) :: betadx,ackl_b

  #ifndef USE_GPU
  real(fp_kind) :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nx-3)
  #endif

  betadx=beta*al
  nrhs=(iend(3)-istart(3)+1)*(iend(2)-istart(2)+1)

  #ifdef USE_GPU
  !$cuf kernel do(1) <<<*,*>>>
  #endif
  do kc=1,nxm
    ackl_b=real(1.0,fp_kind)/(real(1.0,fp_kind)-ac3sk(kc)*betadx)
    amkl(kc)=-am3sk(kc)*betadx*ackl_b
    ackl(kc)=real(1.0,fp_kind)
    apkl(kc)=-ap3sk(kc)*betadx*ackl_b
  enddo

  #ifdef USE_GPU
  call  tepDgtsv_nopivot( nxm, nrhs, amkl, ackl, apkl, rhs, nx) 
  ! Really slow way until we get a proper solver
  !do ic=xstart(3),xend(3)
  !  do jc=xstart(2),xend(2)
  !    istat = cusparseDgtsv_nopivot(cusparseH, nxm, 1, amkl, ackl, apkl, rhs(1,jc,ic), nx) 
  ! end do
  !end do

  #else

  call nvtxStartRangeAsync("Solve CPU")
  !#ifndef USE_HYBRID
  amkT=amkl(2:nxm)
  apkT=apkl(1:(nxm-1))
  ackT=ackl(1:nxm)

  call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
  call dgttrs('N',nxm,nrhs,amkT,ackT,apkT,appk,ipkv,rhs(1, istart(2), istart(3)),nx,info)
  !call dgttrs('N',nxm,nrhs,amkT,ackT,apkT,appk,ipkv,rhs,nx,info)
  !#else
  !!call  tepDgtsv_nopivot_cpu( nxm, nrhs, amkl, ackl, apkl, rhs(1, istart(2), istart(3)), nx) 
  !call  epDgtsv_nopivot_cpu( nxm, nrhs, amkl, ackl, apkl, rhs(1, istart(2), istart(3)), nx) 
  !#endif
  call nvtxEndRangeAsync
  #endif

  #ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
  #else
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(none) &
  !$OMP SHARED(rhs, q, istart, iend, nxm) &
  !$OMP PRIVATE(ic, jc, kc)
  #endif
  do ic=istart(3),iend(3)
    do jc=istart(2),iend(2)
      do kc=1,nxm
        q(kc,jc,ic)=q(kc,jc,ic) + rhs(kc,jc,ic)
      end do
    end do
  end do
  #ifndef USE_GPU
  !$OMP END PARALLEL DO
  #endif

  ! print *,"                     after Solve min(q)", minval(q)
  ! print *,"                                 max(q)", maxval(q)

  return

end subroutine
