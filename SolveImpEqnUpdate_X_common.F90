!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SolveImpEqnUpdate_X_common.F90                 !
!    CONTAINS:                                            !
!                                                         !
!    PURPOSE: Main source code for SolveImpEqnUpdate_X    !
!             subroutines.                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx
real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs
real(fp_kind), dimension(1:nx), intent(IN) :: am3ssk,ac3ck,ap3ssk
#ifdef ON_GPU
attributes(device) :: vx,rhs,am3ssk,ac3ck,ap3ssk
#endif

#ifdef ON_CPU
real(fp_kind) :: amkT(nx-1),apkT(nx-1),appk(nx-2),ackT(nx)
#endif

integer :: jc,kc,info,ic,nrhs,istat
integer :: ipkv(nx)
real(fp_kind) :: betadx,ackl_b

betadx=beta*al
nrhs=(iend(3)-istart(3)+1)*(iend(2)-istart(2)+1)


#ifdef ON_GPU
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

#ifdef ON_GPU
call tepDgtsv_nopivot( nx, nrhs, amkl, ackl, apkl, rhs, nx)      
#endif

#ifdef ON_CPU
call nvtxStartRangeAsync("Solve CPU")
amkT=amkl(2:nx)
apkT=apkl(1:(nx-1))
ackT=ackl(1:nx)

call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
call dgttrs('N',nx,nrhs,amkT,ackT,apkT,appk,ipkv,rhs(1, istart(2), istart(3)),nx,info)
call nvtxEndRangeAsync
#endif

#ifdef ON_GPU
!$cuf kernel do(3) <<<*,*>>>
#endif
#ifdef ON_CPU
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
#ifdef ON_CPU
!$OMP END PARALLEL DO
#endif

#ifdef DEBUG
call compare(rhs,"rhsx")
call compare(vx,"vx")
#endif

return
