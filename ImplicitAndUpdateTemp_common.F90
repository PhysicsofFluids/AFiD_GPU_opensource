!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ImplicitAndUpdateTemp_common.F90               !
!    CONTAINS:                                            !
!                                                         !
!    PURPOSE: Main source code for ImplicitAndUpdateTemp  !
!             subroutines.                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: temp
real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, rutemp, hro
real(fp_kind), dimension(1:ny,1:nz) :: tempbp, temptp
real(fp_kind), dimension(1:nx), intent(IN) :: am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk
#ifdef ON_GPU
attributes(device) :: temp,rhs,rutemp,hro,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk
#endif

integer :: jc,kc,ic
integer :: km,kp
real(fp_kind)    :: alpec,dxxt,betadx,ackl_b
real(fp_kind)    :: app,acc,amm

alpec=al/pec
betadx=real(0.5,fp_kind)*al*dt/pec

#ifdef ON_GPU
!$cuf kernel do(3) <<<*,*>>>
#endif
#ifdef ON_CPU
!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(istart,iend,nxm,temp) &
!$OMP   SHARED(am3ck,ac3ck,ap3ck,ac3ssk) &
!$OMP   SHARED(ga,ro,alpec,dt,betadx) &
!$OMP   SHARED(rhs,rutemp,hro) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app,ackl_b) &
!$OMP   PRIVATE(dxxt)
#endif
do ic=istart(3),iend(3)
  do jc=istart(2),iend(2)
    do kc=2,nxm
      ackl_b=real(1.0,fp_kind)/(real(1.0,fp_kind)-ac3ssk(kc)*betadx)

      !   Calculate second derivative of temperature in the x-direction.
      !   This is the only term calculated implicitly for temperature.

      dxxt= temp(kc+1,jc,ic)*ap3ck(kc) &
           +temp(kc  ,jc,ic)*ac3ck(kc) &
           +temp(kc-1,jc,ic)*am3ck(kc)


      !    Calculate right hand side of Eq. 5 (VO96)
      !    Normalize by ackl_b

      rhs(kc,jc,ic)=(ga*hro(kc,jc,ic)+ro*rutemp(kc,jc,ic) &
              +alpec*dxxt)*dt*ackl_b

      !    Store the non-linear terms for the calculation of 
      !    the next timestep

      rutemp(kc,jc,ic)=hro(kc,jc,ic)

    enddo
  enddo
enddo
#ifdef ON_CPU
!$OMP  END PARALLEL DO
#endif

#ifdef ON_GPU
!$cuf kernel do(2) <<<*,*>>>
#endif
do ic=istart(3),iend(3)
  do jc=istart(2),iend(2)
    rhs( 1,jc,ic)=real(0.0,fp_kind)
    rhs(nx,jc,ic)=real(0.0,fp_kind)
  end do
end do

!  Solve equation and update temperature
call SolveImpEqnUpdate_Temp(temp,rhs,am3ssk,ac3ssk,ap3ssk)

!  Set boundary conditions on the temperature field at top
!  and bottom plates. This seems necessary.

!MF  Using do loops instead of implicit array notation to use CUF kernel on GPU
!JR  Commenting out reset of BCs, doesn't seem necessary anymore. Leaving just in case.

!#ifdef ON_GPU
!!$cuf kernel do(2) <<<*,*>>>
!#endif
!do ic=istart(3),iend(3)
!  do jc=istart(2),iend(2)
!    temp(1 ,jc,ic) = tempbp(jc,ic)
!    temp(nx,jc,ic) = temptp(jc,ic)
!  end do 
!end do 

return
