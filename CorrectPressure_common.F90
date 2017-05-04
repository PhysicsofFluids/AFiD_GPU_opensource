!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CorrectPressure_common.F90                     !
!    CONTAINS:                                            !
!                                                         !
!    PURPOSE: Main source code for CorrectPressure        !
!             subroutines.                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(fp_kind), dimension(1:nx ,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: pr
real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: dphhalo
real(fp_kind), dimension(1:nx), intent(IN) :: amphk,acphk,apphk
integer, dimension(1:nx), intent(IN) :: kmv,kpv
#ifdef ON_GPU
attributes(device) :: pr,dphhalo,amphk,acphk,apphk,kmv,kpv
#endif
integer :: kp,km,jm,jp,jc,kc,ic,ip,im
real(fp_kind)    :: be,amm,acc,app

be=al*beta
#ifdef ON_GPU
!$cuf kernel do(3) <<<*,*>>>
#endif
#ifdef ON_CPU
!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(pr,dphhalo,be,amphk,acphk,apphk) &
!$OMP   SHARED(istart,iend,nxm,kmv,kpv,dzq,dyq) &
!$OMP   PRIVATE(ic,jc,kc) &
!$OMP   PRIVATE(im,jm,km,ip,jp,kp) &
!$OMP   PRIVATE(amm,acc,app)
#endif
do ic=istart(3),iend(3)
  im=ic-1
  ip=ic+1
  do jc=istart(2),iend(2)
    jm=jc-1
    jp=jc+1
    do kc=1,nxm
      kp=kpv(kc)
      km=kmv(kc)
      amm=amphk(kc)
      acc=acphk(kc)
      app=apphk(kc)
      pr(kc,jc,ic)=pr(kc,jc,ic)+dphhalo(kc,jc,ic)-be*( &
           (dphhalo(kc,jc,ip) &
           -real(2.0,fp_kind)*dphhalo(kc,jc,ic) &
           +dphhalo(kc,jc,im))*dzq+ &
           (dphhalo(kc,jp,ic) &
           -real(2.0,fp_kind)*dphhalo(kc,jc,ic) &
           +dphhalo(kc,jm,ic))*dyq+ &
           (dphhalo(kp,jc,ic)*app &
           +dphhalo(kc,jc,ic)*acc &
           +dphhalo(km,jc,ic)*amm))
    enddo
  enddo
enddo
#ifdef ON_CPU
!$OMP END PARALLEL DO
#endif
return
