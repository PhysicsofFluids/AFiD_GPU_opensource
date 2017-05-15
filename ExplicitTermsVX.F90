!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                                         ! 
!    FILE: ExplicitTermsVX.F90                            !
!    CONTAINS: subroutine ExplicitTermsVX                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the x (vertical) dimension          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ExplicitTermsVX)(vx,vy,vz,temp,qcap,udx3c)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz,temp
  real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: qcap
  real(fp_kind), dimension(1:nx), intent(IN) :: udx3c
  #ifdef USE_GPU
  attributes(device)  ::  vx, vy, vz, temp, qcap, udx3c
  #endif

  integer        :: jc,kc
  integer        :: km,kp,jmm,jpp,ic,imm,ipp
  real(fp_kind)  :: hxx,hxy,hxz
  real(fp_kind)  :: udz,udy,tempit
  real(fp_kind)  :: udzq,udyq
  real(fp_kind)  :: dzzvx,dyyvx


  udy=dy*real(0.25,fp_kind)
  udz=dz*real(0.25,fp_kind)

  udyq=dyq/ren
  udzq=dzq/ren

#ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
#endif
  do ic=istart(3), iend(3)
    imm=ic-1
    ipp=ic+1
#ifndef USE_GPU
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,iend,nxm,vz,vy,vx,dz,dy) &
  !$OMP   SHARED(udz, ic, imm, ipp) &
  !$OMP   SHARED(udy,udzq,udyq,udx3c,qcap,temp) &
  !$OMP   PRIVATE(jc,kc,km,kp) &
  !$OMP   PRIVATE(jmm,jpp,tempit) &
  !$OMP   PRIVATE(hxx,hxy,hxz,dzzvx,dyyvx)
#endif
    do jc=istart(2),iend(2)
      jmm=jc-1
      jpp=jc+1
      do kc=2,nxm
        km=kc-1
        kp=kc+1
        !
        !    vx vz term
        !
        !
        !                d  q_x q_t 
        !             -----------
        !                d   t      
        !
        !
        hxz=(((vz(kc,jc,ipp)+vz(km,jc,ipp)) &
        *(vx(kc,jc,ipp)+vx(kc,jc,ic))) &
        -((vz(kc,jc,ic)+vz(km,jc,ic)) &
        *(vx(kc,jc,ic)+vx(kc,jc,imm))))*udz
        !
        !    vx vy term
        !
        !
        !                d  q_x q_r 
        !             -----------
        !                d   r      
        !
        hxy=(((vy(kc,jpp,ic)+vy(km,jpp,ic)) &
        *(vx(kc,jpp,ic)+vx(kc,jc,ic))) &
        -((vy(kc,jc,ic)+vy(km,jc,ic)) &
        *(vx(kc,jc,ic)+vx(kc,jmm,ic))))*udy
        !
        !    vx vx term
        !
        !
        !                 d  q_x q_x 
        !                -----------
        !                 d   x      
        !
        hxx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(vx(kp,jc,ic)+vx(kc,jc,ic)) &
        -(vx(kc,jc,ic)+vx(km,jc,ic))*(vx(kc,jc,ic)+vx(km,jc,ic)) &
        )*udx3c(kc)*real(0.25,fp_kind)
        !
        !  add the buoyancy term
        !
        tempit=temp(kc,jc,ic)
        !
        !   z second derivatives of vx
        !
        dzzvx=(vx(kc,jc,imm) &
             -real(2.0,fp_kind)*vx(kc,jc,ic) &
             +vx(kc,jc,ipp))*udzq
        !
        !   y second derivatives of vx
        !
        dyyvx=(vx(kc,jmm,ic) &
             -real(2.0,fp_kind)*vx(kc,jc,ic) &
             +vx(kc,jpp,ic))*udyq


        qcap(kc,jc,ic) =-(hxx+hxy+hxz)+dyyvx+dzzvx+tempit
          
      enddo
    enddo
#ifndef USE_GPU
    !$OMP  END PARALLEL DO
#endif
  enddo

  return

end subroutine
