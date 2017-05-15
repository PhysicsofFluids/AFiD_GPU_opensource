!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVZ.F90                            !
!    CONTAINS: subroutine ExplicitTermsVZ                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ExplicitTermsVZ)(vx,vy,vz,dq,udx3m,kmv,kpv)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif

  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
  real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: dq
  real(fp_kind), dimension(1:nx), intent(IN) :: udx3m
  integer, dimension(1:nx), intent(IN) :: kmv,kpv
  #ifdef USE_GPU
  attributes(device) :: vx,vy,vz,dq,udx3m,kmv,kpv
  #endif

  integer          :: kc,kp,jpp,jmm,jc,ic,imm,ipp
  integer          :: kmm,kpp
  real(fp_kind)    :: hzx,hzy,hzz,udy,udz
  real(fp_kind)    :: udyq,udzq
  real(fp_kind)    :: dzzvz,dyyvz

  udyq=dyq/ren
  udzq=dzq/ren

  udy=dy*real(0.25,fp_kind)
  udz=dz*real(0.25,fp_kind)

  #ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
  #endif
  do ic=istart(3),iend(3)
    imm=ic-1
    ipp=ic+1
  #ifndef USE_GPU
  !$OMP  PARALLEL DO &
  !$OMP  DEFAULT(none) &
  !$OMP  SHARED(istart,iend,vz,vy,vx,dz,dy,udx3m) &
  !$OMP  SHARED(kmv,kpv,udz) &
  !$OMP  SHARED(ic,imm,ipp) &
  !$OMP  SHARED(udy,udzq,udyq,dq,nxm) &
  !$OMP  PRIVATE(jc,kc,kmm,kp,kpp) &
  !$OMP  PRIVATE(jmm,jpp) &
  !$OMP  PRIVATE(hzz,hzy,hzx,dzzvz,dyyvz)
  #endif
    do jc=istart(2),iend(2)
      jmm=jc-1
      jpp=jc+1
      do kc=1,nxm
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1

        !     vz vz term
        !
        !
        !                 d  q_t q_t 
        !                ------------
        !                 d   t      
        !
        hzz=( (vz(kc,jc,ipp)+vz(kc,jc,ic)) &
         *(vz(kc,jc,ipp)+vz(kc,jc,ic)) &
         -(vz(kc,jc,imm)+vz(kc,jc,ic)) &
         *(vz(kc,jc,imm)+vz(kc,jc,ic)) &
        )*udz

        !     vz vy term
        !
        !
        !                 d  q_t q_r 
        !                ------------
        !                 d   r      
        !
        hzy=( (vy(kc,jpp,ic)+vy(kc,jpp,imm)) &
         *(vz(kc,jpp,ic)+vz(kc,jc,ic)) &
         -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
         *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
        )*udy
        !
        !     vz vx term
        !
        !
        !                 d  q_t q_x 
        !                -----------
        !                 d   x      
        !
        hzx=((vx(kp,jc,ic)+vx(kp,jc,imm))*(vz(kpp,jc,ic)+vz(kc,jc,ic)) &
        -(vx(kc,jc,ic)+vx(kc,jc,imm))*(vz(kc,jc,ic)+vz(kmm,jc,ic)) &
        )*udx3m(kc)*real(0.25,fp_kind)
        !
        !
        !
        !   11 second derivative of vz
        !
        dzzvz=(vz(kc,jc,ipp) &
              -2.0*vz(kc,jc,ic) &
              +vz(kc,jc,imm))*udzq
        !
        !   22 second derivative of vz
        !
        dyyvz=(vz(kc,jpp,ic) &
              -2.0*vz(kc,jc,ic) &
              +vz(kc,jmm,ic))*udyq

        !
        dq(kc,jc,ic)=-(hzx+hzy+hzz)+dyyvz+dzzvz
        !
      enddo
    enddo
  #ifndef USE_GPU
  !$OMP END PARALLEL DO
  #endif
  enddo


  return

end subroutine
