!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVY.F90                            !
!    CONTAINS: subroutine ExplicitTermsVY                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the y (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ExplicitTermsVY)(vx,vy,vz,dph,udx3m,kmv,kpv)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif

  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
  real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(out)  :: dph
  real(fp_kind), dimension(1:nx), intent(IN) :: udx3m
  integer, dimension(1:nx), intent(IN) :: kmv,kpv
  #ifdef USE_GPU
  attributes(device) :: vx,vy,vz,dph,udx3m,kmv,kpv
  #endif

  integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
  integer :: kpp,kmm
  real(fp_kind)    :: udzq,udyq
  real(fp_kind)    :: udy,udz,hyx,hyy,hyz 
  real(fp_kind)    :: dyyvy, dzzvy


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
  !$OMP  SHARED(istart,iend,vz,vy,vx,dz,dy) &
  !$OMP  SHARED(kmv,kpv,udz) &
  !$OMP  SHARED(udy,udzq,udyq,udx3m,dph,nxm) &
  !$OMP  SHARED(ic,imm,ipp) &
  !$OMP  PRIVATE(jc,kc,kmm,kp,kpp) &
  !$OMP  PRIVATE(jmm,jpp) &
  !$OMP  PRIVATE(hyx,hyy,hyz,dyyvy,dzzvy)
  #endif 
    do jc=istart(2),iend(2)
      jmm=jc-1
      jpp=jc+1
      do kc=1,nxm
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1
        !
        !     vy vx term
        !
        !                 d  q_x q_r 
        !                -----------
        !                 d   x      
        !
        hyx=((vx(kp,jc,ic)+vx(kp,jmm,ic))*(vy(kpp,jc,ic)+vy(kc,jc,ic)) &
        -(vx(kc,jc,ic)+vx(kc,jmm,ic))*(vy(kc,jc,ic)+vy(kmm,jc,ic)) &
        )*udx3m(kc)*real(0.25,fp_kind)


        !     
        !     vy vy term
        !
        !                 d  q_r q_r 
        !                ------------
        !                 d   r      
        !
        hyy=( (vy(kc,jpp,ic)+vy(kc,jc,ic)) &
         *(vy(kc,jpp,ic)+vy(kc,jc,ic)) &
         -(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
         *(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
        )*udy
        !
        !     vz vy term
        !
        !                 d  q_t q_r 
        !                ------------
        !                 d   t      
        !
        hyz=( (vy(kc,jc,ipp)+vy(kc,jc,ic)) &
         *(vz(kc,jc,ipp)+vz(kc,jmm,ipp)) &
         -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
         *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
        )*udz

        !
        !   yy second derivative of vy
        !
        dyyvy=(vy(kc,jpp,ic) &
              -real(2.0,fp_kind)*vy(kc,jc,ic) &
              +vy(kc,jmm,ic))*udyq

        !
        !   zz second derivative of vy
        !
        dzzvy=(vy(kc,jc,ipp) &
              -real(2.0,fp_kind)*vy(kc,jc,ic) &
              +vy(kc,jc,imm))*udzq


        dph(kc,jc,ic)=-(hyx+hyy+hyz)+dyyvy+dzzvy
      enddo
    enddo
  #ifndef USE_GPU
    !$OMP  END PARALLEL DO
  #endif
  enddo

  return

end subroutine
