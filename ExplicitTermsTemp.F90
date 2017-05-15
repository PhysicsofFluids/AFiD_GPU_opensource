!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsTemp.F90                          !
!    CONTAINS: subroutine ExplicitTermsTemp               !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(ExplicitTermsTemp)(vx,vy,vz,temp,hro,udx3c)
  use param, only: fp_kind, nxm, dy , dz, pec, dzq, dyq, nx, lvlhalo
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz,temp
  real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: hro
  real(fp_kind), dimension(1:nx), intent(IN) :: udx3c
  #ifdef USE_GPU
  attributes(device) :: vx,vy,vz,temp,hro,udx3c
  #endif

  integer :: jc,kc,ic
  integer :: km,kp,jm,jp,im,ip
  real(fp_kind)    :: htx,hty,htz,udy,udz
  real(fp_kind)    :: udzq,udyq
  real(fp_kind)    :: dyyt,dzzt

  udz=dz*real(0.25,fp_kind)
  udy=dy*real(0.25,fp_kind)
  udzq=dzq/pec
  udyq=dyq/pec

  #ifdef USE_GPU
  !$cuf kernel do(3) <<<*,*>>>
  #endif
  do ic=istart(3),iend(3)
    im=ic-1
    ip=ic+1
  #ifndef USE_GPU
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(istart,iend,vz,vy,vx,nxm) &
  !$OMP   SHARED(udz) &
  !$OMP   SHARED(ic,im,ip) &
  !$OMP   SHARED(udy,udzq,udyq,udx3c,temp,hro) &
  !$OMP   PRIVATE(jc,kc,km,kp,jm,jp) &
  !$OMP   PRIVATE(htx,hty,htz,dyyt,dzzt)
  #endif
    do jc=istart(2),iend(2)
      jm=jc-1
      jp=jc+1
      do kc=2,nxm
        km=kc-1
        kp=kc+1
        !
        !
        !    rho vz term
        !
        !
        !                d  rho q_z
        !             -----------
        !                d   z      
        !
        htz=((vz(km,jc,ip)+vz(kc,jc,ip))*(temp(kc,jc,ip)+temp(kc,jc,ic))- &
        (vz(km,jc,ic)+vz(kc,jc,ic))*(temp(kc,jc,ic)+temp(kc,jc,im)) &
        )*udz
        !
        !
        !    rho vy term
        !
        !
        !                d  rho q_y 
        !             -----------
        !                d   y      
        !
        hty=((vy(kc,jp,ic)+vy(km,jp,ic))*(temp(kc,jp,ic)+temp(kc,jc,ic))- &
        (vy(kc,jc,ic)+vy(km,jc,ic))*(temp(kc,jc,ic)+temp(kc,jm,ic)) &
        )*udy
        !
        !    rho vx term
        !
        !
        !                 d  rho q_x 
        !                -----------
        !                 d   x      
        !
        htx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(temp(kp,jc,ic)+temp(kc,jc,ic))- &
        (vx(kc,jc,ic)+vx(km,jc,ic))*(temp(kc,jc,ic)+temp(km,jc,ic)) &
        )*udx3c(kc)*real(0.25,fp_kind)
        !
        !
        !   zz second derivatives of temp
        !
        dzzt=(temp(kc,jc,ip) &
         -real(2.0,fp_kind)*temp(kc,jc,ic) &
             +temp(kc,jc,im))*udzq

        !
        !   yy second derivatives of temp
        !
        dyyt=(temp(kc,jp,ic) &
         -real(2.0,fp_kind)*temp(kc,jc,ic) &
             +temp(kc,jm,ic))*udyq
        !
        hro(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt
      enddo
    enddo
  #ifndef USE_GPU
  !$OMP  END PARALLEL DO
  #endif
  enddo

  return

end subroutine

