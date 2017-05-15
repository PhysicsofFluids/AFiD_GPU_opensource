!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_GPU
#else
#define MY_ROUTINE(x)  x##_CPU
#endif

subroutine MY_ROUTINE(CorrectVelocity)(vx,vy,vz,dphhalo,udx3c,kmv)
  use param,only: fp_kind, al, dt, dy, dz, nxm, lvlhalo, nx
#ifdef USE_GPU
  use decomp_2d, only: xstart, xend, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#endif
  implicit none

  real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx,vy,vz
  real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: dphhalo
  real(fp_kind), dimension(1:nx), intent(IN) :: udx3c
  integer, dimension(1:nx), intent(IN) :: kmv
#ifdef USE_GPU
  attributes(device) :: vx,vy,vz,dphhalo,udx3c,kmv
#endif

  integer :: jc,jm,kc,km,ic,im
  real(fp_kind)    :: usukm,udy,udz,locdph

  udy = al*dt*dy
  udz = al*dt*dz

#ifdef USE_GPU
  !$cuf kernel do (3) <<<*,*>>>
#else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(vz,vy,vx,dphhalo,udz,udy,udx3c) &
  !$OMP   SHARED(istart,iend,nxm,kmv,dt,al) &
  !$OMP   PRIVATE(ic,jc,kc) &
  !$OMP   PRIVATE(im,jm,km,usukm,locdph)
#endif
  do ic=istart(3),iend(3)
    im=ic-1
    do jc=istart(2),iend(2)
      jm=jc-1
      do kc=1,nxm
        km=kmv(kc)
        usukm = al*dt*udx3c(kc)
        locdph=dphhalo(kc,jc,ic)
        vx(kc,jc,ic)=vx(kc,jc,ic)- &
             (locdph-dphhalo(km,jc,ic))*usukm
        vy(kc,jc,ic)=vy(kc,jc,ic)- &
             (locdph-dphhalo(kc,jm,ic))*udy
        vz(kc,jc,ic)=vz(kc,jc,ic)- &
             (locdph-dphhalo(kc,jc,im))*udz
      enddo
    enddo
  enddo
#ifndef USE_GPU
  !$OMP END PARALLEL DO
#endif
  return

end subroutine
