!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectVelocity_CPU(vx,vy,vz,dphhalo,udx3c,kmv)
  use param,only: fp_kind, al, dt, dy, dz, nxm, lvlhalo, nx
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_CPU
#include "CorrectVelocity_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine CorrectVelocity_GPU(vx,vy,vz,dphhalo,udx3c,kmv)
  use param,only: fp_kind, al, dt, dy, dz, nxm, lvlhalo, nx
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_GPU
#include "CorrectVelocity_common.F90"
#undef ON_GPU

end subroutine
#endif


