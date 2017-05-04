!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVY.F90                            !
!    CONTAINS: subroutine ExplicitTermsVY                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the y (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ExplicitTermsVY_CPU(vx,vy,vz,dph,udx3m,kmv,kpv)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif

  implicit none

#define ON_CPU
#include "ExplicitTermsVY_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine ExplicitTermsVY_GPU(vx,vy,vz,dph,udx3m,kmv,kpv)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif

  implicit none

#define ON_GPU
#include "ExplicitTermsVY_common.F90"
#undef ON_GPU

end subroutine
#endif
