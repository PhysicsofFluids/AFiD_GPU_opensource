!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVX.F90                            !
!    CONTAINS: subroutine ExplicitTermsVX                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the x (vertical) dimension          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVX_CPU(vx,vy,vz,temp,qcap,udx3c)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_CPU
#include "ExplicitTermsVX_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine ExplicitTermsVX_GPU(vx,vy,vz,temp,qcap,udx3c)
  use param, only: fp_kind, ren, nxm, dy, dz, dzq, dyq, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_GPU
#include "ExplicitTermsVX_common.F90"
#undef ON_GPU

end subroutine
#endif


