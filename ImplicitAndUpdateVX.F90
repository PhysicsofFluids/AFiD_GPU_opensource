!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVX_CPU(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
  use param, only: fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use solver_interfaces
  implicit none

#define ON_CPU
#include "ImplicitAndUpdateVX_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine ImplicitAndUpdateVX_GPU(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
  use param, only: fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use solver_interfaces
  implicit none

#define ON_GPU
#include "ImplicitAndUpdateVX_common.F90"
#undef ON_GPU

end subroutine
#endif

