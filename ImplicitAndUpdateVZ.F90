!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVZ_CPU(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
  use param, only : fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, dz, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use solver_interfaces
  implicit none

#define ON_CPU
#include "ImplicitAndUpdateVZ_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine ImplicitAndUpdateVZ_GPU(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
  use param, only : fp_kind, al, beta, ren, nxm, ga, ro, dt, nx, dz, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use solver_interfaces
  implicit none

#define ON_GPU
#include "ImplicitAndUpdateVZ_common.F90"
#undef ON_GPU

end subroutine
#endif


