!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcLocalDivergence.F90                        !
!    CONTAINS: subroutine CalcLocalDivergence             !
!                                                         ! 
!    PURPOSE: Compute the divergence of the intermediate  !
!     velocity at every point for the pressure            !
!     correction step                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcLocalDivergence_CPU(vx,vy,vz,dph,udx3m)
  use param,only:  fp_kind, dt, al, nxm, dz, dy, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_CPU
#include "CalcLocalDivergence_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine CalcLocalDivergence_GPU(vx,vy,vz,dph,udx3m)
  use param,only:  fp_kind, dt, al, nxm, dz, dy, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_GPU
#include "CalcLocalDivergence_common.F90"
#undef ON_GPU

end subroutine
#endif


