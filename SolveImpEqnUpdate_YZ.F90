!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_YZ.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_YZ            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any of the horizontal directions, and updates    !
!     it to time t+dt                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_YZ_CPU(q,rhs,am3sk,ac3sk,ap3sk)
  use param, only: fp_kind, al, beta, nxm, nx, lvlhalo, amkl, apkl, ackl, fkl
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
  use ep_solve
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use nvtx
  implicit none

#define ON_CPU
#include "SolveImpEqnUpdate_YZ_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine SolveImpEqnUpdate_YZ_GPU(q,rhs,am3sk,ac3sk,ap3sk)
  use param, only: fp_kind, al, beta, nxm, nx, lvlhalo, amkl=>amkl_d, apkl=>apkl_d, ackl=>ackl_d, fkl=>fkl_d
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use ep_solve
  use nvtx
  implicit none

#define ON_GPU
#include "SolveImpEqnUpdate_YZ_common.F90"
#undef ON_GPU

end subroutine

#endif


