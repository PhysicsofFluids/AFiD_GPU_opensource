!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_X.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_X             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_X_CPU(vx,rhs,am3ssk,ac3ck,ap3ssk)
  use param,only: fp_kind, nx, nxm, beta, al, lvlhalo, amkl, apkl, ackl
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
  use ep_solve
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use nvtx
  implicit none
 
#define ON_CPU
#include "SolveImpEqnUpdate_X_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine SolveImpEqnUpdate_X_GPU(vx,rhs,am3ssk,ac3ck,ap3ssk)
  use param,only: fp_kind, nx, nxm, beta, al, lvlhalo, amkl=>amkl_d, apkl=>apkl_d, ackl=>ackl_d
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use ep_solve
  use nvtx
  implicit none

#define ON_GPU
#include "SolveImpEqnUpdate_X_common.F90"
#undef ON_GPU

end subroutine

#endif


