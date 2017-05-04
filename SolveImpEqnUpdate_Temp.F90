!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Temp_CPU(temp,rhs,am3ssk,ac3ssk,ap3ssk)
  use param,only: fp_kind, al, dt, pec, nxm, nx, lvlhalo, amkl, apkl, ackl, fkl
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
  use ep_solve
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use nvtx
  implicit none

#define ON_CPU
#include "SolveImpEqnUpdate_Temp_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine SolveImpEqnUpdate_Temp_GPU(temp,rhs,am3ssk,ac3ssk,ap3ssk)
  use param,only: fp_kind, al, dt, pec, nxm, nx, lvlhalo, amkl=>amkl_d, apkl=>apkl_d, ackl=>ackl_d, fkl=>fkl_d
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use ep_solve
  use nvtx
  implicit none

#define ON_GPU
#include "SolveImpEqnUpdate_Temp_common.F90"
#undef ON_GPU

end subroutine

#endif


