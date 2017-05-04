!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateTemp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateTemp           !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and call the implicit solver.       !
!     After this routine, the temperature has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateTemp_CPU(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
  use param, only : fp_kind, al, pec, nxm, ga, ro, dt, nx, ny, nz, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use solver_interfaces
  implicit none

#define ON_CPU
#include "ImplicitAndUpdateTemp_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine ImplicitAndUpdateTemp_GPU(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
  use param, only : fp_kind, al, pec, nxm, ga, ro, dt, nx, ny, nz, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  use solver_interfaces
  implicit none

#define ON_GPU
#include "ImplicitAndUpdateTemp_common.F90"
#undef ON_GPU

end subroutine
#endif


