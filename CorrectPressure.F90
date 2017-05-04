!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectPressure.F90                            !
!    CONTAINS: subroutine CorrectPressure                 !
!                                                         ! 
!    PURPOSE: Apply the pressure correction to the        !
!     pressure                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectPressure_CPU(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
  use param,only: fp_kind, al, beta, nxm, dzq, dyq, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart, xend, istart=>xstart_cpu, iend=>xend_cpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_CPU
#include "CorrectPressure_common.F90"
#undef ON_CPU

end subroutine

#ifdef USE_CUDA
subroutine CorrectPressure_GPU(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
  use param,only: fp_kind, al, beta, nxm, dzq, dyq, nx, lvlhalo
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, istart=>xstart_gpu, iend=>xend_gpu
#else
  use decomp_2d, only: xstart, xend, istart=>xstart, iend=>xend
#endif
  implicit none

#define ON_GPU
#include "CorrectPressure_common.F90"
#undef ON_GPU

end subroutine
#endif


