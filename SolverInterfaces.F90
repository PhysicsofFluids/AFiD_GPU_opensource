!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolverInterfaces.F90                           !
!    CONTAINS: interface ExplicitTermsVX                  !
!              interface ExplicitTermsVY                  !
!              interface ExplicitTermsVZ                  !
!              interface ImplicitAndUpdateVX              !
!              interface ImplicitAndUpdateVY              !
!              interface ImplicitAndUpdateVZ              !
!              interface ImplicitAndUpdateTemp            !
!              interface SolveImpEqnUpdate_X              !
!              interface SolveImpEqnUpdate_YZ             !
!              interface SolveImpEqnUpdate_Temp           !
!              interface CalcLocalDivergence              !
!              interface CorrectVelocity                  !
!              interface CorrectPressure                  !
!                                                         ! 
!    PURPOSE: Declare generic interfaces for solver       !
!             routines for CPU and GPU operation          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module solver_interfaces
!----------------------------------------------------------
  interface ExplicitTermsVX
    subroutine ExplicitTermsVX_CPU(vx,vy,vz,temp,qcap,udx3c)
      use param, only: fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz,temp
      !pgi$ ignore_tkr qcap
      real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: qcap
      real(fp_kind), dimension(1:nx), intent(IN) :: udx3c
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsVX_GPU(vx,vy,vz,temp,qcap,udx3c)
      use param, only: fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz,temp
      !pgi$ ignore_tkr qcap
      real(fp_kind), device, dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: qcap
      real(fp_kind), device, dimension(1:nx), intent(IN) :: udx3c
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ExplicitTermsVY
    subroutine ExplicitTermsVY_CPU(vx,vy,vz,dph,udx3m,kmv,kpv)
      use param, only: fp_kind, nxm, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
      real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(out)  :: dph
      real(fp_kind), dimension(1:nx), intent(IN) :: udx3m
      integer, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsVY_GPU(vx,vy,vz,dph,udx3m,kmv,kpv)
      use param, only: fp_kind, nxm, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
      real(fp_kind), device, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(out)  :: dph
      real(fp_kind), device, dimension(1:nx), intent(IN) :: udx3m
      integer, device, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ExplicitTermsVZ
    subroutine ExplicitTermsVZ_CPU(vx,vy,vz,dq,udx3m,kmv,kpv)
      use param, only: fp_kind, nxm, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
      !pgi$ ignore_tkr dq
      real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(out)  :: dq
      real(fp_kind), dimension(1:nx), intent(IN) :: udx3m
      integer, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsVZ_GPU(vx,vy,vz,dq,udx3m,kmv,kpv)
      use param, only: fp_kind, nxm, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
      !pgi$ ignore_tkr dq
      real(fp_kind), device, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(out)  :: dq
      real(fp_kind), device, dimension(1:nx), intent(IN) :: udx3m
      integer, device, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ExplicitTermsTemp
    subroutine ExplicitTermsTemp_CPU(vx,vy,vz,temp,hro,udx3c)
      use param, only: fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz,temp
      real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: hro
      real(fp_kind), dimension(1:nx), intent(IN) :: udx3c
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsTemp_GPU(vx,vy,vz,temp,hro,udx3c)
      use param, only: fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz,temp
      real(fp_kind), device, dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)),intent(OUT)  :: hro
      real(fp_kind), device, dimension(1:nx), intent(IN) :: udx3c
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateVX
    subroutine ImplicitAndUpdateVX_CPU(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
      use param, only: fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx,pr
      !pgi$ ignore_tkr rhs, qcap
      real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3)) :: rhs, rux, qcap
      real(fp_kind), dimension(1:nx), intent(IN) :: am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateVX_GPU(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
      use param, only: fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, mytype
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx,pr
      !pgi$ ignore_tkr rhs, qcap
      real(fp_kind), device, dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3)) :: rhs, rux, qcap
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateVY
    subroutine ImplicitAndUpdateVY_CPU(vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only: fp_kind, nxm, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vy,pr
      !pgi$ ignore_tkr rhs
      real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, ruy
      real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(IN)  :: dph
      real(fp_kind), dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateVY_GPU(vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only: fp_kind, nx, nxm, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, mytype
#else
      use decomp_2d, only: xstart, xend, mytype
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vy,pr
      !pgi$ ignore_tkr rhs
      real(fp_kind), device, dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, ruy
      real(fp_kind), device, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(IN)  :: dph
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, device, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateVZ
    subroutine ImplicitAndUpdateVZ_CPU(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only : fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vz,pr
      !pgi$ ignore_tkr rhs, dq
      real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, ruz, dq
      real(fp_kind), dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateVZ_GPU(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only : fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, mytype
#else
      use decomp_2d, only: xstart, xend, mytype
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vz,pr
      !pgi$ ignore_tkr rhs, dq
      real(fp_kind), device, dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, ruz, dq
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, device, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateTemp
    subroutine ImplicitAndUpdateTemp_CPU(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
      use param, only : fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: temp
      !pgi$ ignore_tkr rhs
      real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, rutemp, hro
      real(fp_kind), dimension(1:ny,1:nz) :: tempbp, temptp
      real(fp_kind), dimension(1:nx), intent(IN) :: am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateTemp_GPU(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
      use param, only : fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, mytype
#else
      use decomp_2d, only: xstart, xend, mytype
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: temp
      !pgi$ ignore_tkr rhs
      real(fp_kind), device, dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs, rutemp, hro
      real(fp_kind), device, dimension(1:ny,1:nz) :: tempbp, temptp
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface SolveImpEqnUpdate_X
    subroutine SolveImpEqnUpdate_X_CPU(vx,rhs,am3ssk,ac3ck,ap3ssk)
      use param,only: fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx
      real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs
      real(fp_kind), dimension(1:nx), intent(IN) :: am3ssk,ac3ck,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine SolveImpEqnUpdate_X_GPU(vx,rhs,am3ssk,ac3ck,ap3ssk)
      use param,only: fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx
      real(fp_kind), device, dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3ssk,ac3ck,ap3ssk
    end subroutine

#endif
  end interface

!----------------------------------------------------------
  interface SolveImpEqnUpdate_YZ
    subroutine SolveImpEqnUpdate_YZ_CPU(q,rhs,am3sk,ac3sk,ap3sk)
      use param, only: fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: q
      real(fp_kind), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)) :: rhs
      real(fp_kind), dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
    end subroutine

#ifdef USE_CUDA
    subroutine SolveImpEqnUpdate_YZ_GPU(q,rhs,am3sk,ac3sk,ap3sk)
      use param, only: fp_kind, al, beta, nxm, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: q
      real(fp_kind), device, dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)) :: rhs
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3sk,ac3sk,ap3sk
    end subroutine

#endif
  end interface

!----------------------------------------------------------
  interface SolveImpEqnUpdate_Temp
    subroutine SolveImpEqnUpdate_Temp_CPU(temp,rhs,am3ssk,ac3ssk,ap3ssk)
      use param,only: fp_kind, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: temp
      real(fp_kind), dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs
      real(fp_kind), dimension(1:nx), intent(IN) :: am3ssk,ac3ssk,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine SolveImpEqnUpdate_Temp_GPU(temp,rhs,am3ssk,ac3ssk,ap3ssk)
      use param,only: fp_kind, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: temp
      real(fp_kind), device, dimension(1:nx ,xstart(2):xend(2),xstart(3):xend(3))  :: rhs
      real(fp_kind), device, dimension(1:nx), intent(IN) :: am3ssk,ac3ssk,ap3ssk
    end subroutine

#endif
  end interface

!----------------------------------------------------------
  interface CalcLocalDivergence
    subroutine CalcLocalDivergence_CPU(vx,vy,vz,dph,udx3m)
      use param,only:  fp_kind, nxm, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
      real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3))  :: dph
      real(fp_kind), dimension(1:nx), intent(IN) :: udx3m
    end subroutine

#ifdef USE_CUDA
    subroutine CalcLocalDivergence_GPU(vx,vy,vz,dph,udx3m)
      use param,only:  fp_kind, nxm, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu
#else
      use decomp_2d, only: xstart, xend
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: vx,vy,vz
      real(fp_kind), device, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3))  :: dph
      real(fp_kind), device, dimension(1:nx), intent(IN) :: udx3m
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface CorrectVelocity
    subroutine CorrectVelocity_CPU(vx,vy,vz,dphhalo,udx3c,kmv)
      use param,only: fp_kind, nxm, lvlhalo, nx
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx,vy,vz
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: dphhalo
      real(fp_kind), dimension(1:nx), intent(IN) :: udx3c
      integer, dimension(1:nx), intent(IN) :: kmv
    end subroutine

#ifdef USE_CUDA
    subroutine CorrectVelocity_GPU(vx,vy,vz,dphhalo,udx3c,kmv)
      use param,only: fp_kind, nxm, lvlhalo, nx
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, mytype
#else
      use decomp_2d, only: xstart, xend, mytype
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: vx,vy,vz
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), device, dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(IN) :: dphhalo
      real(fp_kind), device, dimension(1:nx), intent(IN) :: udx3c
      integer, device, dimension(1:nx), intent(IN) :: kmv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface CorrectPressure
    subroutine CorrectPressure_CPU(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
      use param,only: fp_kind, nxm, nx, lvlhalo
      use decomp_2d, only: xstart, xend
      implicit none
      real(fp_kind), dimension(1:nx ,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: pr
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo), intent(IN) :: dphhalo
      real(fp_kind), dimension(1:nx), intent(IN) :: amphk,acphk,apphk
      integer, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine CorrectPressure_GPU(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
      use param,only: fp_kind, nxm, nx, lvlhalo
#ifdef USE_HYBRID
      use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, mytype
#else
      use decomp_2d, only: xstart, xend, mytype
#endif
      implicit none
      real(fp_kind), device, dimension(1:nx ,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: pr
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), device, dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo), intent(IN) :: dphhalo
      real(fp_kind), device, dimension(1:nx), intent(IN) :: amphk,acphk,apphk
      integer, device, dimension(1:nx), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface


end module
