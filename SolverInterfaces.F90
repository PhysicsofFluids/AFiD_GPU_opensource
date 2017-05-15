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
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:),intent(IN) :: vx,vy,vz,temp
      !pgi$ ignore_tkr qcap
      real(fp_kind),intent(OUT)  :: qcap
      real(fp_kind), dimension(:), intent(IN) :: udx3c
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsVX_GPU(vx,vy,vz,temp,qcap,udx3c)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:),intent(IN) :: vx,vy,vz,temp
      !pgi$ ignore_tkr qcap
      real(fp_kind), device,intent(OUT)  :: qcap
      real(fp_kind), device, dimension(:), intent(IN) :: udx3c
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ExplicitTermsVY
    subroutine ExplicitTermsVY_CPU(vx,vy,vz,dph,udx3m,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:),intent(IN) :: vx,vy,vz
      real(fp_kind), dimension(:,:,:), intent(OUT)  :: dph
      real(fp_kind), dimension(:), intent(IN) :: udx3m
      integer, dimension(:), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsVY_GPU(vx,vy,vz,dph,udx3m,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:),intent(IN) :: vx,vy,vz
      real(fp_kind), device, dimension(:,:,:), intent(OUT)  :: dph
      real(fp_kind), device, dimension(:), intent(IN) :: udx3m
      integer, device, dimension(:), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ExplicitTermsVZ
    subroutine ExplicitTermsVZ_CPU(vx,vy,vz,dq,udx3m,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:),intent(IN) :: vx,vy,vz
      !pgi$ ignore_tkr dq
      real(fp_kind), intent(OUT)  :: dq
      real(fp_kind), dimension(:), intent(IN) :: udx3m
      integer, dimension(:), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsVZ_GPU(vx,vy,vz,dq,udx3m,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:),intent(IN) :: vx,vy,vz
      !pgi$ ignore_tkr dq
      real(fp_kind), device, intent(OUT)  :: dq
      real(fp_kind), device, dimension(:), intent(IN) :: udx3m
      integer, device, dimension(:), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ExplicitTermsTemp
    subroutine ExplicitTermsTemp_CPU(vx,vy,vz,temp,hro,udx3c)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:),intent(IN) :: vx,vy,vz,temp
      real(fp_kind), dimension(:,:,:),intent(OUT)  :: hro
      real(fp_kind), dimension(:), intent(IN) :: udx3c
    end subroutine

#ifdef USE_CUDA
    subroutine ExplicitTermsTemp_GPU(vx,vy,vz,temp,hro,udx3c)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:),intent(IN) :: vx,vy,vz,temp
      real(fp_kind), device, dimension(:,:,:),intent(OUT)  :: hro
      real(fp_kind), device, dimension(:), intent(IN) :: udx3c
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateVX
    subroutine ImplicitAndUpdateVX_CPU(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: vx, pr, rux
      !pgi$ ignore_tkr rhs, qcap
      real(fp_kind) :: rhs, qcap
      real(fp_kind), dimension(:), intent(IN) :: am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateVX_GPU(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:) :: vx, pr, rux
      !pgi$ ignore_tkr rhs, qcap
      real(fp_kind), device :: rhs, qcap
      real(fp_kind), device, dimension(:), intent(IN) :: am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateVY
    subroutine ImplicitAndUpdateVY_CPU(vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: vy, pr, ruy
      !pgi$ ignore_tkr rhs
      real(fp_kind) :: rhs
      real(fp_kind), dimension(:,:,:), intent(IN)  :: dph
      real(fp_kind), dimension(:), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, dimension(:), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateVY_GPU(vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:) :: vy, pr, ruy
      !pgi$ ignore_tkr rhs
      real(fp_kind), device :: rhs
      real(fp_kind), device, dimension(:,:,:), intent(IN)  :: dph
      real(fp_kind), device, dimension(:), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, device, dimension(:), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateVZ
    subroutine ImplicitAndUpdateVZ_CPU(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only : fp_kind
      real(fp_kind), dimension(:,:,:) :: vz, pr, ruz
      !pgi$ ignore_tkr rhs, dq
      real(fp_kind) :: rhs, dq
      real(fp_kind), dimension(:), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, dimension(:), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateVZ_GPU(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
      use param, only : fp_kind
      real(fp_kind), device, dimension(:,:,:) :: vz, pr, ruz
      !pgi$ ignore_tkr rhs, dq
      real(fp_kind), device :: rhs, dq
      real(fp_kind), device, dimension(:), intent(IN) :: am3sk,ac3sk,ap3sk
      integer, device, dimension(:), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface ImplicitAndUpdateTemp
    subroutine ImplicitAndUpdateTemp_CPU(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
      use param, only : fp_kind
      real(fp_kind), dimension(:,:,:) :: temp, rutemp, hro
      !pgi$ ignore_tkr rhs
      real(fp_kind) :: rhs
      real(fp_kind), dimension(:,:) :: tempbp, temptp
      real(fp_kind), dimension(:), intent(IN) :: am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine ImplicitAndUpdateTemp_GPU(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
      use param, only : fp_kind
      real(fp_kind), device, dimension(:,:,:) :: temp, rutemp, hro
      !pgi$ ignore_tkr rhs
      real(fp_kind), device :: rhs
      real(fp_kind), device, dimension(:,:) :: tempbp, temptp
      real(fp_kind), device, dimension(:), intent(IN) :: am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface SolveImpEqnUpdate_X
    subroutine SolveImpEqnUpdate_X_CPU(vx,rhs,am3ssk,ac3ck,ap3ssk)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: vx, rhs
      real(fp_kind), dimension(:), intent(IN) :: am3ssk,ac3ck,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine SolveImpEqnUpdate_X_GPU(vx,rhs,am3ssk,ac3ck,ap3ssk)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:) :: vx, rhs
      real(fp_kind), device, dimension(:), intent(IN) :: am3ssk,ac3ck,ap3ssk
    end subroutine

#endif
  end interface

!----------------------------------------------------------
  interface SolveImpEqnUpdate_YZ
    subroutine SolveImpEqnUpdate_YZ_CPU(q,rhs,am3sk,ac3sk,ap3sk)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: q, rhs
      real(fp_kind), dimension(:), intent(IN) :: am3sk,ac3sk,ap3sk
    end subroutine

#ifdef USE_CUDA
    subroutine SolveImpEqnUpdate_YZ_GPU(q,rhs,am3sk,ac3sk,ap3sk)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:) :: q, rhs
      real(fp_kind), device, dimension(:), intent(IN) :: am3sk,ac3sk,ap3sk
    end subroutine

#endif
  end interface

!----------------------------------------------------------
  interface SolveImpEqnUpdate_Temp
    subroutine SolveImpEqnUpdate_Temp_CPU(temp,rhs,am3ssk,ac3ssk,ap3ssk)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: temp, rhs
      real(fp_kind), dimension(:), intent(IN) :: am3ssk,ac3ssk,ap3ssk
    end subroutine

#ifdef USE_CUDA
    subroutine SolveImpEqnUpdate_Temp_GPU(temp,rhs,am3ssk,ac3ssk,ap3ssk)
      use param, only: fp_kind
      implicit none
      real(fp_kind), device, dimension(:,:,:) :: temp, rhs
      real(fp_kind), device, dimension(:), intent(IN) :: am3ssk,ac3ssk,ap3ssk
    end subroutine

#endif
  end interface

!----------------------------------------------------------
  interface CalcLocalDivergence
    subroutine CalcLocalDivergence_CPU(vx,vy,vz,dph,udx3m)
      use param, only:  fp_kind
      real(fp_kind), dimension(:,:,:),intent(IN) :: vx,vy,vz
      real(fp_kind), dimension(:,:,:)  :: dph
      real(fp_kind), dimension(:), intent(IN) :: udx3m
    end subroutine

#ifdef USE_CUDA
    subroutine CalcLocalDivergence_GPU(vx,vy,vz,dph,udx3m)
      use param, only:  fp_kind
      real(fp_kind), device, dimension(:,:,:),intent(IN) :: vx,vy,vz
      real(fp_kind), device, dimension(:,:,:)  :: dph
      real(fp_kind), device, dimension(:), intent(IN) :: udx3m
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface CorrectVelocity
    subroutine CorrectVelocity_CPU(vx,vy,vz,dphhalo,udx3c,kmv)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: vx,vy,vz
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), intent(IN) :: dphhalo
      real(fp_kind), dimension(:), intent(IN) :: udx3c
      integer, dimension(:), intent(IN) :: kmv
    end subroutine

#ifdef USE_CUDA
    subroutine CorrectVelocity_GPU(vx,vy,vz,dphhalo,udx3c,kmv)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:) :: vx,vy,vz
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), device, intent(IN) :: dphhalo
      real(fp_kind), device, dimension(:), intent(IN) :: udx3c
      integer, device, dimension(:), intent(IN) :: kmv
    end subroutine
#endif
  end interface

!----------------------------------------------------------
  interface CorrectPressure
    subroutine CorrectPressure_CPU(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), dimension(:,:,:) :: pr
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), intent(IN) :: dphhalo
      real(fp_kind), dimension(:), intent(IN) :: amphk,acphk,apphk
      integer, dimension(:), intent(IN) :: kmv,kpv
    end subroutine

#ifdef USE_CUDA
    subroutine CorrectPressure_GPU(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
      use param, only: fp_kind
      real(fp_kind), device, dimension(:,:,:) :: pr
      !pgi$ ignore_tkr dphhalo
      real(fp_kind), device, intent(IN) :: dphhalo
      real(fp_kind), device, dimension(:), intent(IN) :: amphk,acphk,apphk
      integer, device, dimension(:), intent(IN) :: kmv,kpv
    end subroutine
#endif
  end interface


end module
