!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!                                                         ! 
!    PURPOSE: Initialization routines. Compute the metric !
!     terms and modified wavenumbers for the pressure     !
!     correction                                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitPressureSolver
      use param
      use fftw_params
      use decomp_2d_fft
      implicit none
      integer :: nymh,nymp,j,i,nzmh,nzmp
      integer  :: kc,km,kp
      real(fp_kind) :: ugmmm,a33icc,a33icp
    
!    Initialize wave number definitions

      nymh=nym/2+1
      nymp=nymh+1

      nzmh=nzm/2+1
      nzmp=nzmh+1

      do i=1,nzmh
        ao(i)=(i-1)*real(2.0,fp_kind)*pi
      enddo
      do i=nzmp,nzm
        ao(i)=-(nzm-i+1)*real(2.0,fp_kind)*pi
      enddo
      do i=1,nzm
        ak1(i)=real(2.0,fp_kind)*(real(1.0,fp_kind)-cos(ao(i)/nzm))*(real(nzm,fp_kind)/zlen)**2
      enddo

      do j=1,nymh
        ap(j)=(j-1)*real(2.0,fp_kind)*pi
      enddo
      do j=nymp,nym
        ap(j)=-(nym-j+1)*real(2.0,fp_kind)*pi
      enddo
      do j=1,nym
        ak2(j)=real(2.0,fp_kind)*(real(1.0,fp_kind)-cos(ap(j)/nym))*(real(nym,fp_kind)/ylen)**2
      enddo

!RO   Initialize Tridiagonal matrices for Poisson solver

      do kc=1,nxm
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dxq/g3rc(kc)
        a33icp=kpc(kc)*dxq/g3rc(kp)
        ugmmm=real(1.0,fp_kind)/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo

#ifdef USE_CUDA
       allocate(ak1_d,source=ak1)
       allocate(ak2_d,source=ak2)
       allocate(ao_d,source=ao)
       allocate(ap_d,source=ap)
       allocate(amphk_d,source=amphk)
       allocate(acphk_d,source=acphk)
       allocate(apphk_d,source=apphk)

#endif

!RO   Initialize Pencil transposes for pressure solver

      call decomp_2d_fft_init

      if (.not.planned) then
        call CreateFFTTmpArrays
        call InitPressureSolverPlans
        planned=.true.
      endif

      return
      end
      
