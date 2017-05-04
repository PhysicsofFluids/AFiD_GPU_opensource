!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetTempBCs.F90                                 !
!    CONTAINS: subroutine SetTempBCs                      !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetTempBCs
      use param
      implicit none
      integer :: ic,jc

      do ic=1,nz
       do jc=1,ny
        temptp(jc,ic)=real(0.,fp_kind)
        tempbp(jc,ic)=real(1.,fp_kind)
       enddo
      enddo

#ifdef USE_CUDA
        allocate(temptp_d, source=temptp)
        allocate(tempbp_d, source=tempbp)
#endif
      return
      end
!
