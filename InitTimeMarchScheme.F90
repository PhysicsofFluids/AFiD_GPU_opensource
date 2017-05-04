!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitTimeMarchScheme.F90                        !
!    CONTAINS: subroutine InitTimeMarchScheme             !
!                                                         ! 
!    PURPOSE: Initialize the time-marching constants for  !
!     the integrator                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitTimeMarchScheme
      use param
      implicit none
      integer ns

      if(nsst.gt.1) then   
        gam(1)=real(8.,fp_kind)/real(15.,fp_kind)
        gam(2)=real(5.,fp_kind)/real(12.,fp_kind)
        gam(3)=real(3.,fp_kind)/real(4.,fp_kind)
        rom(1)=real(0.,fp_kind)
        rom(2)=-real(17.,fp_kind)/real(60.,fp_kind)
        rom(3)=-real(5.,fp_kind)/real(12.,fp_kind)
!m======================================================
      if(ismaster) then
        write(6,100) (gam(ns),ns=1,nsst),(rom(ns),ns=1,nsst)
  100   format(/,5x,'The time scheme is a III order Runge-Kutta' &
        ,4x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
      endif
!m======================================================
      else                                                              
        gam(1)=real(1.5,fp_kind)
        gam(2)=real(0.,fp_kind)
        gam(3)=real(0.,fp_kind)
        rom(1)=-real(0.5,fp_kind)
        rom(2)=real(0.,fp_kind)
        rom(3)=real(0.,fp_kind)
!m======================================================
      if(ismaster) then
        write(6,110) gam(1),rom(1)
  110   format(/,5x,'The time scheme is the Adams-Bashfort',4x, &
         'gam= ',f8.3,4x,'ro= ',f8.3)
      endif
     
!m======================================================                                 
      endif                                                             

      do ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
      end do

      return
      end

