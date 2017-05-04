!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcPlateNu.F90                                !
!    CONTAINS: subroutine CalcPlateNu                     !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number at the top     !
!     and bottom plates and output to a file.             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcPlateNu
  use param
#if defined(USE_CUDA) && defined(USE_HYBRID)
  use local_arrays, only: temp_h => temp, temp => temp_d, nusslw, nussup
#elif defined(USE_CUDA)
  use local_arrays, only: temp => temp_d, nusslw, nussup
#else
  use local_arrays, only: temp, nusslw, nussup
#endif
  use mpih
#ifdef USE_HYBRID
  use decomp_2d, only: xstart => xstart_gpu, xend => xend_gpu, xstart_cpu, xend_cpu
#else
  use decomp_2d, only: xstart,xend
#endif
  implicit none
  integer :: j,i
  real(fp_kind) ::  nuslow, nusupp
  real(fp_Kind) :: del,deln
  

  nuslow = real(0.,fp_kind)
  nusupp = real(0.,fp_kind)
  del  = real(1.0,fp_kind)/(xc(2)-xc(1))
  deln = real(1.0,fp_kind)/(xc(nx)-xc(nxm))

#ifdef USE_CUDA
  !$cuf kernel do(2) <<<*,*>>>
#else
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(xstart,xend,temp,del,deln) &
  !$OMP   SHARED(nxm,nx) &
  !$OMP   PRIVATE(i,j) &
  !$OMP   REDUCTION(+:nuslow) &
  !$OMP   REDUCTION(+:nusupp)
#endif
  do i=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        nuslow = nuslow + (temp(1,j,i)-temp(2,j,i))*del
        nusupp = nusupp + (temp(nxm,j,i)-temp(nx,j,i))*deln
     enddo
  end do
#ifndef USE_CUDA
  !$OMP END PARALLEL DO
#endif

#ifdef USE_HYBRID
  !$OMP  PARALLEL DO &
  !$OMP   DEFAULT(none) &
  !$OMP   SHARED(xstart_cpu,xend_cpu,temp_h,del,deln) &
  !$OMP   SHARED(nxm,nx) &
  !$OMP   PRIVATE(i,j) &
  !$OMP   REDUCTION(+:nuslow) &
  !$OMP   REDUCTION(+:nusupp)
  do i=xstart_cpu(3),xend_cpu(3)
     do j=xstart_cpu(2),xend_cpu(2)
        nuslow = nuslow + (temp_h(1,j,i)-temp_h(2,j,i))*del
        nusupp = nusupp + (temp_h(nxm,j,i)-temp_h(nx,j,i))*deln
     enddo
  end do
  !$OMP END PARALLEL DO

#endif

  nuslow = nuslow / real(nzm*nym,fp_kind)
  nusupp = nusupp / real(nzm*nym,fp_kind)


  
  call MpiSumRealScalar(nuslow)
  call MpiSumRealScalar(nusupp)
  nusslw = nuslow
  nussup = nusupp
  if(ismaster) then
     open(97,file="nu_plate.out",status='unknown', &
          access='sequential',position='append')
     write(97,546) time, nuslow, nusupp
546  format(4(1x,e14.6))
     close(97)
  endif
  
  return         
end subroutine CalcPlateNu
