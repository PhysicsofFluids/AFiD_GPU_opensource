!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInitialConditions
  use param
  use local_arrays, only: vy,vx,temp,vz
#ifdef USE_CUDA
  use cudafor
  use local_arrays, only: vy_d,vx_d,temp_d,vz_d
#endif
  use decomp_2d, only: xstart,xend
  use mpih
  implicit none
  integer :: j,k,i
  real(fp_kind) :: xxx,yyy,eps
  
  eps=real(0.01,fp_kind)
  do i=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do k=1,nxm
           vz(k,j,i)=real(0.0,fp_kind)
           yyy=xm(k) 
           xxx=yc(j)            
           vy(k,j,i)=(real(2.0,fp_kind)*yyy-real(6.0,fp_kind)*yyy**2+real(4.0,fp_kind)*yyy**3) &
                &                  *sin(3*xxx)*eps
           yyy=xc(k)          
           xxx=ym(j)
           vx(k,j,i)=-yyy**2*(real(1.0,fp_kind)-yyy)**2*cos(real(3.1,fp_kind)*xxx)*eps
           
        enddo
     enddo
  enddo
  
  do i=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do k=2,nxm
           temp(k,j,i)= tempbp(j,i) - (tempbp(j,i) - temptp(j,i)) &
                *xc(k)
        enddo
        temp(1 ,j,i)=tempbp(j,i)
        temp(nx,j,i)=temptp(j,i)
     end do
  end do
  
  return                                                            
end subroutine CreateInitialConditions
