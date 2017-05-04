!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         ! 
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateGrid
  use param
  use AuxiliaryRoutines
  implicit none
  
  real(fp_kind) :: x1,x2,x3
  real(fp_kind) :: a33, a33m, a33p
  real(fp_kind) :: delet, etain, tstr3
  real(fp_kind) :: z2dp
  
  real(fp_kind), allocatable, dimension(:) :: etaz, etazm
  
  integer :: i, j, kc, km, kp
  integer :: nxmo, nclip
  
  
  do kc=1,nxm
     kmv(kc)=kc-1
     kpv(kc)=kc+1
     if(kc.eq.1) kmv(kc)=kc
     if(kc.eq.nxm) kpv(kc)=kc
  end do
  
  do kc=1,nxm
     kpc(kc)=kpv(kc)-kc
     kmc(kc)=kc-kmv(kc)
  end do

#ifdef USE_CUDA
  allocate(kmv_d, source=kmv)
  allocate(kpv_d, source=kpv)
  allocate(kmc_d, source=kmc)
  allocate(kpc_d, source=kpc)
#endif  
  
  !
  !     UNIFORM (HORIZONTAL DIRECTIONS) GRID
  !
  do  i=1,nz
     x1=real(i-1,fp_kind)/real(nzm,fp_kind)
     zc(i)= zlen*x1
  end do
  
  do i=1,nzm
     zm(i)=(zc(i)+zc(i+1))*real(0.5,fp_kind)
  end do
  
  do j=1,ny
     x2=real(j-1,fp_kind)/real(nym,fp_kind)
     yc(j)= ylen*x2
  end do
  
  do j=1,nym
     ym(j)=(yc(j)+yc(j+1))*real(0.5,fp_kind)
  end do
  
  !
  !     VERTICAL COORDINATE DEFINITION
  !
  !     OPTION 0: UNIFORM CLUSTERING
  !
  call AllocateReal1DArray(etaz,1,nx+500)
  call AllocateReal1DArray(etazm,1,nx+500)
  
  if (istr3.eq.0) then
     do kc=1,nx
        x3=real(kc-1,fp_kind)/real(nxm,fp_kind)
        etaz(kc)=alx3*x3
        xc(kc)=etaz(kc)
     enddo
  endif

  !
  !     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
  !
  
  tstr3=tanh(str3)
  
  if (istr3.eq.4) then
     xc(1)=real(0.0,fp_kind)
     do kc=2,nx
        z2dp=real(2*kc-nx-1,fp_kind)/real(nxm,fp_kind)
        xc(kc)=(real(1.,fp_kind)+tanh(str3*z2dp)/tstr3)*real(0.5,fp_kind)*alx3
        if( xc(kc).lt.real(0.,fp_kind) .or. xc(kc).gt.alx3 )then
           write(*,*)'Grid is too streched: ','zc(',kc,')=',xc(kc)
           stop
        endif
     end do
  end if
  
  !
  !     OPTION 6: CLIPPED CHEBYCHEV-TYPE CLUSTERING
  !
  
  
  if(istr3.eq.6) then
     nclip = int(str3)
     nxmo = nx+nclip+nclip
     do kc=1,nxmo
        etazm(kc)=+cos(pi*(real(kc,fp_kind)-real(0.5,fp_kind))/real(nxmo,fp_kind))
     end do
     do kc=1,nx
        etaz(kc)=etazm(kc+nclip)
     end do
     delet = etaz(1)-etaz(nx)
     etain = etaz(1)
     do kc=1,nx
        etaz(kc)=etaz(kc)/(real(0.5,fp_kind)*delet)
     end do
     xc(1) = real(0.,fp_kind)
     do kc=2,nxm
        xc(kc) = alx3*(real(1.,fp_kind)-etaz(kc))*real(0.5,fp_kind)
     end do
     xc(nx) = alx3
  endif

  call DestroyReal1DArray(etaz)
  call DestroyReal1DArray(etazm)
  
  !m-----------------------------------------
  !
  !     METRIC FOR UNIFORM DIRECTIONS
  !
  
  dx=real(nxm,fp_kind)/alx3
  dy=real(nym,fp_kind)/ylen
  dz=real(nzm,fp_kind)/zlen
  
  dxq=dx*dx                                                      
  dyq=dy*dy                                                      
  dzq=dz*dz                                                      
#ifdef USE_CUDA
  dx_d=dx
  dy_d=dy
  dz_d=dz
  dxq_d=dxq
  dyq_d=dyq
  dzq_d=dzq
#endif
  
  !
  !     STAGGERED COORDINATES AND
  !     METRIC QUANTITIES FOR NON-UNIFORM 
  !     DIRECTIONS
  !
  
  do kc=1,nxm
     xm(kc)=(xc(kc)+xc(kc+1))*real(0.5,fp_kind)
     g3rm(kc)=(xc(kc+1)-xc(kc))*dx
  enddo
  do kc=2,nxm
     g3rc(kc)=(xc(kc+1)-xc(kc-1))*dx*real(0.5,fp_kind)
  enddo
  g3rc(1)=(xc(2)-xc(1))*dx
  g3rc(nx)= (xc(nx)-xc(nxm))*dx

#ifdef USE_CUDA
   allocate(xc_d,source=xc)
   allocate(xm_d,source=xm)
  allocate(g3rm_d, source=g3rm)
  allocate(g3rc_d, source=g3rc)
#endif

  !
  !     WRITE GRID INFORMATION
  !
  do kc=1,nxm
     udx3m(kc) = dx/g3rm(kc)
     udx3c(kc) = dx/g3rc(kc)
  end do
  udx3c(nx) = dx/g3rc(nx)

#ifdef USE_CUDA
  allocate(udx3m_d, source=udx3m)
  allocate(udx3c_d, source=udx3c)  
#endif

  !m====================================================
  if(ismaster) then
     open(unit=78,file='axicor.out',status='unknown')
     do kc=1,nx
        write(78,345) kc,xc(kc),xm(kc),g3rc(kc),g3rm(kc)
     end do
     close(78)
345  format(i4,4(2x,e23.15))
     !m===================================================
     !
     !     QUANTITIES FOR DERIVATIVES
     !
     open(unit=78,file='fact3.out',status='unknown')
     do kc=1,nxm
        write(78,*) kc,udx3m(kc),udx3c(kc)
     end do
     write(78,*) nx,udx3m(nxm),udx3c(nx)
     close(78)
  endif
  
  !
  !    COEFFICIENTS FOR DIFFERENTIATION FOR NON-UNIFORM GRID
  !
  !    Q3 DIFFERENTIATION (CENTERED VARIABLE)
  !
  
  am3ck(1)=real(0.,fp_kind)
  ap3ck(1)=real(0.,fp_kind)
  ac3ck(1)=real(1.,fp_kind)
  am3ck(nx)=real(0.,fp_kind)
  ap3ck(nx)=real(0.,fp_kind)
  ac3ck(nx)=real(1.,fp_kind)
  
  do kc=2,nxm
     km=kc-1
     kp=kc+1
     a33=dxq/g3rc(kc)
     a33p=real(1.0,fp_kind)/g3rm(kc)
     a33m=real(1.0,fp_kind)/g3rm(km)
     ap3ck(kc)=a33*a33p
     am3ck(kc)=a33*a33m
     ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
  enddo

#ifdef USE_CUDA
  allocate(am3ck_d, source=am3ck)
  allocate(ap3ck_d, source=ap3ck)
  allocate(ac3ck_d, source=ac3ck)
#endif

  !
  !    Q1/Q2 DIFFERENTIATION (STAGGERED VARIABLE)
  !
  !
  
  do kc=2,nxm-1
     kp=kc+1
     km=kc-1
     a33=dxq/g3rm(kc)
     a33p= +a33/g3rc(kp)
     a33m= +a33/g3rc(kc)
     ap3sk(kc)=a33p
     am3sk(kc)=a33m
     ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
  enddo
  !    
  !    LOWER WALL BOUNDARY CONDITIONS (INSLWS SETS NO-SLIP vs STRESS-FREE WALL)
  !    
  kc=1
  kp=kc+1
  a33=dxq/g3rm(kc)
  a33p= +a33/g3rc(kp)
  a33m= +a33/g3rc(kc)
  ap3sk(kc)=a33p
  am3sk(kc)=real(0.,fp_kind)
  ac3sk(kc)=-(a33p+inslws*a33m*real(2.,fp_kind))
  
  !    
  !    UPPER WALL BOUNDARY CONDITIONS (INSLWN SETS NO-SLIP vs STRESS-FREE WALL)
  !    
  
  kc=nxm
  kp=kc+1
  a33=dxq/g3rm(kc)
  a33p= +a33/g3rc(kp)
  a33m= +a33/g3rc(kc)
  am3sk(kc)=a33m
  ap3sk(kc)=real(0.,fp_kind)
  ac3sk(kc)=-(a33m+inslwn*a33p*real(2.,fp_kind))
  
  am3ssk(1)=real(0.,fp_kind)
  ap3ssk(1)=real(0.,fp_kind)
  ac3ssk(1)=real(1.,fp_kind)
  
  !
  !    TEMPERATURE DIFFERENTIATION (CENTERED VARIABLE)
  !
  
  
  do kc=2,nxm
     kp=kc+1
     km=kc-1
     a33=dxq/g3rc(kc)
     a33p=real(1.,fp_kind)/g3rm(kc)
     a33m=real(1.,fp_kind)/g3rm(km)
     ap3ssk(kc)=a33*a33p
     am3ssk(kc)=a33*a33m
     ac3ssk(kc)=-(ap3ssk(kc)+am3ssk(kc))
  enddo
  
  am3ssk(nx)=real(0.,fp_kind)
  ap3ssk(nx)=real(0.,fp_kind)
  ac3ssk(nx)=real(1.,fp_kind)

#ifdef USE_CUDA
  allocate(am3sk_d, source=am3sk)
  allocate(ap3sk_d, source=ap3sk)
  allocate(ac3sk_d, source=ac3sk)
  allocate(am3ssk_d, source=am3ssk)
  allocate(ap3ssk_d, source=ap3ssk)
  allocate(ac3ssk_d, source=ac3ssk)
#endif
  
  return                                                            
end subroutine CreateGrid

