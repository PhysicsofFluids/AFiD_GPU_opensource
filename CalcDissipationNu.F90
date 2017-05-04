!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcDissipationNu.F90                          !
!    CONTAINS: subroutine CalcDissipationNu               !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number through the    !
!     global balance equations relating dissipation and   !
!     heat transport.                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_CUDA

#DEFINE BLKX (64)
#DEFINE BLKY (1)

module diss_cuda
contains
      attributes(global) subroutine diss_kernel(vx,vy,vz,tp,xc,xm,disste,dissth, &
                                                gnu_mu,gnu_th,nx,ny,nz,ldx,ldy,dx,dy,dz)
      use precision
      implicit none
      !args
      integer, value, intent(IN) :: nx,ny,nz,ldx,ldy
      real(fp_kind), value, intent(IN) :: dx,dy,dz
      real(fp_kind), device, dimension(1:ldx,1:ldy,nz), intent(IN) :: vx,vy,vz,tp
      real(fp_kind), device, dimension(*), intent(IN) :: xc,xm
      real(fp_kind), device, dimension(1:nx,1:ny) :: disste, dissth
      real(fp_kind), device, dimension(*) :: gnu_mu, gnu_th
      !shared vars
      real(fp_kind), shared :: snu_mu(BLKX*BLKY), snu_th(BLKX*BLKY)
      !local vars
      integer :: i,j,js,k,k1,k2,hid,tidx,tidy
      real(fp_kind) :: dxc,dxm,lec,lem
      real(fp_kind) :: hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz
      real(fp_kind) :: tx,ty,tz
      real(fp_kind) :: nu_th,nu_mu,volt
      real(fp_kind) :: dissipte,dissipte2,dissipth,sumte,sumth
      real(fp_kind) :: vx1,vx2,vx3,vx4,vx5,vx6,vy1,vy2,vy3,vy4,vy5,vy6
      real(fp_kind) :: vz1,vz2,vz3,vz4,vz5,vz6,tp1,tp2,tp3,tp4,tp5,tp6

      nu_mu = real(0.0,fp_kind); nu_th = real(0.0,fp_kind)
      sumte = real(0.0,fp_kind); sumth = real(0.0,fp_kind)

      tidx = threadIdx%x; tidy = threadIdx%y
      k  = tidx + (blockIdx%x-1)*BLKX
      js = tidy + (blockIdx%y-1)*BLKY
      if( k<=nx .and. js<=ny ) then
           k1 = k+1; k2 = k
           dxc=real(1.0,fp_kind)/(xc(k1)-xc(k2))
           lem=(xc(k1)-xc(k2))*dx
           if(k==nx) then
             k1 = k; k2 = k-1
           endif
           dxm=real(1.0,fp_kind)/(xm(k1)-xm(k2))
           lec=(xm(k1)-xm(k2))*dx

           do j = js, ny, gridDim%y*BLKY
 
             tp4=tp(k  ,j,1); tp5=tp(k+1,j,1); tp6=tp(k,j+1,1);
             do i = 1, nz
               tp1=tp4;         tp2=tp5;           tp3=tp6;
               tp4=tp(k,j,i+1); tp5=tp(k+1,j,i+1); tp6=tp(k,j+1,i+1);
               tx = (tp2 - tp1)*dxc
               ty = (tp3 - tp1)*dy
               tz = (tp4 - tp1)*dz
               dissipth = tx*tx + ty*ty + tz*tz
               nu_th = nu_th + dissipth*lem
               sumth = sumth + dissipth
             end do

             vx4=vx(k,j,1); vx5=vx(k+1,j,1); vx6=vx(k,j+1,1);
             do i = 1, nz
               vx1=vx4;         vx2=vx5;           vx3=vx6;
               vx4=vx(k,j,i+1); vx5=vx(k+1,j,i+1); vx6=vx(k,j+1,i+1);
               hxx = (vx2 - vx1)*dxc
               hxy = (vx3 - vx1)*dy
               hxz = (vx4 - vx1)*dz
               dissipte  = (hxx*hxx+hxy*hxy+hxz*hxz)
               nu_mu = nu_mu + (dissipte*lem)
               sumte = sumte + (dissipte)
             end do

             vy4=vy(k,j,1); vy5=vy(k+1,j,1); vy6=vy(k,j+1,1);
             do i = 1, nz
               vy1=vy4;         vy2=vy5;           vy3=vy6;
               vy4=vy(k,j,i+1); vy5=vy(k+1,j,i+1); vy6=vy(k,j+1,i+1);
               hyx = (vy2 - vy1)*dxm
               hyy = (vy3 - vy1)*dy
               hyz = (vy4 - vy1)*dz
               dissipte2 = (hyx*hyx+hyy*hyy+hyz*hyz)
               nu_mu = nu_mu + (dissipte2*lec)
               sumte = sumte + (dissipte2)
             end do

             vz4=vz(k,j,1); vz5=vz(k+1,j,1); vz6=vz(k,j+1,1);
             do i = 1, nz
               vz1=vz4;         vz2=vz5;           vz3=vz6;
               vz4=vz(k,j,i+1); vz5=vz(k+1,j,i+1); vz6=vz(k,j+1,i+1);
               hzx = (vz2 - vz1)*dxm
               hzy = (vz3 - vz1)*dy
               hzz = (vz4 - vz1)*dz
               dissipte2 = (hzx*hzx+hzy*hzy+hzz*hzz)
               nu_mu = nu_mu + (dissipte2*lec)
               sumte = sumte + (dissipte2)
             end do

           end do

           disste(k,js) = disste(k,js) + sumte
           dissth(k,js) = dissth(k,js) + sumth

        end if
        tidx = threadIdx%x+(threadIdx%y-1)*BLKX
        hid = (BLKX*BLKY)/2
        snu_mu(tidx) = nu_mu
        snu_th(tidx) = nu_th
        call syncthreads()
        do while(hid >= 1)
          if(tidx<=hid) then
            snu_mu(tidx) = snu_mu(tidx) + snu_mu(tidx+hid)
            snu_th(tidx) = snu_th(tidx) + snu_th(tidx+hid)
          endif
          call syncthreads()
          hid = hid/2
        enddo
        if(tidx==1) then
          gnu_mu(blockIdx%x + (blockIdx%y-1)*gridDim%x) = snu_mu(1)
          gnu_th(blockIdx%x + (blockIdx%y-1)*gridDim%x) = snu_th(1)
        endif

      end subroutine diss_kernel
end module diss_cuda
#endif

      subroutine CalcDissipationNu
      use mpih
#if defined(USE_CUDA)
      use cudafor
      use diss_cuda
      use param, only: fp_kind,dx,dy,dz,xm_d,xc_d,nxm,nym,nzm,pra,time,ismaster
      use local_arrays,only: vz_d,vy_d,vx_d,temp_d
      use stat_arrays,only: disste_d,dissth_d, stat_columns
      use decomp_2d, only: gnu_mu=>work1_r_d, gnu_th=>work2_r_d, hnu_mu=>work1_r, hnu_th=>work2_r, xstart_gpu, xend_gpu
#endif

#if defined(USE_HYBRID) || !defined(USE_CUDA)
      use param, only: fp_kind,dx,dy,dz,xm,xc,nxm,nym,nzm,pra,time,ismaster
      use local_arrays,only: vz,vy,vx,temp
      use stat_arrays,only: disste,dissth
      use decomp_2D, only: xstart_cpu, xend_cpu
      use nvtx
#endif

      implicit none
      integer :: i,j,k
      real(fp_kind) :: nu_th,nu_mu,volt
#ifdef USE_CUDA
      TYPE(dim3) :: tBlock, grid
#endif

#if defined(USE_HYBRID) || !defined(USE_CUDA)
      integer :: ip,jp,kp
      real(fp_kind) :: dxm(nxm),dxc(nxm),lem(nxm),lec(nxm)
      real(fp_kind) :: hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz
      real(fp_kind) :: tx,ty,tz
      real(fp_kind) :: dissipte,dissipte2,dissipth
#endif
      
      nu_mu = real(0.0,fp_kind)
      nu_th = real(0.0,fp_kind)

#if defined(USE_CUDA) 
     tBlock = dim3(BLKX,BLKY,1)
     grid   = dim3(ceiling(real(nxm)/BLKX),ceiling(real(stat_columns)/BLKY),1)

     call diss_kernel<<<grid,tBlock>>>(vx_d(1,xstart_gpu(2),xstart_gpu(3)),vy_d(1,xstart_gpu(2),xstart_gpu(3)),vz_d(1,xstart_gpu(2),xstart_gpu(3)),&
       temp_d(1,xstart_gpu(2),xstart_gpu(3)), xc_d,xm_d,disste_d,dissth_d,gnu_mu,gnu_th,nxm,(xend_gpu(2)-xstart_gpu(2)+1),(xend_gpu(3)-xstart_gpu(3)+1),&
       SIZE(vx_d,1),SIZE(vx_d,2), dx, dy, dz )

#endif

#if defined(USE_HYBRID) || !defined(USE_CUDA)
      call nvtxStartRangeAsync("dis_cpu")
! RS: Calculated geometrical quantities outside the main loop.
! RS: Could potentially be moved to initialization routine 
      do k=1,nxm
        kp=k+1
        dxm(k)=1.0/(xm(kp)-xm(k))
        dxc(k)=1.0/(xc(kp)-xc(k))
        lem(k)=(xc(kp)-xc(k))*dx
        lec(k)=(xm(kp)-xm(k))*dx
      enddo
      lec(nxm)=(xm(nxm)-xm(nxm-1))*dx
      dxm(nxm)=1.0/(xm(nxm)-xm(nxm-1))

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart_cpu,xend_cpu,vz,vy,vx,temp) &
!$OMP  SHARED(nxm) &
!$OMP  SHARED(dxc,dxm,dy,dz,lec,lem) &
!$OMP  PRIVATE(i,j,k,ip,jp,kp) &
!$OMP  PRIVATE(dissipte,dissipte2,dissipth) &
!$OMP  PRIVATE(hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz) &
!$OMP  PRIVATE(tx,ty,tz) &
!$OMP  REDUCTION(+:nu_mu, nu_th, disste, dissth)
      do i=xstart_cpu(3),xend_cpu(3)
       ip= i+1
        do j=xstart_cpu(2),xend_cpu(2)
        jp=j+1

        do k=1,nxm
        kp=k+1

!       Viscous dissipation rate
!                       1  |         | 2
!                     ---- | nabla  u|
!                      Re  |         |
       hxx=(vx(kp,j,i)-vx(k,j,i))*dxc(k)
       hxy=(vx(k,jp,i)-vx(k,j,i))*dy
       hxz=(vx(k,j,ip)-vx(k,j,i))*dz

       hyx=(vy(kp,j,i)-vy(k,j,i))*dxm(k)
       hyy=(vy(k,jp,i)-vy(k,j,i))*dy
       hyz=(vy(k,j,ip)-vy(k,j,i))*dz

       hzx=(vz(kp,j,i)-vz(k,j,i))*dxm(k)
       hzy=(vz(k,jp,i)-vz(k,j,i))*dy
       hzz=(vz(k,j,ip)-vz(k,j,i))*dz

       dissipte  = (hxx*hxx+hxy*hxy+hxz*hxz)
       dissipte2 = (hyx*hyx+hyy*hyy+hyz*hyz)+(hzx*hzx+hzy*hzy+hzz*hzz)

       nu_mu = nu_mu + dissipte*lem(k)+dissipte2*lec(k)

!      Thermal gradient dissipation rate
!                       1  |         | 2
!                     ---- | nabla  T|
!                      Pe  |         |

       tx=(temp(kp,j,i)-temp(k,j,i))*dxc(k)
       ty=(temp(k,jp,i)-temp(k,j,i))*dy
       tz=(temp(k,j,ip)-temp(k,j,i))*dz

       dissipth  = tx*tx + ty*ty + tz*tz
       nu_th = nu_th+dissipth*lem(k)

!JR Swapping critical section with a reduction
!!$OMP CRITICAL
       disste(k,1) =  disste(k,1) + (dissipte + dissipte2) 
       dissth(k,1) =  dissth(k,1) + dissipth               
!!$OMP END CRITICAL

       end do
       end do
       end do
!$OMP  END PARALLEL DO

       call nvtxEndRangeAsync

#endif

#ifdef USE_CUDA
!JR  Add GPU contribution
     hnu_mu(1:grid%x*grid%y) = gnu_mu(1:grid%x*grid%y)
     hnu_th(1:grid%x*grid%y) = gnu_th(1:grid%x*grid%y)

     do i=1,grid%x*grid%y
       nu_mu = nu_mu + hnu_mu(i)
       nu_th = nu_th + hnu_th(i)
     end do
#endif


       nu_mu = nu_mu*pra

       call MpiSumRealScalar(nu_th)
       call MpiSumRealScalar(nu_mu)
      
       volt = real(1.0,fp_kind)/(real(nxm,fp_kind)*real(nzm,fp_kind)*real(nym,fp_kind))

      if(ismaster) then
      nu_mu = nu_mu*volt + 1
      nu_th = nu_th*volt 
      open(92,file='nu_diss.out',status='unknown',access='sequential', &
       position='append')
      write(92,*) time,nu_mu,nu_th
      close(92)
      endif

      return   
      end
