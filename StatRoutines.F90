!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_CUDA
module stat_cuda
contains
      attributes(global) subroutine stats_kernel(vx,vy,vz,tp,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,t1,t2,t3,t4,tx,nx,ny,nz,ldx,ldy,factor)
      use precision
      implicit none
      !args
      integer, value, intent(IN) :: nx,ny,nz,ldx,ldy
      real(fp_kind), value, intent(IN) :: factor
      real(fp_kind), device, dimension(1:ldx,1:ldy*nz), intent(IN) :: vx,vy,vz,tp
      real(fp_kind), device, dimension(1:nx,1:ny) :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,t1,t2,t3,t4,tx
      !local vars
      integer :: i,ij,ijp,j,k,tidy,tidx
      real(fp_kind) :: vx1,vy1,vz1,tp1
      real(fp_kind) :: vx2,vy2,vz2,tp2
      real(fp_kind) :: sx1,sy1,sz1,st1,stx
      real(fp_kind) :: sx2,sy2,sz2,st2
      real(fp_kind) :: sx3,sy3,sz3,st3
      real(fp_kind) :: sx4,sy4,sz4,st4
      k = threadIdx%x + (blockIdx%x-1)*blockDim%x
      ij = threadIdx%y + (blockIdx%y-1)*blockDim%y
      if( k<=nx .and. ij<=ny*nz ) then
           sx1=0.0; sy1=0.0; sz1=0.0; st1=0.0; stx=0.0
           sx2=0.0; sy2=0.0; sz2=0.0; st2=0.0
           sx3=0.0; sy3=0.0; sz3=0.0; st3=0.0
           sx4=0.0; sy4=0.0; sz4=0.0; st4=0.0
           do while(ij <= ny*nz)
             j = mod(ij-1,ny) + 1
             i = (ij+ny-1)/ny
             ijp = (i-1)*ldy + j

             vx1 = vx(k,ijp)
             vy1 = vy(k,ijp)
             vz1 = vz(k,ijp)
             tp1 = tp(k,ijp)

             vx2 = vx1*vx1
             vy2 = vy1*vy1
             vz2 = vz1*vz1
             tp2 = tp1*tp1

             sx1 = sx1 + vx1*factor
             sy1 = sy1 + vy1*factor
             sz1 = sz1 + vz1*factor
             st1 = st1 + tp1*factor
             stx = stx + tp1*vx1*factor

             sx2 = sx2 + vx2*factor
             sy2 = sy2 + vy2*factor
             sz2 = sz2 + vz2*factor
             st2 = st2 + tp2*factor

             sx3 = sx3 + vx1*vx2*factor
             sy3 = sy3 + vy1*vy2*factor
             sz3 = sz3 + vz1*vz2*factor
             st3 = st3 + tp1*tp2*factor

             sx4 = sx4 + vx2*vx2*factor
             sy4 = sy4 + vy2*vy2*factor
             sz4 = sz4 + vz2*vz2*factor
             st4 = st4 + tp2*tp2*factor

             ij = ij + gridDim%y * blockDim%y
           end do
           j = threadIdx%y + (blockIdx%y-1)*blockDim%y
           x1(k,j) = x1(k,j) + sx1
           y1(k,j) = y1(k,j) + sy1
           z1(k,j) = z1(k,j) + sz1
           t1(k,j) = t1(k,j) + st1
           tx(k,j) = tx(k,j) + stx
           x2(k,j) = x2(k,j) + sx2
           y2(k,j) = y2(k,j) + sy2
           z2(k,j) = z2(k,j) + sz2
           t2(k,j) = t2(k,j) + st2
           x3(k,j) = x3(k,j) + sx3
           y3(k,j) = y3(k,j) + sy3
           z3(k,j) = z3(k,j) + sz3
           t3(k,j) = t3(k,j) + st3
           x4(k,j) = x4(k,j) + sx4
           y4(k,j) = y4(k,j) + sy4
           z4(k,j) = z4(k,j) + sz4
           t4(k,j) = t4(k,j) + st4
        end if
      end subroutine stats_kernel
end module stat_cuda
#endif

      subroutine CalcStats
      use param
#if defined(USE_CUDA)
      use cudafor
      use stat_cuda
      use local_arrays, only: vz_d,vy_d,vx_d,temp_d
      use stat_arrays,  only: vx_m1_d,vy_m1_d,vz_m1_d,tp_m1_d &
                             ,vx_m2_d,vy_m2_d,vz_m2_d,tp_m2_d &
                             ,vx_m3_d,vy_m3_d,vz_m3_d,tp_m3_d &
                             ,vx_m4_d,vy_m4_d,vz_m4_d,tp_m4_d &
                             ,tpvx_m1_d, nstatsamples, tstat, stat_columns
      use decomp_2d, only: xstart_gpu,xend_gpu
#endif
#if defined(USE_HYBRID) || !defined(USE_CUDA)
      use local_arrays, only: vz,vy,vx,temp
      use stat_arrays
      use nvtx
      use decomp_2d, only: xstart_cpu,xend_cpu
#endif

      use mpih
      implicit none
      real(fp_kind) :: usnzm,usnym,factor
      real(fp_kind) :: vx1,vy1,vz1,tp1
      real(fp_kind) :: vx2,vy2,vz2,tp2
      integer :: i,j,k
#if defined(USE_CUDA)
      type(dim3) :: tBlock, grid
#endif

      nstatsamples = nstatsamples + 1
      tstat        = tstat        + dt

      usnym = real(1.0,fp_kind)/nym
      usnzm = real(1.0,fp_kind)/nzm
      factor= usnym*usnzm*dt

#ifdef USE_CUDA
     tBlock = dim3(64,1,1)
     grid   = dim3(ceiling(real(nxm)/tBlock%x),stat_columns,1)
     call stats_kernel<<<grid,tBlock>>>(vx_d(1,xstart_gpu(2),xstart_gpu(3)),vy_d(1,xstart_gpu(2),xstart_gpu(3)),vz_d(1,xstart_gpu(2),xstart_gpu(3)),&
                                   temp_d(1,xstart_gpu(2),xstart_gpu(3)),vx_m1_d,vx_m2_d,vx_m3_d,vx_m4_d,vy_m1_d,vy_m2_d,vy_m3_d,vy_m4_d,vz_m1_d, &
                                   vz_m2_d,vz_m3_d,vz_m4_d,tp_m1_d,tp_m2_d,tp_m3_d,tp_m4_d,tpvx_m1_d,(nxm),(xend_gpu(2)-xstart_gpu(2)+1), &
                                   (xend_gpu(3)-xstart_gpu(3)+1),SIZE(vx_d,1), SIZE(vx_d,2), factor )
#endif

#if defined(USE_HYBRID) || !defined(USE_CUDA)
      call nvtxStartRangeAsync("Stats_cpu")
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) &
      !$OMP PRIVATE(i,j,k) &
      !$OMP PRIVATE(vx1,vx2,vy1,vy2,vz1,vz2,tp1,tp2) &
      !$OMP SHARED(vx, vy, vz, temp, nxm, factor,xstart_cpu,xend_cpu)&
      !$OMP REDUCTION(+:vx_m1, vx_m2, vx_m3, vx_m4) &
      !$OMP REDUCTION(+:vy_m1, vy_m2, vy_m3, vy_m4) &
      !$OMP REDUCTION(+:vz_m1, vz_m2, vz_m3, vz_m4) &
      !$OMP REDUCTION(+:tp_m1, tp_m2, tp_m3, tp_m4) &
      !$OMP REDUCTION(+:tpvx_m1)

      do i=xstart_cpu(3),xend_cpu(3)
       do j=xstart_cpu(2),xend_cpu(2)
         do k=1,nxm
           vx1 = vx(k,j,i)
           vx2 = vx1*vx1
           vx_m1(k,1)   = vx_m1(k,1)   +   vx1*factor
           vx_m2(k,1)   = vx_m2(k,1)   +   vx2*factor
           vx_m3(k,1)   = vx_m3(k,1)   +   vx1*vx2*factor
           vx_m4(k,1)   = vx_m4(k,1)   +   vx2*vx2*factor

           vy1 = vy(k,j,i)
           vy2 = vy1*vy1
           vy_m1(k,1)   = vy_m1(k,1)   +   vy1*factor
           vy_m2(k,1)   = vy_m2(k,1)   +   vy2*factor
           vy_m3(k,1)   = vy_m3(k,1)   +   vy1*vy2*factor
           vy_m4(k,1)   = vy_m4(k,1)   +   vy2*vy2*factor

           vz1 = vz(k,j,i)
           vz2 = vz1*vz1
           vz_m1(k,1)   = vz_m1(k,1)   +   vz1*factor
           vz_m2(k,1)   = vz_m2(k,1)   +   vz2*factor
           vz_m3(k,1)   = vz_m3(k,1)   +   vz1*vz2*factor
           vz_m4(k,1)   = vz_m4(k,1)   +   vz2*vz2*factor

           tp1 = temp(k,j,i)
           tp2 = tp1*tp1
           tp_m1(k,1) = tp_m1(k,1) +   tp1*factor
           tp_m2(k,1) = tp_m2(k,1) +   tp2*factor
           tp_m3(k,1) = tp_m3(k,1) +   tp1*tp2*factor
           tp_m4(k,1) = tp_m4(k,1) +   tp2*tp2*factor

           tpvx_m1(k,1)= tpvx_m1(k,1)+ tp1*vx1*factor
 
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO
      call nvtxEndRangeAsync
#endif

      return  
      end
!    
!***********************************************************************
      subroutine WriteStats
      use mpih
      use param
      use stat_arrays
      use hdf5
      use decomp_2d, only: xstart,xend

      implicit none

      integer :: nstatsamples_old, i,j,k
      real    :: tstat_old

      character*30 dsetname_vxm1
      character*30 dsetname_vym1
      character*30 dsetname_vzm1

      character*30 dsetname_vxm2
      character*30 dsetname_vym2
      character*30 dsetname_vzm2

      character*30 dsetname_vxm3
      character*30 dsetname_vym3
      character*30 dsetname_vzm3

      character*30 dsetname_vxm4
      character*30 dsetname_vym4
      character*30 dsetname_vzm4

      character*30 dsetname_tempm1
      character*30 dsetname_tempm2
      character*30 dsetname_tempm3
      character*30 dsetname_tempm4

      character*30 dsetname_tempvxm1
      
      character*30 dsetname_dissth
      character*30 dsetname_disste
      character*30 filnam,dsetname,dsetname2
      logical :: fexist

      filnam = trim('stafield_master.h5')

      dsetname_vxm1 = trim('vx_m1')
      dsetname_vym1 = trim('vy_m1')
      dsetname_vzm1 = trim('vz_m1')

      dsetname_vxm2 = trim('vx_m2')
      dsetname_vym2 = trim('vy_m2')
      dsetname_vzm2 = trim('vz_m2')

      dsetname_vxm3 = trim('vx_m3')
      dsetname_vym3 = trim('vy_m3')
      dsetname_vzm3 = trim('vz_m3')

      dsetname_vxm4 = trim('vx_m4')
      dsetname_vym4 = trim('vy_m4')
      dsetname_vzm4 = trim('vz_m4')

      dsetname_tempm1 = trim('temp_m1')
      dsetname_tempm2 = trim('temp_m2')
      dsetname_tempm3 = trim('temp_m3')
      dsetname_tempm4 = trim('temp_m4')

      dsetname_tempvxm1 = trim('tempvx_m1')

      dsetname_dissth = trim('dissth')
      dsetname_disste = trim('disste')

      dsetname = trim('averaging_time')
      dsetname2 = trim('averaging_time2')

      inquire(file=filnam,exist=fexist)
      if (.not.fexist) then 
        if(ismaster) write(6,*) 'Unable to read statistical files'
        if(ismaster) write(6,*) 'Restarting statistics from zero' 
        readstats=.false.
      end if
       

      if (ismaster) then
       if(readstats) then
        call HdfSerialReadIntScalar(dsetname,filnam,nstatsamples_old)
        call HdfSerialReadRealScalar(dsetname2,filnam,tstat_old)
        nstatsamples = nstatsamples + nstatsamples_old
        tstat       = tstat       + tstat_old
       else 
        call HdfCreateBlankFile(filnam)
       endif
      end if

#if defined(USE_CUDA)
#ifdef USE_HYBRID
! Add first column of GPU data to CPU data
  vx_m1(:,2) = vx_m1_d(:,1)
  vy_m1(:,2) = vy_m1_d(:,1)
  vz_m1(:,2) = vz_m1_d(:,1)
  tp_m1(:,2) = tp_m1_d(:,1)

  vx_m2(:,2) = vx_m2_d(:,1)
  vy_m2(:,2) = vy_m2_d(:,1)
  vz_m2(:,2) = vz_m2_d(:,1)
  tp_m2(:,2) = tp_m2_d(:,1)

  vx_m3(:,2) = vx_m3_d(:,1)
  vy_m3(:,2) = vy_m3_d(:,1)
  vz_m3(:,2) = vz_m3_d(:,1)
  tp_m3(:,2) = tp_m3_d(:,1)

  vx_m4(:,2) = vx_m4_d(:,1)
  vy_m4(:,2) = vy_m4_d(:,1)
  vz_m4(:,2) = vz_m4_d(:,1)
  tp_m4(:,2) = tp_m4_d(:,1)

  tpvx_m1(:,2) = tpvx_m1_d(:,1)

  do k = 1, nxm
    vx_m1(k,1) = vx_m1(k,1) + vx_m1(k,2)
    vy_m1(k,1) = vy_m1(k,1) + vy_m1(k,2)
    vz_m1(k,1) = vz_m1(k,1) + vz_m1(k,2)
    tp_m1(k,1) = tp_m1(k,1) + tp_m1(k,2)

    vx_m2(k,1) = vx_m2(k,1) + vx_m2(k,2)
    vy_m2(k,1) = vy_m2(k,1) + vy_m2(k,2)
    vz_m2(k,1) = vz_m2(k,1) + vz_m2(k,2)
    tp_m2(k,1) = tp_m2(k,1) + tp_m2(k,2)

    vx_m4(k,1) = vx_m3(k,1) + vx_m3(k,2)
    vy_m4(k,1) = vy_m3(k,1) + vy_m3(k,2)
    vz_m4(k,1) = vz_m3(k,1) + vz_m3(k,2)
    tp_m4(k,1) = tp_m3(k,1) + tp_m3(k,2)
    
    vx_m4(k,1) = vx_m1(k,1) + vx_m1(k,2)
    vy_m4(k,1) = vy_m1(k,1) + vy_m1(k,2)
    vz_m4(k,1) = vz_m1(k,1) + vz_m1(k,2)
    tp_m4(k,1) = tp_m1(k,1) + tp_m1(k,2)

    tpvx_m1(k,1) = tpvx_m1(k,1) + tpvx_m1(k,2)
  end do

  ! Copy remaining GPU columns to CPU
  vx_m1(:,2:stat_columns) = vx_m1_d(:,2:stat_columns)
  vy_m1(:,2:stat_columns) = vy_m1_d(:,2:stat_columns)
  vz_m1(:,2:stat_columns) = vz_m1_d(:,2:stat_columns)
  tp_m1(:,2:stat_columns) = tp_m1_d(:,2:stat_columns)

  vx_m2(:,2:stat_columns) = vx_m2_d(:,2:stat_columns)
  vy_m2(:,2:stat_columns) = vy_m2_d(:,2:stat_columns)
  vz_m2(:,2:stat_columns) = vz_m2_d(:,2:stat_columns)
  tp_m2(:,2:stat_columns) = tp_m2_d(:,2:stat_columns)

  vx_m3(:,2:stat_columns) = vx_m3_d(:,2:stat_columns)
  vy_m3(:,2:stat_columns) = vy_m3_d(:,2:stat_columns)
  vz_m3(:,2:stat_columns) = vz_m3_d(:,2:stat_columns)
  tp_m3(:,2:stat_columns) = tp_m3_d(:,2:stat_columns)

  vx_m4(:,2:stat_columns) = vx_m4_d(:,2:stat_columns)
  vy_m4(:,2:stat_columns) = vy_m4_d(:,2:stat_columns)
  vz_m4(:,2:stat_columns) = vz_m4_d(:,2:stat_columns)
  tp_m4(:,2:stat_columns) = tp_m4_d(:,2:stat_columns)

  tpvx_m1(:,2:stat_columns) = tpvx_m1_d(:,2:stat_columns)

#else
!   Copy the arrays back to CPU
  vx_m1 = vx_m1_d
  vy_m1 = vy_m1_d
  vz_m1 = vz_m1_d
  tp_m1 = tp_m1_d

  vx_m2 = vx_m2_d
  vy_m2 = vy_m2_d
  vz_m2 = vz_m2_d
  tp_m2 = tp_m2_d

  vx_m3 = vx_m3_d
  vy_m3 = vy_m3_d
  vz_m3 = vz_m3_d
  tp_m3 = tp_m3_d

  vx_m4 = vx_m4_d
  vy_m4 = vy_m4_d
  vz_m4 = vz_m4_d
  tp_m4 = tp_m4_d

  tpvx_m1 = tpvx_m1_d
#endif


!SUM COLUMNS INTO FIRST COLUMN
        do j=2,stat_columns
          do k=1,nxm
               vx_m1(k,1) = vx_m1(k,1) + vx_m1(k,j)
               vy_m1(k,1) = vy_m1(k,1) + vy_m1(k,j)
               vz_m1(k,1) = vz_m1(k,1) + vz_m1(k,j)
               tp_m1(k,1) = tp_m1(k,1) + tp_m1(k,j)

               vx_m2(k,1) = vx_m2(k,1) + vx_m2(k,j)
               vy_m2(k,1) = vy_m2(k,1) + vy_m2(k,j)
               vz_m2(k,1) = vz_m2(k,1) + vz_m2(k,j)
               tp_m2(k,1) = tp_m2(k,1) + tp_m2(k,j)

               vx_m3(k,1) = vx_m3(k,1) + vx_m3(k,j)
               vy_m3(k,1) = vy_m3(k,1) + vy_m3(k,j)
               vz_m3(k,1) = vz_m3(k,1) + vz_m3(k,j)
               tp_m3(k,1) = tp_m3(k,1) + tp_m3(k,j)

               vx_m4(k,1) = vx_m4(k,1) + vx_m4(k,j)
               vy_m4(k,1) = vy_m4(k,1) + vy_m4(k,j)
               vz_m4(k,1) = vz_m4(k,1) + vz_m4(k,j)
               tp_m4(k,1) = tp_m4(k,1) + tp_m4(k,j)

               tpvx_m1(k,1) = tpvx_m1(k,1) + tpvx_m1(k,j)
            end do
         end do

      if(disscal) then
#ifdef USE_HYBRID
        disste(:,2) = disste_d(:,1)
        dissth(:,2) = dissth_d(:,1)

        do k=1,nxm
          disste(k,1) = disste(k,1) + disste(k,2)
          dissth(k,1) = dissth(k,1) + dissth(k,2)
        end do

        disste(:,2:stat_columns) = disste_d(:,2:stat_columns)
        dissth(:,2:stat_columns) = dissth_d(:,2:stat_columns)

#else
        disste = disste_d
        dissth = dissth_d
#endif

        do j=2,stat_columns
          do k=1,nxm
               disste(k,1) = disste(k,1) + disste(k,j)
               dissth(k,1) = dissth(k,1) + dissth(k,j)
            end do
         end do
       endif
#endif
      call StatReadReduceWrite(vx_m1,filnam,dsetname_vxm1)
      call StatReadReduceWrite(vy_m1,filnam,dsetname_vym1)
      call StatReadReduceWrite(vz_m1,filnam,dsetname_vzm1)

      call StatReadReduceWrite(vx_m2,filnam,dsetname_vxm2)
      call StatReadReduceWrite(vy_m2,filnam,dsetname_vym2)
      call StatReadReduceWrite(vz_m2,filnam,dsetname_vzm2)

      call StatReadReduceWrite(vx_m3,filnam,dsetname_vxm3)
      call StatReadReduceWrite(vy_m3,filnam,dsetname_vym3)
      call StatReadReduceWrite(vz_m3,filnam,dsetname_vzm3)

      call StatReadReduceWrite(vx_m4,filnam,dsetname_vxm4)
      call StatReadReduceWrite(vy_m4,filnam,dsetname_vym4)
      call StatReadReduceWrite(vz_m4,filnam,dsetname_vzm4)

      call StatReadReduceWrite(tp_m1,filnam,dsetname_tempm1)
      call StatReadReduceWrite(tp_m2,filnam,dsetname_tempm2)
      call StatReadReduceWrite(tp_m3,filnam,dsetname_tempm3)
      call StatReadReduceWrite(tp_m4,filnam,dsetname_tempm4)

      call StatReadReduceWrite(tpvx_m1,filnam,dsetname_tempvxm1)

      if(disscal) then

       do k=1,nxm
               disste(k,1) = disste(k,1) / (ren*real(nym)*real(nzm))
               dissth(k,1) = dissth(k,1) / (pec*real(nym)*real(nzm))
       end do
 
       call StatReadReduceWrite(dissth,filnam,dsetname_dissth)
       call StatReadReduceWrite(disste,filnam,dsetname_disste)
      end if

      if (ismaster) then

       write(*,*) 'tstat',tstat
       call HdfSerialWriteIntScalar(dsetname,filnam,nstatsamples)
       call HdfSerialWriteRealScalar(dsetname2,filnam,tstat)

       dsetname = trim('X_cordin')
       call HdfSerialWriteReal1D(dsetname,filnam,xm,nxm)

       dsetname = trim('Rayleigh Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,ray)

       dsetname = trim('Prandtl Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,pra)


      endif

      return  
      end
