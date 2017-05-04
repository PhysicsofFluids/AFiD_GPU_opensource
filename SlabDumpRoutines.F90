!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SlabDumper.F90                                 !
!    CONTAINS: subroutines SlabDumper, InitializeSlabDump !
!     DumpSingleSlab                                      ! 
!                                                         ! 
!    PURPOSE: Auxiliary routines used for memory allocs   !
!     and memory freeing                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SlabDumper
      use param
#ifdef USE_CUDA
      use local_arrays, only: vx_d,vy_d,vz_d,temp_d
#endif
      use local_arrays, only: temp,vx,vy,vz
      use stat3_param
      use decomp_2d, only: xstart,xend
#ifdef USE_HYBRID
      use decomp_2d, only: xsize_gpu
#endif
      implicit none
      integer :: i,j,m
      real(fp_kind),dimension(xstart(2):xend(2),xstart(3):xend(3)) :: &
     &      vxcc,vycc,vzcc,tempcc
      character*70 :: filnam
      character*1 :: charm


#if defined(USE_CUDA) && defined(USE_HYBRID)
!JR Copy data from GPU subdomain to CPU before writing to disk
     istat = cudaMemcpy(vx(1, xstart(2)-lvlhalo, xstart(3)), vx_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
     istat = cudaMemcpy(vy(1, xstart(2)-lvlhalo, xstart(3)), vy_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
     istat = cudaMemcpy(vz(1, xstart(2)-lvlhalo, xstart(3)), vz_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
     istat = cudaMemcpy(temp(1, xstart(2)-lvlhalo, xstart(3)), temp_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
#elif defined(USE_CUDA)
!MF Copy data from GPU to CPU before writing to disk
     vx   = vx_d
     vy   = vy_d
     vz   = vz_d
     temp = temp_d
#endif


!EP   Slabs
!EP   cell center only vx

      do m=1,9
!$OMP  PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(i,j)
        do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           vxcc(j,i) = (vx(kslab(m),j,i)+vx(kslab(m)+1,j,i))*real(0.5,fp_kind)
           vycc(j,i) = vy(kslab(m),j,i)
           vzcc(j,i) = vz(kslab(m),j,i)
           tempcc(j,i) = temp(kslab(m),j,i)
          enddo
         enddo
!$OMP  END PARALLEL DO
      write(charm,28) m
   28 format(i1.1)
      filnam='slab'//charm//'vx_'
      call DumpSingleSlab(vxcc,filnam)
      filnam='slab'//charm//'vy_'
      call DumpSingleSlab(vycc,filnam)
      filnam='slab'//charm//'vz_'
      call DumpSingleSlab(vzcc,filnam)
      filnam='slab'//charm//'temp_'
      call DumpSingleSlab(tempcc,filnam)
      enddo

      return
      end subroutine SlabDumper

!===========================================================================

      subroutine InitializeSlabDump
      use param
      use stat3_param
      implicit none
      integer :: i,k,j
      real(fp_kind) :: xmloc
      character(len=4) :: dummy

!EP   Read from stst3.in
      
      open(unit=19,file='stst3.in',status='old')
        read(19,301) dummy
        read(19,*) (xslab(i),i=2,9)
301     format(a4)                
      close(19)

!EP   Compute which kslab corresponds to which xslab
      
      kslab = 2
      
        do k=2,nxm
          xmloc=xm(k)
          do j=2,9
            if(xm(k).gt.xslab(j).and.xm(k-1).lt.xslab(j)) then
             kslab(j) = k
            endif
          enddo
        enddo


!EP   Write probe and slab locations
      
      if (ismaster) then
      open(unit=23,file='stst3locs.out',status='unknown')
        rewind(23)
        write(23,*) (kslab(i),i=1,9)
      close(23)
      endif

      return
      end subroutine InitializeSlabDump

!==================================================================
      
      subroutine DumpSingleSlab(var,filnam)
      USE param
      use mpih
      USE hdf5
      use decomp_2d, only: xstart,xend
      IMPLICIT none

      real(fp_kind), intent(in) :: var(xstart(2):xend(2) &
     &                  ,xstart(3):xend(3))

      real(fp_kind) :: tprfi
      integer :: itime

      character*70,intent(in) :: filnam
      character*70 :: namfile,dsetname
      character*8 :: ipfi

      tprfi = real(1.,fp_kind)/tout
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i8.8)

      namfile=trim('./stst3/'//trim(filnam)//trim(ipfi)//'.h5')
      dsetname = trim('var')

      call HdfWriteReal2D(dsetname,namfile,var)

      if(ismaster) then
       dsetname = trim('time')
       call HdfSerialWriteRealScalar(dsetname,namfile,time)
      endif

      return                                                          
      end subroutine DumpSingleSlab
