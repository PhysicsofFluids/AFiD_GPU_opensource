!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteFlowField.F90                             !
!    CONTAINS: subroutine WriteFlowField                  !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine WriteFlowField
      use param
#ifdef USE_CUDA
      use local_arrays, only: vz_d,vy_d,vx_d,temp_d
#endif
      use local_arrays, only: vz,vy,vx,temp
#ifdef USE_HYBRID
      use decomp_2d, only: xstart, xsize_gpu
#endif
      implicit none
      character*30 :: filnam1,dsetname

#if defined(USE_CUDA) && defined(USE_HYBRID)
!JR Copy data from GPU subdomain to CPU before writing to disk
     istat = cudaMemcpy(vx(1, xstart(2)-lvlhalo, xstart(3)), vx_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
     istat = cudaMemcpy(vy(1, xstart(2)-lvlhalo, xstart(3)), vy_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
     istat = cudaMemcpy(vz(1, xstart(2)-lvlhalo, xstart(3)), vz_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
     istat = cudaMemcpy(temp(1, xstart(2)-lvlhalo, xstart(3)), temp_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2) * (xsize_gpu(3)))
#elif defined(USE_CUDA)
      ! Copy solution back to CPU memory
      vx=vx_d
      vy=vy_d
      vz=vz_d
      temp=temp_d
#endif

      filnam1 = trim('continua_temp.h5')
      call HdfWriteRealHalo3D(filnam1,temp)
      filnam1 = trim('continua_vx.h5')
      call HdfWriteRealHalo3D(filnam1,vx)
      filnam1 = trim('continua_vy.h5')
      call HdfWriteRealHalo3D(filnam1,vy)
      filnam1 = trim('continua_vz.h5')
      call HdfWriteRealHalo3D(filnam1,vz)

      if (ismaster) then !EP only write once
       filnam1 = trim('continua_master.h5')
       call HdfCreateBlankFile(filnam1)
 
       dsetname = trim('nx')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nx)
       dsetname = trim('ny')
       call HdfSerialWriteIntScalar(dsetname,filnam1,ny)
       dsetname = trim('nz')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nz)
       dsetname = trim('ylen')
       call HdfSerialWriteRealScalar(dsetname,filnam1,ylen)
       dsetname = trim('zlen')
       call HdfSerialWriteRealScalar(dsetname,filnam1,zlen)
       dsetname = trim('time')
       call HdfSerialWriteRealScalar(dsetname,filnam1,time)
       dsetname = trim('istr3')
       call HdfSerialWriteIntScalar(dsetname,filnam1,istr3)
       dsetname = trim('str3')
       call HdfSerialWriteRealScalar(dsetname,filnam1,str3)

      endif

      end subroutine WriteFlowField
