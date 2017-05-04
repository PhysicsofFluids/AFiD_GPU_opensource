!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to X pencil

  subroutine transpose_z_to_x_real(src, dst, decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%COL_INFO%SND_P
    call mem_split_zx_real(src, s1, s2, s3, work1, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_zx_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%y1dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%COL_INFO%RCV_P
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y1cnts_s, decomp%y1disp_s, &
            real_type, work2, decomp%x1cnts_s, decomp%x1disp_s, &
            real_type, decomp%COL_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%y1count, &
         real_type, work2_r, decomp%x1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    call mem_merge_zx_real(work2, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_zx_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif
    
    return
  end subroutine transpose_z_to_x_real


#ifdef OCC
  subroutine transpose_z_to_x_real_start(handle, src, dst, sbuf, rbuf)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    
    ! rearrange source array as send buffer
    call mem_split_zx_real(src, s1, s2, s3, sbuf, dims(1), &
         decomp%y1dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y1count, real_type, &
         rbuf, decomp%x1count, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y1cnts, decomp%y1disp, real_type, &
         rbuf, decomp%x1cnts, decomp%x1disp, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_z_to_x_real_start


  subroutine transpose_z_to_x_real_wait(handle, src, dst, sbuf, rbuf)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_zx_real(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

    return
  end subroutine transpose_z_to_x_real_wait
#endif


  subroutine transpose_z_to_x_complex(src, dst, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: ierror

      call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz, &
        dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz,MPI_COMM_WORLD,ierror)

    return

  end subroutine transpose_z_to_x_complex

#ifdef USE_CUDA
  subroutine transpose_z_to_x_real_d(src, dst, decomp)

    implicit none

    real(mytype), device, dimension(:,:,:), intent(IN) :: src
    real(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat
    integer :: i,j,k, m,i1,i2,pos

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
    !call mem_split_zx_complex(src, s1, s2, s3, work1_c, dims(1), &
    !     decomp%y1dist, decomp)

    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y1dist(m)-1
       end if

       pos = decomp%y1disp(m) + 1

       istat = cudaMemcpy2D( work1_r_d(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2,s1*(i2-i1+1), s3 )

    end do

    work1_r = work1_r_d

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
            real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
            real_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    !call mem_merge_zx_complex(work2_c, d1, d2, d3, dst, dims(1), &
    !     decomp%x1dist, decomp)

    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if

       pos = decomp%x1disp(m) + 1

       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_r(pos), i2-i1+1, i2-i1+1,d2*d3 )

    end do


    return
  end subroutine transpose_z_to_x_real_d


  subroutine transpose_z_to_x_complex_d(src, dst, decomp)

    implicit none

    complex(mytype), device, dimension(:,:,:), intent(IN) :: src
    complex(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat
    integer :: i,j,k, m,i1,i2,pos,iter,dest,sorc

call nvtxStartRange("tranZX",2)

#ifdef EPA2A
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    if(IAND(dims(1),dims(1)-1)==0) then

    ! rearrange source array as send buffer
    do iter=0,dims(1)-1
       m = IEOR(col_rank,iter)
       if(iter==0) then ! if exchanging with self
         pos = decomp%w1disp(m) + 1
         istat = cudaMemcpyAsync( work1_c_d(pos), src(1,1,decomp%z1idx(m)), decomp%z1cnts(m), a2a_comp )
       else
         pos = decomp%z1disp(m) + 1
         !istat = cudaMemcpy2DAsync( work1_r_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, 0 )
         istat = cudaMemcpyAsync( work1_c(pos), src(1,1,decomp%z1idx(m)), decomp%z1cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(1)-1
       m = IEOR(col_rank,iter)
      call MPI_IRECV( work2_c(decomp%w1disp(m)+1), decomp%w1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(1)-1
       m = IEOR(col_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_c(decomp%z1disp(m)+1), decomp%z1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, ierror)
call nvtxEndRangeAsync
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%w1disp(m) + 1
       istat = cudaMemcpyAsync( work1_c_d(pos), work2_c(pos), decomp%w1cnts(m), a2a_h2d )
       istat = cudaEventRecord( a2a_event(iter), a2a_h2d )
    end do

    ! rearrange receive buffer
    do iter=0,dims(1)-1
       m = IEOR(col_rank,iter)
       pos = decomp%w1disp(m) + 1
       istat = cudaStreamWaitEvent( a2a_comp, a2a_event(iter), 0 )
       istat = cudaMemcpy2DAsync( dst(decomp%x1idx(m),1,1), d1, work1_c_d(pos), decomp%x1dist(m), decomp%x1dist(m), d2*d3, stream=a2a_comp )
    end do


   else

    ! rearrange source array as send buffer
    do iter=0,dims(1)-1
       dest = mod(col_rank + iter, dims(1))
       m = dest
       !m = IEOR(col_rank,iter)
       if(iter==0) then ! if exchanging with self
         pos = decomp%w1disp(m) + 1
         istat = cudaMemcpyAsync( work1_c_d(pos), src(1,1,decomp%z1idx(m)), decomp%z1cnts(m), a2a_comp )
       else
         pos = decomp%z1disp(m) + 1
         !istat = cudaMemcpy2DAsync( work1_r_d(pos), decomp%x1dist(m),
         !src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, 0 )
         istat = cudaMemcpyAsync( work1_c(pos), src(1,1,decomp%z1idx(m)), decomp%z1cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(1)-1
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = sorc
       !m = IEOR(col_rank,iter)
       call MPI_IRECV( work2_c(decomp%w1disp(m)+1), decomp%w1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(1)-1
       dest = mod(col_rank + iter, dims(1))
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = dest
       !m = IEOR(col_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter) 
       call MPI_SEND( work1_c(decomp%z1disp(m)+1), decomp%z1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, ierror)
call nvtxEndRangeAsync
       m = sorc
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%w1disp(m) + 1
       istat = cudaMemcpyAsync( work1_c_d(pos), work2_c(pos), decomp%w1cnts(m), a2a_h2d )
       istat = cudaEventRecord( a2a_event(iter), a2a_h2d )
    end do

    ! rearrange receive buffer
    do iter=0,dims(1)-1
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = sorc
       !m = IEOR(col_rank,iter)
       pos = decomp%w1disp(m) + 1
       istat = cudaStreamWaitEvent( a2a_comp, a2a_event(iter), 0 )
       istat = cudaMemcpy2DAsync( dst(decomp%x1idx(m),1,1), d1, work1_c_d(pos), decomp%x1dist(m), decomp%x1dist(m), d2*d3, stream=a2a_comp )
    end do


#if 0
    istat = cudaMemcpy( work1_c, src, s1*s2*s3 )

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_c, decomp%z1cnts, decomp%z1disp, &
         complex_type, work2_c, decomp%w1cnts, decomp%w1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if

       pos = decomp%w1disp(m) + 1

       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_c(pos), i2-i1+1, i2-i1+1, d2*d3 )

    end do
#endif

    end if
#else
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
    !call mem_split_zx_complex(src, s1, s2, s3, work1_c, dims(1), &
    !     decomp%y1dist, decomp)

    istat = cudaMemcpy( work1_c, src, s1*s2*s3 )

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_c, decomp%z1cnts, decomp%z1disp, &
         complex_type, work2_c, decomp%w1cnts, decomp%w1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    !call mem_merge_zx_complex(work2_c, d1, d2, d3, dst, dims(1), &
    !     decomp%x1dist, decomp)

    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if

       pos = decomp%w1disp(m) + 1

       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_c(pos), i2-i1+1, i2-i1+1, d2*d3 )

    end do
#endif

call nvtxEndRange
    return
  end subroutine transpose_z_to_x_complex_d
#endif

#ifdef OCC
  subroutine transpose_z_to_x_complex_start(handle, src, dst, sbuf, &
       rbuf)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    
    ! rearrange source array as send buffer
    call mem_split_zx_complex(src, s1, s2, s3, sbuf, dims(1), &
         decomp%y1dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y1count, &
         complex_type, rbuf, decomp%x1count, &
         complex_type, DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y1cnts, decomp%y1disp, &
         complex_type, rbuf, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_z_to_x_complex_start


  subroutine transpose_z_to_x_complex_wait(handle, src, dst, sbuf, &
       rbuf, decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_zx_complex(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

    return
  end subroutine transpose_z_to_x_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_zx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_zx_real


  subroutine mem_split_zx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_zx_complex


  subroutine mem_merge_zx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%x1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_zx_real


  subroutine mem_merge_zx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%x1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_zx_complex
