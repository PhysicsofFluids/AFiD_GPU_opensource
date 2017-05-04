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

! This file contains the routines that transpose data from Z to Y pencil

  subroutine transpose_z_to_y_real(src, dst, decomp)

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
    work1_p = decomp%ROW_INFO%SND_P
    call mem_split_zy_real(src, s1, s2, s3, work1, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call mem_split_zy_real(src, s1, s2, s3, work1_r, dims(2), &
            decomp%z2dist, decomp)
    end if
#else
    ! note the src array is suitable to be a send buffer
    ! so no split operation needed
#endif
#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%z2cnts_s, decomp%z2disp_s, &
            real_type, work2, decomp%y2cnts_s, decomp%y2disp_s, &
            real_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    if (decomp%even) then
       call MPI_ALLTOALL(src, decomp%z2count, &
            real_type, work2_r, decomp%y2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_r, decomp%z2count, &
            real_type, work2_r, decomp%y2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, &
         real_type, work2_r, decomp%y2cnts, decomp%y2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_zy_real(work2, d1, d2, d3, dst, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_merge_zy_real(work2_r, d1, d2, d3, dst, dims(2), &
         decomp%y2dist, decomp)
#endif
    
    return
  end subroutine transpose_z_to_y_real


#ifdef OCC
  subroutine transpose_z_to_y_real_start(handle, src, dst, sbuf, rbuf)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: ierror

    sbuf = src

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%z2count, real_type, &
         rbuf, decomp%y2count, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%z2cnts, decomp%z2disp, real_type, &
         rbuf, decomp%y2cnts, decomp%y2disp, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_z_to_y_real_start


  subroutine transpose_z_to_y_real_wait(handle, src, dst, sbuf, rbuf)

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
    call mem_merge_zy_real(rbuf, d1, d2, d3, dst, dims(2), &
         decomp%y2dist, decomp)

    return
  end subroutine transpose_z_to_y_real_wait
#endif


  subroutine transpose_z_to_y_complex(src, dst, decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: decomp

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
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

    !print *,nrank,"z->y src ",s1,s2,s3
    !print *,nrank,"z->y dst ",d1,d2,d3

    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%ROW_INFO%SND_P_c
    call mem_split_zy_complex(src, s1, s2, s3, work1, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call mem_split_zy_complex(src, s1, s2, s3, work1_c, dims(2), &
            decomp%z2dist, decomp)
    end if
#else
    ! note the src array is suitable to be a send buffer
    ! so no split operation needed
#endif
#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P_c
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%z2cnts_s, decomp%z2disp_s, &
            complex_type, work2, decomp%y2cnts_s, decomp%y2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    if (decomp%even) then
       call MPI_ALLTOALL(src, decomp%z2count, &
            complex_type, work2_c, decomp%y2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_c, decomp%z2count, &
            complex_type, work2_c, decomp%y2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, &
         complex_type, work2_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_zy_complex(work2, d1, d2, d3, dst, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_merge_zy_complex(work2_c, d1, d2, d3, dst, dims(2), &
         decomp%y2dist, decomp)
#endif

    return
  end subroutine transpose_z_to_y_complex

#ifdef USE_CUDA
  subroutine transpose_z_to_y_complex_d(src, dst, decomp)

    implicit none

    complex(mytype), device, dimension(:,:,:), intent(IN) :: src
    complex(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i,j,k, m,i1,i2, pos, iter, sorc, dest

call nvtxStartRange("tranZY",4)

#ifdef EPA2A
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    if(IAND(dims(2),dims(2)-1)==0) then

    ! rearrange source array as send buffer
    do iter=0,dims(2)-1
       m = IEOR(row_rank,iter)
       if(iter==0) then ! if exchanging with self
         pos = decomp%y2disp(m) + 1
         istat = cudaMemcpyAsync( work1_c_d(pos), src(1,1,decomp%z2idx(m)), decomp%z2cnts(m), a2a_comp )
       else
         pos = decomp%z2disp(m) + 1
         istat = cudaMemcpyAsync( work1_c(pos), src(1,1,decomp%z2idx(m)), decomp%z2cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(2)-1
       m = IEOR(row_rank,iter)
      call MPI_IRECV( work2_c(decomp%y2disp(m)+1), decomp%y2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(2)-1
       m = IEOR(row_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_c(decomp%z2disp(m)+1), decomp%z2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
call nvtxEndRangeAsync
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpyAsync( work1_c_d(pos), work2_c(pos), decomp%y2cnts(m), a2a_h2d )
       istat = cudaEventRecord( a2a_event(iter), a2a_h2d )
    end do

    ! rearrange receive buffer
    do iter=0,dims(2)-1
       m = IEOR(row_rank,iter)
       pos = decomp%y2disp(m) + 1
       istat = cudaStreamWaitEvent( a2a_comp, a2a_event(iter), 0 )
       istat = cudaMemcpy2DAsync( dst(1,decomp%y2idx(m),1), d1*d2, work1_c_d(pos), d1*(decomp%y2dist(m)), d1*(decomp%y2dist(m)), d3, stream=a2a_comp )
    end do

    else

    ! rearrange source array as send buffer
    do iter=0,dims(2)-1
       dest = mod(row_rank + iter, dims(2))
       m = dest
       !m = IEOR(row_rank,iter)
       if(iter==0) then ! if exchanging with self
         pos = decomp%y2disp(m) + 1
         istat = cudaMemcpyAsync( work1_c_d(pos), src(1,1,decomp%z2idx(m)), decomp%z2cnts(m), a2a_comp )
       else
         pos = decomp%z2disp(m) + 1
         istat = cudaMemcpyAsync( work1_c(pos), src(1,1,decomp%z2idx(m)), decomp%z2cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(2)-1
       sorc = mod(row_rank - iter + dims(2), dims(2))
       m = sorc
       !m = IEOR(row_rank,iter)
      call MPI_IRECV( work2_c(decomp%y2disp(m)+1), decomp%y2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(2)-1
       dest = mod(row_rank + iter, dims(2))
       sorc = mod(row_rank - iter + dims(2), dims(2))
       m = dest
       !m = IEOR(row_rank,iter) 
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_c(decomp%z2disp(m)+1), decomp%z2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
call nvtxEndRangeAsync
       m = sorc
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpyAsync( work1_c_d(pos), work2_c(pos), decomp%y2cnts(m), a2a_h2d )
       istat = cudaEventRecord( a2a_event(iter), a2a_h2d )
    end do

    ! rearrange receive buffer
    do iter=0,dims(2)-1
       sorc = mod(row_rank - iter + dims(2), dims(2))
       m = sorc
       !m = IEOR(row_rank,iter) 
       pos = decomp%y2disp(m) + 1
       istat = cudaStreamWaitEvent( a2a_comp, a2a_event(iter), 0 )
       istat = cudaMemcpy2DAsync( dst(1,decomp%y2idx(m),1), d1*d2, work1_c_d(pos), d1*(decomp%y2dist(m)), d1*(decomp%y2dist(m)), d3, stream=a2a_comp )
    end do

#if 0
    istat = cudaMemcpy( work1_c, src, s1*s2*s3 )

    call MPI_ALLTOALLV(work1_c, decomp%z2cnts, decomp%z2disp, &
         complex_type, work2_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)

    ! rearrange receive buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpy2D( dst(1,i1,1), d1*d2, work2_c(pos), d1*(i2-i1+1), d1*(i2-i1+1), d3 )
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

    istat = cudaMemcpy( work1_c, src, s1*s2*s3 )

    call MPI_ALLTOALLV(work1_c, decomp%z2cnts, decomp%z2disp, &
         complex_type, work2_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)

    ! rearrange receive buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpy2D( dst(1,i1,1), d1*d2, work2_c(pos), d1*(i2-i1+1), d1*(i2-i1+1), d3 )
    end do

#endif

call nvtxEndRange
    return
  end subroutine transpose_z_to_y_complex_d
#endif

#ifdef OCC
  subroutine transpose_z_to_y_complex_start(handle, src, dst, sbuf, &
       rbuf)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: ierror

    sbuf = src

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%z2count, &
         complex_type, rbuf, decomp%y2count, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%z2cnts, decomp%z2disp, &
         complex_type, rbuf, decomp%y2cnts, decomp%y2disp, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_z_to_y_complex_start


  subroutine transpose_z_to_y_complex_wait(handle, src, dst, sbuf, &
       rbuf )

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
    call mem_merge_zy_complex(rbuf, d1, d2, d3, dst, dims(2), &
         decomp%y2dist, decomp)

    return
  end subroutine transpose_z_to_y_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_zy_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
       pos = decomp%z2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_zy_real


  subroutine mem_split_zy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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
       pos = decomp%z2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_zy_complex


  subroutine mem_merge_zy_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
       pos = decomp%y2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_zy_real


  subroutine mem_merge_zy_complex(in,n1,n2,n3,out,iproc,dist,decomp)
    
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
       pos = decomp%y2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_zy_complex
