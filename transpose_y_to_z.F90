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

! This file contains the routines that transpose data from Y to Z pencil

  subroutine transpose_y_to_z_real(src, dst, decomp)

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
    call mem_split_yz_real(src, s1, s2, s3, work1, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_split_yz_real(src, s1, s2, s3, work1_r, dims(2), &
         decomp%y2dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            real_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            real_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    if (decomp%even) then
       call MPI_ALLTOALL(work1_r, decomp%y2count, &
            real_type, dst, decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_r, decomp%y2count, &
            real_type, work2_r, decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
         real_type, dst, decomp%z2cnts, decomp%z2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_yz_real(work2, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call mem_merge_yz_real(work2_r, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif
    
    return
  end subroutine transpose_y_to_z_real


#ifdef OCC
  subroutine transpose_y_to_z_real_start(handle, src, dst, sbuf, rbuf, &
       decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN)  :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)

    ! rearrange source array as send buffer
    call mem_split_yz_real(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y2count, real_type, &
         rbuf, decomp%z2count, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, real_type, &
         rbuf, decomp%z2cnts, decomp%z2disp, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_real_start


  subroutine transpose_y_to_z_real_wait(handle, src, dst, sbuf, rbuf, &
       decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_real_wait
#endif


  subroutine transpose_y_to_z_complex(src, dst, decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp


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

    !print *,nrank,"y->z src ",s1,s2,s3
    !print *,nrank,"y->z dst ",d1,d2,d3
    
    ! rearrange source array as send buffer
    call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
         decomp%y2dist, decomp)
    
    ! define receive buffer
    
    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, dst, decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)

    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed

    return
  end subroutine transpose_y_to_z_complex

#ifdef USE_CUDA
  subroutine transpose_y_to_z_complex_d(src, dst, decomp)
#ifdef DEBUG
    use ep_debug
#endif
    implicit none

    complex(mytype), device, dimension(:,:,:), intent(IN) :: src
    complex(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i,j,k,m,i1,i2,pos,dest,sorc,iter

call nvtxStartRange("tranYZ",1)

#ifdef EPA2A

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)


    if(IAND(dims(2),dims(2)-1)==0) then
      !print *,"dims(2)= ",dims(2)," is a power of 2!"

    do iter=0,dims(2)-1
       m = IEOR(row_rank,iter)
       if(iter==0) then ! if exchanging with self
         istat = cudaMemcpy2DAsync( dst(1,1,decomp%z2idx(m)), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_comp )
       else
         pos = decomp%y2disp(m) + 1
         istat = cudaMemcpy2DAsync( work1_c_d(pos), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_comp )
         istat = cudaEventRecord( a2a_event(iter+dims(1)), a2a_comp )
         istat = cudaStreamWaitEvent( a2a_d2h, a2a_event(iter+dims(1)), 0)
         istat = cudaMemcpyAsync( work1_c(pos), work1_c_d(pos), decomp%y2cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(2)-1
       m = IEOR(row_rank,iter)
      call MPI_IRECV( work2_c(decomp%z2disp(m)+1), decomp%z2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(2)-1
       m = IEOR(row_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_c(decomp%y2disp(m)+1), decomp%y2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
call nvtxEndRangeAsync
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%z2disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z2idx(m)), work2_c(pos), decomp%z2cnts(m), a2a_h2d )
    end do

    istat = cudaEventRecord( a2a_event(0), 0)

    else

    do iter=0,dims(2)-1
       dest = mod(row_rank + iter, dims(2))
       m = dest       
       !m = IEOR(row_rank,iter)
       if(iter==0) then ! if exchanging with self
         istat = cudaMemcpy2DAsync( dst(1,1,decomp%z2idx(m)), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_comp )
       else
         pos = decomp%y2disp(m) + 1
         istat = cudaMemcpy2DAsync( work1_c_d(pos), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_comp )
         istat = cudaEventRecord( a2a_event(iter+dims(1)), a2a_comp )
         istat = cudaStreamWaitEvent( a2a_d2h, a2a_event(iter+dims(1)), 0)
         istat = cudaMemcpyAsync( work1_c(pos), work1_c_d(pos), decomp%y2cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(2)-1
       sorc = mod(row_rank - iter + dims(2), dims(2))
       m = sorc
       !m = IEOR(row_rank,iter)
       call MPI_IRECV( work2_c(decomp%z2disp(m)+1), decomp%z2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(2)-1
       dest = mod(row_rank + iter, dims(2))
       sorc = mod(row_rank - iter + dims(2), dims(2))
       m = dest
       !m = IEOR(row_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_c(decomp%y2disp(m)+1), decomp%y2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
call nvtxEndRangeAsync
       m = sorc
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%z2disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z2idx(m)), work2_c(pos), decomp%z2cnts(m), a2a_h2d )
    end do

    istat = cudaEventRecord( a2a_event(0), 0)


#if 0
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if

       pos = decomp%y2disp(m) + 1

       istat = cudaMemcpy2D( work1_c_d(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2, s1*(i2-i1+1), s3 )
       if(istat) print *,nrank,"ERROR a2a memcpy2D!"
    end do

#ifdef DEBUG
    if(row_rank/=999) call compare(work1_c_d,"wk1")
#endif
    istat = cudaMemcpy( work1_c, work1_c_d, s1*s2*s3)
    if(istat) print *,nrank,"ERROR a2a memcpy!"
    ! define receive buffer

    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, work2_c, decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)
    if(ierror) print *,nrank,"ERROR A2A"
    istat = cudaMemcpy( dst, work2_c, d1*d2*d3 )
    if(istat) print *,nrank,"ERROR a2a memcpy!"
#endif
    end if

#else

#if 0
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpy2D( work1_c_d(pos), s1*(decomp%y2dist(m)), src(1,i1,1), s1*s2, s1*(decomp%y2dist(m)), s3 )
       istat = cudaMemcpy( work1_c(pos), work1_c_d(pos), decomp%y2cnts(m))
    end do

    !istat = cudaMemcpy( work1_c, work1_c_d, s1*s2*s3)

    ! define receive buffer

    do m=0,dims(2)-1
      call MPI_IRECV( work2_c(decomp%z2disp(m)+1), decomp%z2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(m), ierror)
    end do

    do m=0,dims(2)-1
      call MPI_ISEND( work1_c(decomp%y2disp(m)+1), decomp%y2cnts(m), complex_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(m+dims(2)), ierror)
    end do

    call MPI_WAITALL(2*dims(2), a2a_requests, a2a_status, ierror)



!    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
!         complex_type, work2_c, decomp%z2cnts, decomp%z2disp, &
!         complex_type, DECOMP_2D_COMM_ROW, ierror)

    do m=0,dims(2)-1
      pos = decomp%z2disp(m) + 1
      istat = cudaMemcpy( dst(1,1,decomp%z2idx(m)), work2_c(pos), decomp%z2cnts(m) )
    end do

#else
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if

       pos = decomp%y2disp(m) + 1

       istat = cudaMemcpy2D( work1_c_d(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2, s1*(i2-i1+1), s3 )
       if(istat) print *,nrank,"ERROR!",istat
    end do
  
    istat = cudaMemcpy( work1_c, work1_c_d, s1*s2*s3)

    ! define receive buffer
if(row_rank<1000) then
    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, work2_c, decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)

       if(ierror) print *,nrank,"ERROR! A2A",ierror
endif
    istat = cudaMemcpy( dst, work2_c, d1*d2*d3 )
       if(istat) print *,nrank,"ERROR!",istat
#endif
#endif

call nvtxEndRange

    return
  end subroutine transpose_y_to_z_complex_d

#endif

#ifdef OCC
  subroutine transpose_y_to_z_complex_start(handle, src, dst, sbuf, &
       rbuf, decomp)

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
    call mem_split_yz_complex(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y2count, &
         complex_type, rbuf, decomp%z2count, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, &
         complex_type, rbuf,decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_complex_start


  subroutine transpose_y_to_z_complex_wait(handle, src, dst, sbuf, &
       rbuf, decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN) :: decomp


    integer :: ierror


    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yz_real


  subroutine mem_split_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yz_complex


  subroutine mem_merge_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_real


  subroutine mem_merge_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)
    
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
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_complex
