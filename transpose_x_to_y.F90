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

! This file contains the routines that transpose data from X to Y pencil

!#ifdef USE_CUDA
#if 0
  attributes(global) subroutine mem_split_xy_real_kernel(in,n1,n2,n3,out,iproc,dist,disp)

    implicit none

    integer, value, intent(IN) :: n1,n2,n3
    real(mytype), device, dimension(1:n1,1:n2,1:n3), intent(IN) :: in
    real(mytype), device, dimension(*), intent(OUT) :: out
    integer, value, intent(IN) :: iproc
    integer, device, dimension(0:iproc-1), intent(IN) :: dist
    !TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, device, dimension(0:iproc-1), intent(IN) :: disp
    integer :: i,j,k, m,i1,i2,pos, im

    i = threadIdx%x + (blockIdx%x-1)*blockDim%x;
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y;
    k = blockIdx%z;

    m = 0
    i1 = 1
    i2 = dist(0)

    do im = 0, iproc-2
       if( i > i2 ) then
          m = m + 1
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if
    end do

    pos = 1+disp(m) + (i-i1) + (j-1)*dist(m) + (k-1)*dist(m)*n2


    if(i<=n1 .and. j<=n2 .and. k<=n3) then

    !if(i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) then
    ! write(*,*) "I am thread",i,j,k,"pos = ",pos
    !end if

!    if(i<n1 .and. j<n2 .and. k<n3) then
       out(pos) = in(i,j,k)
    end if
  end subroutine mem_split_xy_real_kernel

  subroutine transpose_x_to_y_real_d(src, dst, decomp)

    implicit none

    real(mytype), device, dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i

    TYPE(dim3) :: tBlock, grid

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(s1)/tBlock%x), ceiling(real(s2)/tBlock%y), s3)
    call mem_split_xy_real_kernel<<<grid,tBlock>>>(src, s1, s2, s3, work1_r_d, dims(1), &
         decomp%x1dist_d, decomp%x1disp_d)

    istat = cudaDeviceSynchronize()
    if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)

    work1_r = work1_r_d

    !do i=5200,5210
    !  write(*,*) 'work1(',i,')= ',work1_r(i)
    !end do

    ! transpose using MPI_ALLTOALL(V) 
    call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, & 
         real_type, work2_r, decomp%y1cnts, decomp%y1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    call mem_merge_xy_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp) 
    
    return
  end subroutine transpose_x_to_y_real_d
#endif

  subroutine transpose_x_to_y_real(src, dst, decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    !print *,nrank, "x->y src ",s1,s2,s3
    !print *,nrank, "x->y dst ",d1,d2,d3


    ! rearrange source array as send buffer
    call mem_split_xy_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%x1dist, decomp)


    !compare(work1_r,"work1")
    !call MPI_BARRIER(MPI_COMM_WORLD)
    !stop
    !return

      !do i=5200,5210
      !  write(*,*) 'work1(',i,')= ',work1_r(i)
      !end do

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, &
         real_type, work2_r, decomp%y1cnts, decomp%y1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    call mem_merge_xy_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)
    
    return
  end subroutine transpose_x_to_y_real

#ifdef USE_CUDA
  subroutine transpose_x_to_y_real_d(src, dst, decomp)

    implicit none

    real(mytype), device, dimension(:,:,:), intent(IN) :: src
    real(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i,j,k,m,i1,i2,pos,dest,sorc,iter

call nvtxStartRange("tranXY",0)

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
         pos = decomp%y1disp(m) + 1
         istat = cudaMemcpy2DAsync( work2_r_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
       else
         pos = decomp%x1disp(m) + 1
         istat = cudaMemcpy2DAsync( work1_r_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
         istat = cudaEventRecord( a2a_event(iter+dims(1)), a2a_comp )
         istat = cudaStreamWaitEvent( a2a_d2h, a2a_event(iter+dims(1)), 0)
         istat = cudaMemcpyAsync( work1_r(pos), work1_r_d(pos), decomp%x1cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(1)-1
       m = IEOR(col_rank,iter)
      call MPI_IRECV( work2_r(decomp%y1disp(m)+1), decomp%y1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(1)-1
       m = IEOR(col_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_r(decomp%x1disp(m)+1), decomp%x1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, ierror)
call nvtxEndRangeAsync
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpyAsync( work2_r_d(pos), work2_r(pos), decomp%y1cnts(m), a2a_h2d )
       istat = cudaEventRecord( a2a_event(iter), a2a_h2d )
    end do

    ! rearrange receive buffer
    do iter=0,dims(1)-1
       m = IEOR(col_rank,iter)
       pos = decomp%y1disp(m) + 1
       istat = cudaStreamWaitEvent( a2a_comp, a2a_event(iter), 0 )
       istat = cudaMemcpy2DAsync( dst(1,decomp%y1idx(m),1), d1*d2, work2_r_d(pos), d1*(decomp%y1dist(m)), d1*(decomp%y1dist(m)), d3, stream=a2a_comp )
    end do

    else

    ! rearrange source array as send buffer
    do iter=0,dims(1)-1
       dest = mod(col_rank + iter, dims(1))
       m = dest
       !m = IEOR(col_rank,iter)
       if(iter==0) then ! if exchanging with self
         pos = decomp%y1disp(m) + 1
         istat = cudaMemcpy2DAsync( work2_r_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
       else
         pos = decomp%x1disp(m) + 1
         istat = cudaMemcpy2DAsync( work1_r_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
         istat = cudaEventRecord( a2a_event(iter+dims(1)), a2a_comp )
         istat = cudaStreamWaitEvent( a2a_d2h, a2a_event(iter+dims(1)), 0)
         istat = cudaMemcpyAsync( work1_r(pos), work1_r_d(pos), decomp%x1cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(1)-1
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = sorc
       !m = IEOR(col_rank,iter)
       call MPI_IRECV( work2_r(decomp%y1disp(m)+1), decomp%y1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(1)-1
       dest = mod(col_rank + iter, dims(1))
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = dest
       !m = IEOR(col_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_r(decomp%x1disp(m)+1), decomp%x1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, ierror)
call nvtxEndRangeAsync
       m = sorc
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpyAsync( work2_r_d(pos), work2_r(pos), decomp%y1cnts(m), a2a_h2d )
       istat = cudaEventRecord( a2a_event(iter), a2a_h2d )
    end do

    ! rearrange receive buffer
    do iter=0,dims(1)-1
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = sorc
       !m = IEOR(col_rank,iter)
       pos = decomp%y1disp(m) + 1
       istat = cudaStreamWaitEvent( a2a_comp, a2a_event(iter), 0 )
       istat = cudaMemcpy2DAsync( dst(1,decomp%y1idx(m),1), d1*d2, work2_r_d(pos), d1*(decomp%y1dist(m)), d1*(decomp%y1dist(m)), d3, stream=a2a_comp )
    end do


#if 0
    ! rearrange source array as send buffer
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if
       pos = decomp%x1disp(m) + 1
       istat = cudaMemcpy2DAsync( work1_r_d(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1, s2*s3, 0 )
       istat = cudaMemcpyAsync( work1_r(pos), work1_r_d(pos), decomp%x1cnts(m), 0 )
       istat = cudaEventRecord( a2a_event(m), 0 ) 
    end do

    !istat = cudaMemcpy( work1_r, work1_r_d, s1*s2*s3)

    do m=0,dims(1)-1
      call MPI_IRECV( work2_r(decomp%y1disp(m)+1), decomp%y1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(m), ierror)
    end do

    do m=0,dims(1)-1
      istat = cudaEventSynchronize( a2a_event(m) )
      call MPI_ISEND( work1_r(decomp%x1disp(m)+1), decomp%x1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(m+dims(1)), ierror)
    end do

    call MPI_WAITALL(2*dims(1), a2a_requests, a2a_status, ierror)


    ! transpose using MPI_ALLTOALL(V)
!    call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, &
!            real_type, work2_r, decomp%y1cnts, decomp%y1disp, &
!            real_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y1dist(m)-1
       end if
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpy2D( dst(1,i1,1), d1*d2, work2_r(pos), d1*(i2-i1+1), d1*(i2-i1+1), d3 )
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
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if
       pos = decomp%x1disp(m) + 1
       istat = cudaMemcpy2D( work1_r_d(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1, s2*s3 )
       if(istat) print *,nrank,"ERROR!",istat
    end do

    istat = cudaMemcpy( work1_r, work1_r_d, s1*s2*s3)
       if(istat) print *,nrank,"ERROR!",istat

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, &
         real_type, work2_r, decomp%y1cnts, decomp%y1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if(ierror) print *,nrank,"A2A err=",ierror

    ! rearrange receive buffer
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y1dist(m)-1
       end if
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpy2D( dst(1,i1,1), d1*d2, work2_r(pos), d1*(i2-i1+1), d1*(i2-i1+1), d3 )
       if(istat) print *,nrank,"ERROR!",istat
    end do
#endif

call nvtxEndRange

    return
  end subroutine transpose_x_to_y_real_d

  subroutine transpose_x_to_y_complex_d(src, dst, decomp)

    implicit none

    complex(mytype), device, dimension(:,:,:), intent(IN) :: src
    complex(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i,j,k,m,i1,i2,pos

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
    !call mem_split_xy_complex(src, s1, s2, s3, work1_c, dims(1), &
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

       istat = cudaMemcpy2D( work1_c_d(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1,s2*s3 )

    end do

    work1_c = work1_c_d


    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, work2_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    !call mem_merge_xy_complex(work2_c, d1, d2, d3, dst, dims(1), &
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

       istat = cudaMemcpy2D( dst(1,i1,1), d1*d2, work2_c(pos), d1*(i2-i1+1),d1*(i2-i1+1), d3 )
       if(istat) print *,"ERROR, istat =",istat
    end do


    return
  end subroutine transpose_x_to_y_complex_d

#endif

  subroutine transpose_x_to_y_complex(src, dst, decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
    call mem_split_xy_complex(src, s1, s2, s3, work1_c, dims(1), &
         decomp%x1dist, decomp)
    
    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, work2_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    call mem_merge_xy_complex(work2_c, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)

    return
  end subroutine transpose_x_to_y_complex



  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_xy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(1:n1,1:n2,1:n3), intent(IN) :: in
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

       pos = decomp%x1disp(m) + 1

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_xy_real


  subroutine mem_split_xy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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

       pos = decomp%x1disp(m) + 1

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_xy_complex


  subroutine mem_merge_xy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(1:n1,1:n2,1:n3), intent(OUT) :: out
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

       pos = decomp%y1disp(m) + 1

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
  end subroutine mem_merge_xy_real


  subroutine mem_merge_xy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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

       pos = decomp%y1disp(m) + 1

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
  end subroutine mem_merge_xy_complex
