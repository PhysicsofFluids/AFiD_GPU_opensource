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

  subroutine transpose_x_to_z_real(src, dst, decomp)

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
    call mem_split_xz_real(src, s1, s2, s3, work1_r, dims(1), &
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
    call mem_merge_xz_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)
    
    return
  end subroutine transpose_x_to_z_real

#ifdef USE_CUDA
  subroutine transpose_x_to_z_real_d(src, dst, decomp)

    implicit none

    real(mytype), device, dimension(:,:,:), intent(IN) :: src
    real(mytype), device, dimension(:,:,:), intent(OUT) :: dst
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

    !print *,nrank, "x->y src ",s1,s2,s3
    !print *,nrank, "x->y dst ",d1,d2,d3


    ! rearrange source array as send buffer
    !call mem_split_xz_real(src, s1, s2, s3, work1_r, dims(1), &
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

       istat = cudaMemcpy2D( work1_r_d(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1, s2*s3 )

    end do

    work1_r = work1_r_d
 
    !compare(work1_r_d,"work1")
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
    !call mem_merge_xz_real(work2_r, d1, d2, d3, dst, dims(1), &
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

       istat = cudaMemcpy2D( dst(1,i1,1), d1*d2, work2_r(pos), d1*(i2-i1+1), d1*(i2-i1+1), d3 )
       if(istat) print *,"ERROR, istat =",istat
    end do


    return
  end subroutine transpose_x_to_z_real_d

  subroutine transpose_x_to_z_complex_d(src, dst, decomp)

    implicit none

    complex(mytype), device, dimension(:,:,:), intent(IN) :: src
    complex(mytype), device, dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror,istat
    integer :: i,j,k,m,i1,i2,pos,iter,dest,sorc

call nvtxStartRange("tranXZ",3)

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
         istat = cudaMemcpy2DAsync( dst(1,1,decomp%z1idx(m)), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
       else
         pos = decomp%w1disp(m) + 1
         istat = cudaMemcpy2DAsync( work1_c_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
         istat = cudaEventRecord( a2a_event(iter+dims(1)), a2a_comp )
         istat = cudaStreamWaitEvent( a2a_d2h, a2a_event(iter+dims(1)), 0)
         istat = cudaMemcpyAsync( work1_c(pos), work1_c_d(pos), decomp%w1cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(1)-1
       m = IEOR(col_rank,iter)
      call MPI_IRECV( work2_c(decomp%z1disp(m)+1), decomp%z1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(1)-1
       m = IEOR(col_rank,iter)
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter)
       call MPI_SEND( work1_c(decomp%w1disp(m)+1), decomp%w1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, ierror)
call nvtxEndRangeAsync
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%z1disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z1idx(m)), work2_c(pos), decomp%z1cnts(m), a2a_h2d )
    end do

    istat = cudaEventRecord( a2a_event(0), 0 )

    else

    ! rearrange source array as send buffer
    do iter=0,dims(1)-1
       dest = mod(col_rank + iter, dims(1))
       m = dest
       !m = IEOR(col_rank,iter)
       if(iter==0) then ! if exchanging with self
         istat = cudaMemcpy2DAsync( dst(1,1,decomp%z1idx(m)), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
       else
         pos = decomp%w1disp(m) + 1
         istat = cudaMemcpy2DAsync( work1_c_d(pos), decomp%x1dist(m), src(decomp%x1idx(m),1,1), s1, decomp%x1dist(m), s2*s3, stream=a2a_comp )
         istat = cudaEventRecord( a2a_event(iter+dims(1)), a2a_comp )
         istat = cudaStreamWaitEvent( a2a_d2h, a2a_event(iter+dims(1)), 0)
         istat = cudaMemcpyAsync( work1_c(pos), work1_c_d(pos), decomp%w1cnts(m), a2a_d2h )
         istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
       end if
    end do

    do iter=1,dims(1)-1
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = sorc
       !m = IEOR(col_rank,iter) 
       call MPI_IRECV( work2_c(decomp%z1disp(m)+1), decomp%z1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(1)-1
       dest = mod(col_rank + iter, dims(1))
       sorc = mod(col_rank - iter + dims(1), dims(1))
       m = dest
       !m = IEOR(col_rank,iter) 
       istat = cudaEventSynchronize( a2a_event(iter) )
call nvtxStartRangeAsync("MPI",iter) 
       call MPI_SEND( work1_c(decomp%w1disp(m)+1), decomp%w1cnts(m), complex_type, m, 0, DECOMP_2D_COMM_COL, ierror)
call nvtxEndRangeAsync
       m = sorc
       call MPI_WAIT(a2a_requests(iter), a2a_status(1,iter), ierror)
       pos = decomp%z1disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z1idx(m)), work2_c(pos), decomp%z1cnts(m), a2a_h2d )
    end do

    istat = cudaEventRecord( a2a_event(0), 0 ) 

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
       pos = decomp%w1disp(m) + 1
       istat = cudaMemcpy2D( work1_c_d(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1,s2*s3 )
    end do

    istat = cudaMemcpy( work1_c, work1_c_d, s1*s2*s3 )

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_c, decomp%w1cnts, decomp%w1disp, &
         complex_type, work2_c, decomp%z1cnts, decomp%z1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)


    istat = cudaMemcpy( dst, work2_c, d1*d2*d3 )
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
       pos = decomp%w1disp(m) + 1
       istat = cudaMemcpy2D( work1_c_d(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1,s2*s3 )
    end do

    istat = cudaMemcpy( work1_c, work1_c_d, s1*s2*s3 )

    ! transpose using MPI_ALLTOALL(V)
    call MPI_ALLTOALLV(work1_c, decomp%w1cnts, decomp%w1disp, &
         complex_type, work2_c, decomp%z1cnts, decomp%z1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)


    istat = cudaMemcpy( dst, work2_c, d1*d2*d3 )
#endif

call nvtxEndRange
    return
  end subroutine transpose_x_to_z_complex_d

#endif

  subroutine transpose_x_to_z_complex(src, dst, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: ierror

      call MPI_Alltoallw(src,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz, &
        dst,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz,MPI_COMM_WORLD,ierror)

    return

  end subroutine transpose_x_to_z_complex



  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_xz_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
  end subroutine mem_split_xz_real


  subroutine mem_split_xz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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
  end subroutine mem_split_xz_complex


  subroutine mem_merge_xz_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
  end subroutine mem_merge_xz_real


  subroutine mem_merge_xz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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
  end subroutine mem_merge_xz_complex
