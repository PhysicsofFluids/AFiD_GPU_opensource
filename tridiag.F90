#ifdef DEBUG
module ep_debug
  use precision
  use decomp_2d
  use mpih
#ifdef USE_CUDA
  use cudafor
#endif

  implicit none
  integer :: counter=0

  interface compare
     module procedure compare_cpu_1d
     module procedure compare_cpu_3d
     module procedure compare_cpu_complex_3d
#ifdef USE_CUDA
     module procedure compare_gpu_1d
     module procedure compare_cpugpu_1d
     module procedure compare_gpu_3d
     module procedure compare_gpu_complex_3d
     module procedure compare_cpugpu_3d
#endif
  end interface compare

  interface write_plane_jk
    module procedure write_cpu_plane_jk
#ifdef USE_CUDA
    module procedure write_gpu_plane_jk
#endif
  end interface write_plane_jk

contains
  subroutine compare_cpu_1d(A,filename)
    implicit none
    real(fp_kind), dimension(:), intent(IN) :: A
    character (len=*), intent(IN) :: filename
    character (len=4) :: itcount
    character (len=1) :: proc 
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    !print *,"open file..." 
    open(unit=10, status='replace',file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)), form='unformatted')
    !print *,"write file..."
    write(10) A
    !print *,"close file..."
    close(10)
  end subroutine compare_cpu_1d

  subroutine compare_cpu_3d(A,filename)
    implicit none
    real(fp_kind), dimension(:,:,:), intent(IN) :: A
    character (len=*), intent(IN) :: filename
    character (len=4) :: itcount
    character (len=1) :: proc
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    !print *,"open file..." 
    open(unit=10, status='replace', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)), form='unformatted')
    !print *,"write file..."
    write(10) A
    !print *,"close file..."
    close(10)
  end subroutine compare_cpu_3d

  subroutine compare_cpu_complex_3d(A,filename)
    implicit none
    complex(fp_kind), dimension(:,:,:), intent(IN) :: A
    character (len=*), intent(IN) :: filename
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    !print *,"open file..." 
    open(unit=10, status='replace',file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)), form='unformatted')
    !print *,"write file..."
    write(10) A
    !print *,"close file..."
    close(10)
  end subroutine compare_cpu_complex_3d

#ifdef USE_CUDA
  subroutine compare_gpu_1d(A,filename)
    implicit none
    real(fp_kind), dimension(:), device :: A
    character (len=*), intent(IN) :: filename
    real(fp_kind), dimension(:), allocatable :: A_d
    real(fp_kind), dimension(:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf
    integer :: i,imax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    allocate(A_h, source=A)
    allocate(A_d, source=A) 
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
        do i=lbound(A,1),ubound(A,1)
          !print *,nrank,"compare1d i: ",i," cpu: ",A_h(i)
          !print *,nrank,"compare1d i: ",i," gpu: ",A_d(i)
          if(abs(A_h(i)) > 1e-10) then 
            perr = abs(A_h(i) - A_d(i))/abs(A_h(i))*100.d0
            norm = norm + A_h(i)*A_h(i);
            l2normerr = l2normerr + (A_h(i) - A_d(i))*(A_h(i) - A_d(i))
          else
            perr = 0.0d0
          endif
          if(perr>maxerr .and. A_h(i) /= 0.d0) then
            maxerr = perr
            imax = i
          endif
        enddo

 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 l2normerr = buf
 call MPI_ALLREDUCE(maxerr,buf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
if(maxerr == buf) then
  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  write(*,"(I4,A8,A16,2X,ES10.3,A12,ES10.3,A6,I5,A6,2X,E20.14,2X,A6,2X,E20.14)")nrank,filename,"l2norm error",l2normerr,"max error",maxerr,"% at",imax,"cpu=",A_h(imax),"gpu=",A_d(imax)
   else
     if(nrank==0) write(*,"(A8,A16)") filename, "EXACT MATCH"
   endif
endif
  !if(counter<=900) A = A_h
  deallocate(A_h,A_d)
  end subroutine compare_gpu_1d

  subroutine compare_gpu_complex_3d(A,filename)
    implicit none
    complex(fp_kind), dimension(:,:,:), device :: A
    character (len=*), intent(IN) :: filename
    complex(fp_kind), dimension(:,:,:), allocatable :: A_d
    complex(fp_kind), dimension(:,:,:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    allocate(A_h, source=A)
    allocate(A_d, source=A)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    jmax=1
    kmax=1
    do k=lbound(A,3),ubound(A,3)
      do j=lbound(A,2),ubound(A,2)
        do i=lbound(A,1),ubound(A,1)
          if(abs(A_h(i,j,k)) >= 1e-10) then
            perr = abs(A_h(i,j,k) - A_d(i,j,k))/abs(A_h(i,j,k))*100.d0
            norm = norm + abs(A_h(i,j,k)*A_h(i,j,k));
            l2normerr = l2normerr + abs((A_h(i,j,k) - A_d(i,j,k))*(A_h(i,j,k) - A_d(i,j,k)))
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. a_h(i,j,k)/=0.0d0 .and. a_d(i,j,k)/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
            kmax = k
          endif

          !if(counter==1504 .and. perr /= 0.d0) then
#if 0
           if(k<=2.and.j<=2.and.i<3) then
            write(*,"(I4,1X,I3,I3,I3,1X,A4,1X,E20.14,1X,A4,E20.14)") nrank,i,j,k,"cpu=",REAL(A_h(i,j,k))," ",AIMAG(A_h(i,j,k))
            write(*,"(I4,1X,I3,I3,I3,1X,A4,1X,E20.14,1X,A4,E20.14)") nrank,i,j,k,"gpu=",REAL(A_d(i,j,k))," ",AIMAG(A_d(i,j,k))
            !pause
          endif
#endif
        enddo
     enddo
  enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 l2normerr = buf
 !print *,"maxerr=",maxerr
 call MPI_ALLREDUCE(maxerr,buf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 !print *,"merr,buf=",maxerr,buf
if(maxerr == buf) then

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  write(*,"(I4,A8,A16,2X,ES10.3,A12,ES10.3,A6,I5,I5,I5,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
  nrank,filename,"l2norm error",l2normerr,"max error",maxerr,"% at",kmax,jmax,imax,"cpu=",REAL(A_h(imax,jmax,kmax)),AIMAG(A_h(imax,jmax,kmax)),"gpu=",REAL(A_d(imax,jmax,kmax)),AIMAG(A_d(imax,jmax,kmax))
   else
     if(nrank==0) write(*,"(A8,A16)") filename, "EXACT MATCH"
   endif
endif
  !if(counter<=900) A = A_h
  deallocate(A_h,A_d)
  end subroutine compare_gpu_complex_3d

  subroutine compare_cpugpu_1d(A1,A2)
    implicit none
    real(fp_kind), dimension(:), device :: A2
    real(fp_kind), dimension(:) :: A1
    real(fp_kind), dimension(:), allocatable :: A_d
    real(fp_kind), dimension(:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf,tl2n,yl2n,cl2n
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    allocate(A_h, source=A1)
    allocate(A_d, source=A2)
    l2normerr = 0.d0
    cl2n = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
        do i=lbound(A1,1),ubound(A1,1)
          if(abs(A_h(i)) >= 1e-10) then
            perr = abs(A_h(i) - A_d(i))/abs(A_h(i))*100.d0
            norm = norm + A_h(i)*A_h(i);
            yl2n = (A_h(i) - A_d(i))*(A_h(i) - A_d(i)) - cl2n
            tl2n = l2normerr + yl2n
            cl2n = (tl2n-l2normerr) - yl2n
            l2normerr = tl2n
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. a_h(i)/=0.0d0 .and. a_d(i)/=0.0d0) then
            maxerr = perr
            imax = i
          endif
        enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 l2normerr = buf
 !print *,"maxerr=",maxerr
 call MPI_ALLREDUCE(maxerr,buf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 !print *,"merr,buf=",maxerr,buf

if(maxerr == buf) then

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  write(*,"(I4,A16,2X,ES10.3,A12,ES10.3,A6,I5,I5,I5,A6,2X,E20.14,2X,A6,2X,E20.14)") nrank,"l2norm error",l2normerr,"max error",maxerr,"% at",imax,"cpu=",A_h(imax),"gpu=",A_d(imax)
   else
     if(nrank==0) write(*,"(A8,A16)") ,"EXACT MATCH"
   endif
endif
  !if(counter<=900) A = A_h
  deallocate(A_h,A_d)

  !A1 = A2 
  end subroutine compare_cpugpu_1d


  subroutine compare_gpu_3d(A,filename)
    implicit none
    real(fp_kind), dimension(:,:,:), device :: A
    character (len=*), intent(IN) :: filename
    real(fp_kind), dimension(:,:,:), allocatable :: A_d
    real(fp_kind), dimension(:,:,:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf,tl2n,yl2n,cl2n
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    allocate(A_h, source=A)
    allocate(A_d, source=A)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)), form='unformatted')
    read(10) A_h
    close(10)
    l2normerr = 0.d0
    cl2n = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    jmax=1
    kmax=1
    do k=lbound(A,3),ubound(A,3)
      do j=lbound(A,2),ubound(A,2)
        do i=lbound(A,1),ubound(A,1)
          if(abs(A_h(i,j,k)) >= 1e-10) then
            perr = abs(A_h(i,j,k) - A_d(i,j,k))/abs(A_h(i,j,k))*100.d0
            norm = norm + A_h(i,j,k)*A_h(i,j,k);
            yl2n = (A_h(i,j,k) - A_d(i,j,k))*(A_h(i,j,k) - A_d(i,j,k)) - cl2n
            tl2n = l2normerr + yl2n
            cl2n = (tl2n-l2normerr) - yl2n
            l2normerr = tl2n
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. a_h(i,j,k)/=0.0d0 .and. a_d(i,j,k)/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
            kmax = k
          endif
        enddo
     enddo
  enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 l2normerr = buf
 !print *,"maxerr=",maxerr
 call MPI_ALLREDUCE(maxerr,buf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 !print *,"merr,buf=",maxerr,buf

if(maxerr == buf) then

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  write(*,"(I4,A8,A16,2X,ES10.3,A12,ES10.3,A6,I5,I5,I5,A6,2X,E20.14,2X,A6,2X,E20.14)") nrank,filename,"l2norm error",l2normerr,"max error",maxerr,"% at",kmax,jmax,imax,"cpu=",A_h(imax,jmax,kmax),"gpu=",A_d(imax,jmax,kmax)
   else
     if(nrank==0) write(*,"(A8,A16)") filename, "EXACT MATCH"
   endif
endif
  !if(counter<=900) A = A_h
  deallocate(A_h,A_d)
  end subroutine compare_gpu_3d

  subroutine compare_cpugpu_3d(A1,A2)
    implicit none
    real(fp_kind), dimension(:,:,:), device :: A2
    real(fp_kind), dimension(:,:,:) :: A1
    real(fp_kind), dimension(:,:,:), allocatable :: A_d
    real(fp_kind), dimension(:,:,:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf,tl2n,yl2n,cl2n
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank

    allocate(A_h, source=A1)
    allocate(A_d, source=A2)
    l2normerr = 0.d0
    cl2n = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    jmax=1
    kmax=1
    do k=lbound(A1,3),ubound(A1,3)
      do j=lbound(A1,2),ubound(A1,2)
        do i=lbound(A1,1),ubound(A1,1)
          if(abs(A_h(i,j,k)) >= 1e-10) then
            perr = abs(A_h(i,j,k) - A_d(i,j,k))/abs(A_h(i,j,k))*100.d0
            norm = norm + A_h(i,j,k)*A_h(i,j,k);
            yl2n = (A_h(i,j,k) - A_d(i,j,k))*(A_h(i,j,k) - A_d(i,j,k)) - cl2n
            tl2n = l2normerr + yl2n
            cl2n = (tl2n-l2normerr) - yl2n
            l2normerr = tl2n
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. a_h(i,j,k)/=0.0d0 .and. a_d(i,j,k)/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
            kmax = k
          endif
        enddo
     enddo
  enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 l2normerr = buf
 !print *,"maxerr=",maxerr
 call MPI_ALLREDUCE(maxerr,buf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 !print *,"merr,buf=",maxerr,buf

if(maxerr == buf) then

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)
  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  write(*,"(I4,A16,2X,ES10.3,A12,ES10.3,A6,I5,I5,I5,A6,2X,E20.14,2X,A6,2X,E20.14)") nrank,"l2norm error",l2normerr,"max error",maxerr,"% at",imax,jmax,kmax,"cpu=",A_h(imax,jmax,kmax),"gpu=",A_d(imax,jmax,kmax)
   else
     if(nrank==0) write(*,"(A8,A16)") ,"EXACT MATCH"
   endif
endif
  !if(counter<=900) A = A_h
  deallocate(A_h,A_d)

  !A1 = A2 
  end subroutine compare_cpugpu_3d

  subroutine write_cpu_plane_jk(A1,iw,filename)
    implicit none
    real(fp_kind), dimension(:,:,:) :: A1
    integer, intent(IN) ::  iw
    character (len=*), intent(IN) :: filename
    real(fp_kind), dimension(:,:,:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf,tl2n,yl2n,cl2n
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter

    allocate(A_h, source=A1)
    open(unit=10, status='replace',file=filename//".planejk.cpu.dat")
    do k=lbound(A1,3),ubound(A1,3)
      do j=lbound(A1,2),ubound(A1,2)
        write(10,"(ES30.16,1x)", advance='no') A_h(iw,j,k)
      end do
        write(10,*)
    end do
    close(10)

    deallocate(A_h)
  end subroutine write_cpu_plane_jk

  subroutine write_gpu_plane_jk(A1,iw,filename)
    implicit none
    real(fp_kind), dimension(:,:,:), device :: A1
    integer, intent(IN) ::  iw
    character (len=*), intent(IN) :: filename
    real(fp_kind), dimension(:,:,:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf,tl2n,yl2n,cl2n
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter

    allocate(A_h, source=A1)
    open(unit=10, status='replace',file=filename//".planejk.gpu.dat")
    do k=lbound(A1,3),ubound(A1,3)
      do j=lbound(A1,2),ubound(A1,2)
        write(10,"(ES30.16,1x)", advance='no') A_h(iw,j,k)
      end do
        write(10,*)
    end do
    close(10)

    deallocate(A_h)
  end subroutine write_gpu_plane_jk


  subroutine write_cpugpu_2d_ik(A1,A2,jw,filename)
    implicit none
    real(fp_kind), dimension(:,:,:), device :: A2
    real(fp_kind), dimension(:,:,:) :: A1
    integer, intent(IN) ::  jw
    character (len=*), intent(IN) :: filename
    real(fp_kind), dimension(:,:,:), allocatable :: A_d
    real(fp_kind), dimension(:,:,:), allocatable :: A_h
    real(fp_kind) :: maxerr,perr,l2normerr,norm,buf,tl2n,yl2n,cl2n
    integer :: i,j,k,imax,jmax,kmax
    character (len=4) :: itcount
    character (len=1) :: proc

    counter = counter + 1
    write(itcount,'(i4)') counter

    allocate(A_h, source=A1)
    allocate(A_d, source=A2)
    open(unit=10, status='replace',file=filename//".planeik.cpu.dat")
    do k=lbound(A1,3),ubound(A1,3)
      do i=lbound(A1,1),ubound(A1,1)
        write(10,"(ES12.5,1x)", advance='no') A_h(i,jw,k)
      end do
        write(10,*)
    end do
    close(10)

    open(unit=10, status='replace',file=filename//".planeik.gpu.dat")
    do k=lbound(A1,3),ubound(A1,3)
      do i=lbound(A1,1),ubound(A1,1)
        write(10,"(ES12.5, 1x)", advance='no') A_d(i,jw,k)
      end do
        write(10,*)
    end do
    close(10)

  !A1 = A2 
  deallocate(A_h,A_d)
  end subroutine write_cpugpu_2d_ik



#endif

end module ep_debug
#endif

#ifdef USE_CUDA
!#define USE_DBLDBL

module ep_solve
  use precision
  use cudafor
  use decomp_2d, only: a2a_comp
!#ifdef USE_DBLDBL
  use libm
  type dbldbl
    sequence
    real(fp_kind) hi,lo
  end type dbldbl
!#endif

contains

  attributes(global) subroutine tepDgtsv_nopivot_kernel(n,nrhs,a,b,c,x,ldx,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), device, intent(in) :: a,b,c
    real(fp_kind), dimension(nrhs,ldx), device :: x,tmp

    real(fp_kind) :: m
    integer irhs,i

    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x

    if ( irhs <= nrhs) then
    tmp(irhs,1) = c(1)     /b(1)
      x(irhs,1) = x(irhs,1)/b(1)

    do i=2,N 
      m = b(i) - tmp(irhs,i-1)*a(i)
      tmp(irhs,i) = c(i)/m
        x(irhs,i) = (x(irhs,i) - x(irhs,i-1)*a(i))/m
    end do

    do i=n-1,1,-1
      x(irhs,i) = x(irhs,i) - tmp(irhs,i)*x(irhs,i+1)
    end do
    end if

  end subroutine tepDgtsv_nopivot_kernel

  #ifdef USE_HYBRID
  !JR TODO: Using a slow carbon-copy of cuda kernel. Can remove transponse and use openmp to improve.
  subroutine tepDgtsv_nopivot_kernel_cpu(n,nrhs,a,b,c,x,ldx,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), intent(in) :: a,b,c
    real(fp_kind), dimension(nrhs,ldx) :: x,tmp

    real(fp_kind) :: m
    integer irhs,i

    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(irhs, i, m) SHARED(n, nrhs, tmp, a, b, c, x)
    do irhs = 1, nrhs
      tmp(irhs,1) = c(1)     /b(1)
        x(irhs,1) = x(irhs,1)/b(1)

      do i=2,N 
        m = b(i) - tmp(irhs,i-1)*a(i)
        tmp(irhs,i) = c(i)/m
          x(irhs,i) = (x(irhs,i) - x(irhs,i-1)*a(i))/m
      end do

      do i=n-1,1,-1
        x(irhs,i) = x(irhs,i) - tmp(irhs,i)*x(irhs,i+1)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine tepDgtsv_nopivot_kernel_cpu

  subroutine epDgtsv_nopivot_kernel_cpu(n,nrhs,a,b,c,x,ldx,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), intent(in) :: a,b,c
    real(fp_kind), dimension(ldx, nrhs) :: x,tmp

    real(fp_kind) :: m
    integer irhs,i

    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(irhs, i, m) SHARED(n, nrhs, tmp, a, b, c, x)
    do irhs = 1, nrhs
      tmp(1,irhs) = c(1)     /b(1)
        x(1,irhs) = x(1,irhs)/b(1)

      do i=2,N 
        m = b(i) - tmp(i-1,irhs)*a(i)
        tmp(i,irhs) = c(i)/m
          x(i,irhs) = (x(i,irhs) - x(i-1,irhs)*a(i))/m
      end do

      do i=n-1,1,-1
        x(i,irhs) = x(i,irhs) - tmp(i,irhs)*x(i+1,irhs)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine epDgtsv_nopivot_kernel_cpu

  #endif



  attributes(global) subroutine tran_rhs_kernel(idata,odata,n1,n2)
    implicit none
    integer, value, intent(IN) :: n1,n2
    real(fp_kind), dimension(n2,n1), device, intent(IN)  :: idata
    real(fp_kind), dimension(n1,n2), device, intent(OUT) :: odata
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    if(i<=n1 .and. j<=n2) then
      odata(i,j) = idata(j,i)
    end if
  end subroutine tran_rhs_kernel

  subroutine tepDgtsv_nopivot(n,nrhs,a,b,c,x,ldx)
    use decomp_2d, only: rhs_t_d=>work2_r_d
    implicit none
    integer, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), device, intent(IN) :: a,b,c
    real(fp_kind), dimension(LDX,NRHS), device :: x
    type(dim3) :: tBlock, grid, tBlock_tran, grid_tran
 
    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(nrhs)/tBlock_tran%x), ceiling(real(n)/tBlock_tran%y), 1)

#ifdef USE_HYBRID
    call tran_rhs_kernel<<<grid_tran,tBlock_tran,0,a2a_comp>>>(x,rhs_t_d,nrhs,ldx)
#else
    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(x,rhs_t_d,nrhs,ldx)
#endif

    tBlock = dim3(128,1,1)
    grid = dim3(ceiling(real(nrhs)/tBlock%x), 1, 1)

#ifdef USE_HYBRID
    call tepDgtsv_nopivot_kernel<<<grid,tBlock,0,a2a_comp>>>(n,nrhs,a,b,c,rhs_t_d,ldx,x)
#else
    call tepDgtsv_nopivot_kernel<<<grid,tBlock>>>(n,nrhs,a,b,c,rhs_t_d,ldx,x)
#endif

    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(n)/tBlock_tran%x),ceiling(real(nrhs)/tBlock_tran%y), 1)
 
#ifdef USE_HYBRID
    call tran_rhs_kernel<<<grid_tran,tBlock_tran,0,a2a_comp>>>(rhs_t_d,x,ldx,nrhs)
#else
    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(rhs_t_d,x,ldx,nrhs)
#endif

  end subroutine tepDgtsv_nopivot

#ifdef USE_HYBRID
  subroutine tepDgtsv_nopivot_cpu(n,nrhs,a,b,c,x,ldx)
    implicit none
    integer, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), intent(IN) :: a,b,c
    real(fp_kind), dimension(LDX,NRHS) :: x
    real(fp_kind), dimension(NRHS,LDX) :: rhs_t
    integer :: irhs, i

    do irhs = 1,nrhs
      do i = 1, LDX
        rhs_t(irhs, i) = x(i, irhs)
      end do 
    end do

    call tepDgtsv_nopivot_kernel_cpu(n,nrhs,a,b,c,rhs_t,ldx,x)

    do irhs = 1,nrhs
      do i = 1, LDX
        x(i, irhs) = rhs_t(irhs, i)
      end do 
    end do


  end subroutine tepDgtsv_nopivot_cpu

  subroutine epDgtsv_nopivot_cpu(n,nrhs,a,b,c,x,ldx)
    implicit none
    integer, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), intent(IN) :: a,b,c
    real(fp_kind), dimension(NRHS,LDX) :: x
    real(fp_kind), dimension(NRHS,LDX) :: tmp
    integer :: irhs, i

    call epDgtsv_nopivot_kernel_cpu(n,nrhs,a,b,c,x,ldx,tmp)


  end subroutine epDgtsv_nopivot_cpu

#endif

  attributes(global) subroutine epDgtsv_nopivot_kernel(n,nrhs,a,b,c,x,ldx,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), device, intent(in) :: a,b,c
    real(fp_kind), dimension(ldx,nrhs), device :: x,tmp

    real(fp_kind) :: m
    integer irhs,i

    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x

    if ( irhs <= nrhs) then
    tmp(1,irhs) = c(1)     /b(1)
      x(1,irhs) = x(1,irhs)/b(1)

    do i=2,N
      m = b(i) - tmp(i-1,irhs)*a(i)
      tmp(i,irhs) = c(i)/m
        x(i,irhs) = (x(i,irhs) - x(i-1,irhs)*a(i))/m
    end do

    do i=n-1,1,-1
      x(i,irhs) = x(i,irhs) - tmp(i,irhs)*x(i+1,irhs)
    end do
    end if

  end subroutine epDgtsv_nopivot_kernel

  subroutine epDgtsv_nopivot(n,nrhs,a,b,c,x,ldx,tmp)
    implicit none
    integer, intent(in) :: n,nrhs,ldx
    real(fp_kind), dimension(*), device, intent(IN) :: a,b,c
    real(fp_kind), dimension(LDX,NRHS), device :: x,tmp
    type(dim3) :: tBlock, grid

    tBlock = dim3(128,1,1)
    grid = dim3(ceiling(real(nrhs)/tBlock%x), 1, 1)

    call epDgtsv_nopivot_kernel<<<grid,tBlock>>>(n,nrhs,a,b,c,x,ldx,tmp)

  end subroutine epDgtsv_nopivot

  attributes(global) subroutine tepZgtsv_pressure_kernel(n,nrhs,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,tmp,joff,koff)
    implicit none
    integer, value, intent(in) :: n,nrhs,ny,nz,ldx,joff,koff
    real(fp_kind), dimension(:), device, intent(in) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(nrhs,ldx), device :: x
!pgi$ ignore_tkr (t) tmp
    real(fp_kind), dimension(nrhs,ldx), device :: tmp

    real(fp_kind) :: m,acphT_b,a,c
    integer :: irhs,i,j,k

    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = mod(irhs-1,ny) + 1 + joff - 1
    k = (irhs+ny-1)/ny + koff - 1

    if ( irhs <= nrhs) then

    acphT_b = 1.0d0/(acphk(1)-ak2(j)-ak1(k))
    tmp(irhs,1) = apphk(1)*acphT_b
      x(irhs,1) = x(irhs,1)*acphT_b

    do i=2,N
      acphT_b = 1.0d0/(acphk(i)-ak2(j)-ak1(k))
      a = amphk(i)*acphT_b
      c = apphk(i)*acphT_b
      ! When using the FMA instruction, m could be zero
      !m = 1.0d0 - tmp(irhs,i-1)*a
      ! explicit use of rounded multiply
      m = 1.0d0 - __dmul_rn( tmp(irhs,i-1),a)
      tmp(irhs,i) = c/m
        x(irhs,i) = (x(irhs,i)*acphT_b - x(irhs,i-1)*a)/m
    end do

    do i=n-1,1,-1
      x(irhs,i) = x(irhs,i) - tmp(irhs,i)*x(irhs,i+1)
    end do
    end if

  end subroutine tepZgtsv_pressure_kernel

  attributes(global) subroutine tran_complex_rhs_kernel(idata,odata,n1,n2)
    implicit none
    integer, value, intent(IN) :: n1,n2
    complex(fp_kind), dimension(n2,n1), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n1,n2), device, intent(OUT) :: odata
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    if(i<=n1 .and. j<=n2) then
      odata(i,j) = idata(j,i)
    end if
  end subroutine tran_complex_rhs_kernel


  subroutine tepZgtsv_pressure(n, ny, nz, amphk, acphk, apphk, ak1, ak2, x, ldx, joff, koff, work)
    !use fftw_params, only: dphc_t_d
    implicit none
    integer, intent(in) :: n,ny,nz,ldx,joff,koff
    real(fp_kind), dimension(:), device, intent(IN) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(ny*nz,n), device :: work
    complex(fp_kind), dimension(n,ny*nz), device :: x
    type(dim3) :: tBlock, grid, tBlock_tran, grid_tran
    integer :: nrhs

    nrhs = ny*nz

    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(nrhs)/tBlock_tran%x),ceiling(real(n)/tBlock_tran%y), 1)

    call tran_complex_rhs_kernel<<<grid_tran,tBlock_tran>>>(x,work,nrhs,ldx)

    tBlock = dim3(128,1,1)
    grid = dim3(ceiling(real(nrhs)/tBlock%x), 1, 1)

    call tepZgtsv_pressure_kernel<<<grid,tBlock>>>(n,nrhs,ny,nz,amphk,acphk,apphk,ak1,ak2,work,ldx,x,joff,koff)

    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(n)/tBlock_tran%x),ceiling(real(nrhs)/tBlock_tran%y), 1)

    call tran_complex_rhs_kernel<<<grid_tran,tBlock_tran>>>(work,x,ldx,nrhs)

  end subroutine tepZgtsv_pressure

  attributes(device) type(dbldbl) function add_dbldbl(a,b) !add=a+b
    implicit none
    type(dbldbl) a, b
    real(fp_kind) t1, t2, t3, t4, t5, e
    t1 = __dadd_rn (a%hi, b%hi)
    t2 = __dadd_rn (t1, -a%hi)
    t3 = __dadd_rn (__dadd_rn (a%hi, t2 - t1), __dadd_rn (b%hi, -t2))
    t4 = __dadd_rn (a%lo, b%lo)
    t2 = __dadd_rn (t4, -a%lo)
    t5 = __dadd_rn (__dadd_rn (a%lo, t2 - t4), __dadd_rn (b%lo, -t2))
    t3 = __dadd_rn (t3, t4)
    t4 = __dadd_rn (t1, t3)
    t3 = __dadd_rn (t1 - t4, t3)
    t3 = __dadd_rn (t3, t5)
    e = __dadd_rn (t4, t3)
    add_dbldbl%hi = e
    add_dbldbl%lo = __dadd_rn (t4 - e, t3)
  end function add_dbldbl

  attributes(device) type(dbldbl) function sub_dbldbl(a,b) !sub=a-b 
    implicit none
    type(dbldbl) a, b 
    real(fp_kind) t1, t2, t3, t4, t5, e
    t1 = __dadd_rn (a%hi, -b%hi)
    t2 = __dadd_rn (t1, -a%hi);
    t3 = __dadd_rn (__dadd_rn (a%hi, t2 - t1), - __dadd_rn (b%hi, t2))
    t4 = __dadd_rn (a%lo, -b%lo)
    t2 = __dadd_rn (t4, -a%lo)
    t5 = __dadd_rn (__dadd_rn (a%lo, t2 - t4), - __dadd_rn (b%lo, t2))
    t3 = __dadd_rn (t3, t4)
    t4 = __dadd_rn (t1, t3)
    t3 = __dadd_rn (t1 - t4, t3)
    t3 = __dadd_rn (t3, t5)
    e = __dadd_rn (t4, t3)
    sub_dbldbl%hi = e
    sub_dbldbl%lo = __dadd_rn (t4 - e, t3)
  end function sub_dbldbl

  attributes(device) type(dbldbl) function mul_dbldbl(a,b) !mul=a*b 
    implicit none
    type(dbldbl) a, b
    type(dbldbl) t
    real(fp_kind) e
    t%hi = __dmul_rn (a%hi, b%hi)
    t%lo = fma (a%hi, b%hi,-t%hi)
    t%lo = fma (a%lo, b%lo, t%lo)
    t%lo = fma (a%hi, b%lo, t%lo)
    t%lo = fma (a%lo, b%hi, t%lo)
    e    = __dadd_rn (t%hi, t%lo)
    mul_dbldbl%hi = e 
    mul_dbldbl%lo = __dadd_rn (t%hi - e, t%lo)
  end function mul_dbldbl

  attributes(device) type(dbldbl) function div_dbldbl(a,b) !div=a/b 
    implicit none
    type(dbldbl) a, b
    type(dbldbl) t
    real(fp_kind) e, r
    r = 1.0d0 / b%hi;
    t%hi = __dmul_rn (a%hi, r)
    e    = fma (b%hi, -t%hi, a%hi)
    t%hi = fma (r, e, t%hi)
    t%lo = fma (b%hi, -t%hi, a%hi)
    t%lo = __dadd_rn (a%lo, t%lo)
    t%lo = fma (b%lo, -t%hi, t%lo)
    e    = __dmul_rn (r, t%lo)
    t%lo = fma (b%hi, -e, t%lo)
    t%lo = fma (r, t%lo, e)
    e    = __dadd_rn (t%hi, t%lo)
    div_dbldbl%hi = e
    div_dbldbl%lo = __dadd_rn (t%hi - e, t%lo)
  end function div_dbldbl

  attributes(device) type(dbldbl) function neg_dbldbl(a) !neg=-a 
    implicit none
    type(dbldbl) a
    neg_dbldbl%hi = -a%hi
    neg_dbldbl%lo = -a%lo
  end function neg_dbldbl

  attributes(device) type(dbldbl) function mul_double_to_dbldbl(a,b) !mul=a*b 
    implicit none
    real(fp_kind) a,b
    real(fp_kind) c
    c = __dmul_rn(a,b)
    mul_double_to_dbldbl%hi = c
    mul_double_to_dbldbl%lo = fma(a, b, -c)
  end function mul_double_to_dbldbl

  attributes(device) type(dbldbl) function add_double_to_dbldbl(a,b) !add=a+b 
    implicit none
    real(fp_kind) a,b
    real(fp_kind) t1, t2, c
    c  = __dadd_rn (a,  b )
    t1 = __dadd_rn (c, -a )
    t2 = __dadd_rn (c, -t1)
    t1 = __dadd_rn (b, -t1)
    t2 = __dadd_rn (a, -t2)
    add_double_to_dbldbl%hi = c
    add_double_to_dbldbl%lo = __dadd_rn (t1, t2)
  end function add_double_to_dbldbl

  attributes(device) real(fp_kind) function get_dbldbl_hi(a) !get=a%hi 
    implicit none
    type(dbldbl) a
    get_dbldbl_hi = a%hi
  end function get_dbldbl_hi

  attributes(device) real(fp_kind) function get_dbldbl_lo(a) !get=a%lo
    implicit none
    type(dbldbl) a
    get_dbldbl_lo = a%lo
  end function get_dbldbl_lo

  attributes(device) type(dbldbl) function make_dbldbl(a,b) !make%hi=a make%lo=b
    implicit none
    real(fp_kind) a, b
    make_dbldbl%hi = a
    make_dbldbl%lo = b
  end function make_dbldbl


  attributes(global) subroutine epZgtsv_pressure_piv_kernel_dd(n,nrhs,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,DLHI,DLLO,DHI,DLO,DUHI,DULO,XRHI,XRLO,XIHI,XILO,joff,koff)
    implicit none
    integer, value, intent(in) :: n,nrhs,ny,nz,ldx,joff,koff
    real(fp_kind), dimension(:), device, intent(in) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(ldx,nrhs), device :: x
!    real(fp_kind), dimension(ldx,nrhs), device :: DLHI,DLLO,DHI,DLO,DUHI,DULO,XRHI,XRLO,XIHI,XILO
    real(fp_kind), dimension(nrhs,ldx), device :: DLHI,DLLO,DHI,DLO,DUHI,DULO,XRHI,XRLO,XIHI,XILO
    type(dbldbl) :: acphT_b,fact,duip,dlip
    real(fp_kind) :: tempxrhi,tempxrlo,tempxihi,tempxilo
    integer :: irhs,i,j,k
!    j = threadIdx%x + (blockIdx%x-1)*blockDim%x
!    k = threadIdx%y + (blockIdx%y-1)*blockDim%y
!    irhs = (k-1)*ny + j
    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = mod(irhs-1,ny) + 1 
    k = (irhs+ny-1)/ny

    if ( irhs <= nrhs .and. j<=ny .and. k<=nz) then
        j = j + joff - 1
        k = k + koff - 1
        DHI(irhs,1)=1.0d0
        DLO(irhs,1)=0.0d0
        acphT_b%hi = get_dbldbl_hi( div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl( make_dbldbl(acphk(1),0.0d0), add_double_to_dbldbl(ak2(j),ak1(k)) ) ) )
        acphT_b%lo = get_dbldbl_lo( div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl( make_dbldbl(acphk(1),0.0d0), add_double_to_dbldbl(ak2(j),ak1(k)) ) ) )
        !acphT_b = make_dbldbl(temphi,templo)
!1.0d0/(acphk(1)-ak2(j)-ak1(k)) !div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl( make_dbldbl(acphk(1),0.0d0), add_double_to_dbldbl(ak2(j),ak1(k)) ) )
        !if(j==1.and.k==1) print *,"acphT_b at 1,1,1=",get_dbldbl_hi(acphT_b),get_dbldbl_lo(acphT_b)
        DUHI(irhs,1) = get_dbldbl_hi( mul_dbldbl( make_dbldbl(apphk(1),0.0d0), acphT_b) )!apphk(1)*acphT_b !get_dbldbl_hi( mul_dbldbl( make_dbldbl(apphk(1),0.0d0), acphT_b) )
        DULO(irhs,1) = get_dbldbl_lo( mul_dbldbl( make_dbldbl(apphk(1),0.0d0), acphT_b) )
        XRHI(irhs,1) = get_dbldbl_hi( mul_dbldbl( acphT_b, make_dbldbl(REAL (x(1,irhs)),0.0d0) ) )!acphT_b*REAL(x(1,irhs))!get_dbldbl_hi( mul_dbldbl( acphT_b, make_dbldbl(REAL (x(1,irhs)),0.0d0) ) )
        XRLO(irhs,1) = get_dbldbl_lo( mul_dbldbl( acphT_b, make_dbldbl(REAL (x(1,irhs)),0.0d0) ) )
        XIHI(irhs,1) = get_dbldbl_hi( mul_dbldbl( acphT_b, make_dbldbl(AIMAG(x(1,irhs)),0.0d0) ) )!acphT_b*AIMAG(x(1,irhs))!get_dbldbl_hi( mul_dbldbl( acphT_b, make_dbldbl(AIMAG(x(1,irhs)),0.0d0) ) )
        XILO(irhs,1) = get_dbldbl_lo( mul_dbldbl( acphT_b, make_dbldbl(AIMAG(x(1,irhs)),0.0d0) ) )
        !x (1,irhs) = DCMPLX( XRHI(irhs,1), XIHI(irhs,1) )
        !if(j==1.and.k==1) print *,"x 1,1,1 =",REAL(x(1,irhs)),AIMAG(x(1,irhs))
        do i=1,n-1

!#define BUG_REPRO
#ifdef BUG_REPRO
          acphT_b    =                div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl(add_double_to_dbldbl(acphk(i+1),-ak2(j)),make_dbldbl(ak1(k),0.d0) ) )
#else
          acphT_b%hi = get_dbldbl_hi( div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl(add_double_to_dbldbl(acphk(i+1),-ak2(j)),make_dbldbl(ak1(k),0.d0) ) ) )
          acphT_b%lo = get_dbldbl_lo( div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl(add_double_to_dbldbl(acphk(i+1),-ak2(j)),make_dbldbl(ak1(k),0.d0) ) ) )
#endif

          !acphT_b = make_dbldbl(temphi, templo)
          !acphT_b%hi = temphi
          !acphT_b%lo = templo
          !if(j==1.and.k==1) print*,"ac",get_dbldbl_hi(acphT_b),get_dbldbl_lo(acphT_b),acphT_b%hi,acphT_b%lo
          !acphT_b = tempxrhi!div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl(make_dbldbl(acphk(i+1),0.0d0), add_double_to_dbldbl(ak2(j),ak1(k)) ) ) )
!1.0d0/(acphk(i+1)-ak2(j)-ak1(k))!div_dbldbl( make_dbldbl(1.0d0,0.0d0), sub_dbldbl(make_dbldbl(acphk(i+1),0.0d0), add_double_to_dbldbl(ak2(j),ak1(k)) ) )
          duip%hi = get_dbldbl_hi( mul_dbldbl( make_dbldbl(apphk(i+1),0.0d0), acphT_b ) )!apphk(i+1)*acphT_b!mul_dbldbl( make_dbldbl(apphk(i+1),0.0d0), acphT_b )
          duip%lo = get_dbldbl_lo( mul_dbldbl( make_dbldbl(apphk(i+1),0.0d0), acphT_b ) )
          dlip%hi = get_dbldbl_hi( mul_dbldbl( make_dbldbl(amphk(i+1),0.0d0), acphT_b ) )!amphk(i+1)*acphT_b!mul_dbldbl( make_dbldbl(amphk(i+1),0.0d0), acphT_b )
          dlip%lo = get_dbldbl_lo( mul_dbldbl( make_dbldbl(amphk(i+1),0.0d0), acphT_b ) )
          if( abs( DHI(irhs,i) ) >= abs( dlip%hi ) ) then
            fact%hi   = get_dbldbl_hi( div_dbldbl( dlip, make_dbldbl( DHI(irhs,i),DLO(irhs,i)) ) )!dlip / DHI(irhs,i)!div_dbldbl( dlip, make_dbldbl( DHI(irhs,i),DLO(irhs,i)) )
            fact%lo   = get_dbldbl_lo( div_dbldbl( dlip, make_dbldbl( DHI(irhs,i),DLO(irhs,i)) ) )
            DHI(irhs,i+1)  =  get_dbldbl_hi( sub_dbldbl( make_dbldbl(1.0d0,0.0d0), mul_dbldbl( fact, make_dbldbl(DUHI(irhs,i),DULO(irhs,i))) ) )
!1.0d0 - fact*DUHI(irhs,i)!get_dbldbl_hi( sub_dbldbl( make_dbldbl(1.0d0,0.0d0), mul_dbldbl(fact,make_dbldbl(DUHI(irhs,i),DULO(irhs,i))) ) )
            DLO(irhs,i+1)  =  get_dbldbl_lo( sub_dbldbl( make_dbldbl(1.0d0,0.0d0), mul_dbldbl( fact, make_dbldbl(DUHI(irhs,i),DULO(irhs,i))) ) )
            !if(j==1.and.k==1) print *,"D",i,DHI(irhs,i+1),DLO(irhs,i+1)
            
            tempxrhi = get_dbldbl_hi(sub_dbldbl( mul_dbldbl(make_dbldbl(REAL (x(i+1,irhs)),0.0d0), acphT_b), mul_dbldbl(fact,make_dbldbl(XRHI(irhs,i),XRLO(irhs,i))) ))
            tempxrlo = get_dbldbl_lo(sub_dbldbl( mul_dbldbl(make_dbldbl(REAL (x(i+1,irhs)),0.0d0), acphT_b), mul_dbldbl(fact,make_dbldbl(XRHI(irhs,i),XRLO(irhs,i))) ))
            tempxihi = get_dbldbl_hi(sub_dbldbl( mul_dbldbl(make_dbldbl(AIMAG(x(i+1,irhs)),0.0d0), acphT_b), mul_dbldbl(fact,make_dbldbl(XIHI(irhs,i),XILO(irhs,i))) ))
            tempxilo = get_dbldbl_lo(sub_dbldbl( mul_dbldbl(make_dbldbl(AIMAG(x(i+1,irhs)),0.0d0), acphT_b), mul_dbldbl(fact,make_dbldbl(XIHI(irhs,i),XILO(irhs,i))) ))
            XRHI(irhs,i+1) = tempxrhi
            XRLO(irhs,i+1) = tempxrlo
            XIHI(irhs,i+1) = tempxihi
            XILO(irhs,i+1) = tempxilo
            !x (i+1,irhs) = DCMPLX( XRHI(irhs,i+1), XIHI(irhs,i+1) )
            dlip%hi = 0.0d0!make_dbldbl( 0.d0, 0.d0 )
            dlip%lo = 0.0d0
          else
!#if 0
            !print *,"pivot in row: ",i,j,k
            fact%hi = get_dbldbl_hi( div_dbldbl( make_dbldbl(DHI(irhs,i),DLO(irhs,i)), dlip ) )
            fact%lo = get_dbldbl_lo( div_dbldbl( make_dbldbl(DHI(irhs,i),DLO(irhs,i)), dlip ) )
            DHI(irhs,i) = dlip%hi
            DLO(irhs,i) = dlip%lo
            DHI(irhs,i+1) = get_dbldbl_hi( sub_dbldbl( make_dbldbl( DUHI(irhs,i),DULO(irhs,i)), fact ) )
            DLO(irhs,i+1) = get_dbldbl_lo( sub_dbldbl( make_dbldbl( DUHI(irhs,i),DULO(irhs,i)), fact ) )
            dlip%hi = duip%hi
            dlip%lo = duip%lo
            duip%hi = get_dbldbl_hi( mul_dbldbl( neg_dbldbl(fact), duip ) )
            duip%lo = get_dbldbl_lo( mul_dbldbl( neg_dbldbl(fact), duip ) )
            DUHI(irhs,i) = 1.0d0
            DULO(irhs,i) = 0.0d0
            tempxrhi = XRHI(irhs,i)
            tempxrlo = XRLO(irhs,i)
            tempxihi = XIHI(irhs,i)
            tempxilo = XILO(irhs,i)
            XRHI(irhs,i) = get_dbldbl_hi( mul_dbldbl(make_dbldbl(REAL (x(i+1,irhs)),0.0d0),acphT_b) )
            XRLO(irhs,i) = get_dbldbl_lo( mul_dbldbl(make_dbldbl(REAL (x(i+1,irhs)),0.0d0),acphT_b) )
            XIHI(irhs,i) = get_dbldbl_hi( mul_dbldbl(make_dbldbl(AIMAG(x(i+1,irhs)),0.0d0),acphT_b) )
            XILO(irhs,i) = get_dbldbl_lo( mul_dbldbl(make_dbldbl(AIMAG(x(i+1,irhs)),0.0d0),acphT_b) )
            !x (i,irhs) = DCMPLX( XRHI(irhs,i), XIHI(irhs,i) )
            XRHI(irhs,i+1) = get_dbldbl_hi( sub_dbldbl( make_dbldbl(tempxrhi,tempxrlo), mul_dbldbl( fact, make_dbldbl( XRHI(irhs,i), XRLO(irhs,i) ) ) ) )
            XRLO(irhs,i+1) = get_dbldbl_lo( sub_dbldbl( make_dbldbl(tempxrhi,tempxrlo), mul_dbldbl( fact, make_dbldbl( XRHI(irhs,i), XRLO(irhs,i) ) ) ) )
            XIHI(irhs,i+1) = get_dbldbl_hi( sub_dbldbl( make_dbldbl(tempxihi,tempxilo), mul_dbldbl( fact, make_dbldbl( XIHI(irhs,i), XILO(irhs,i) ) ) ) )
            XILO(irhs,i+1) = get_dbldbl_lo( sub_dbldbl( make_dbldbl(tempxihi,tempxilo), mul_dbldbl( fact, make_dbldbl( XIHI(irhs,i), XILO(irhs,i) ) ) ) ) 
            !x (i+1,irhs) = DCMPLX( XRHI(irhs,i+1), XIHI(irhs,i+1) )
         endif
!#endif
         DLHI(irhs,i+1) = dlip%hi
         DLLO(irhs,i+1) = dlip%lo
         DUHI(irhs,i+1) = duip%hi
         DULO(irhs,i+1) = duip%lo
        enddo

        !if(j==1.and.k==1) print*,"XoDo",XRHI(irhs,n),XRLO(irhs,n),DHI(irhs,n),DLO(irhs,n)
        tempxrhi = get_dbldbl_hi( div_dbldbl( make_dbldbl(XRHI(irhs,n),XRLO(irhs,n)), make_dbldbl(DHI(irhs,n),DLO(irhs,n)) ) )
!XRHI(n,irhs)/DHI(n,irhs)!get_dbldbl_hi( div_dbldbl( make_dbldbl(XRHI(n,irhs),XRLO(n,irhs)), make_dbldbl(DHI(n,irhs),DLO(n,irhs)) ) )
        tempxrlo = get_dbldbl_lo( div_dbldbl( make_dbldbl(XRHI(irhs,n),XRLO(irhs,n)), make_dbldbl(DHI(irhs,n),DLO(irhs,n)) ) )
        tempxihi = get_dbldbl_hi( div_dbldbl( make_dbldbl(XIHI(irhs,n),XILO(irhs,n)), make_dbldbl(DHI(irhs,n),DLO(irhs,n)) ) ) 
!XIHI(n,irhs)/DHI(n,irhs)!get_dbldbl_hi( div_dbldbl( make_dbldbl(XIHI(n,irhs),XILO(n,irhs)), make_dbldbl(DHI(n,irhs),DLO(n,irhs)) ) )
        tempxilo = get_dbldbl_lo( div_dbldbl( make_dbldbl(XIHI(irhs,n),XILO(irhs,n)), make_dbldbl(DHI(irhs,n),DLO(irhs,n)) ) )
        XRHI(irhs,n) = tempxrhi
        XRLO(irhs,n) = tempxrlo
        XIHI(irhs,n) = tempxihi
        XILO(irhs,n) = tempxilo
        x (n,irhs) = DCMPLX( tempxrhi, tempxihi )

        !if(j==1.and.k==1) print *,"x n,1,1=",REAL(x(n,irhs)),AIMAG(x(n,irhs)),DHI(irhs,n),DLO(irhs,n)

tempxrhi = get_dbldbl_hi(div_dbldbl(sub_dbldbl(make_dbldbl(XRHI(irhs,n-1),XRLO(irhs,n-1)),mul_dbldbl(make_dbldbl(DUHI(irhs,n-1),DULO(irhs,n-1)), & 
                               make_dbldbl(XRHI(irhs,n),XRLO(irhs,n)))),make_dbldbl(DHI(irhs,n-1),DLO(irhs,n-1)))) 
!(XRHI(irhs,n-1)-DUHI(irhs,n-1)*XRHI(irhs,n))/DHI(irhs,n-1)!get_dbldbl_hi(div_dbldbl(sub_dbldbl(make_dbldbl(XRHI(irhs,n-1),XRLO(irhs,n-1)),mul_dbldbl(make_dbldbl(DUHI(irhs,n-1),DULO(irhs,n-1)), &
                              !make_dbldbl(XRHI(irhs,n),XRLO(irhs,n)))),make_dbldbl(DHI(irhs,n-1),DLO(irhs,n-1)))) 
tempxrlo = get_dbldbl_lo(div_dbldbl(sub_dbldbl(make_dbldbl(XRHI(irhs,n-1),XRLO(irhs,n-1)),mul_dbldbl(make_dbldbl(DUHI(irhs,n-1),DULO(irhs,n-1)), &
                               make_dbldbl(XRHI(irhs,n),XRLO(irhs,n)))),make_dbldbl(DHI(irhs,n-1),DLO(irhs,n-1))))
tempxihi = get_dbldbl_hi(div_dbldbl(sub_dbldbl(make_dbldbl(XIHI(irhs,n-1),XILO(irhs,n-1)),mul_dbldbl(make_dbldbl(DUHI(irhs,n-1),DULO(irhs,n-1)), & 
                               make_dbldbl(XIHI(irhs,n),XILO(irhs,n)))),make_dbldbl(DHI(irhs,n-1),DLO(irhs,n-1))))
!(XIHI(irhs,n-1)-DUHI(irhs,n-1)*XIHI(irhs,n))/DHI(irhs,n-1)!get_dbldbl_hi(div_dbldbl(sub_dbldbl(make_dbldbl(XIHI(irhs,n-1),XILO(irhs,n-1)),mul_dbldbl(make_dbldbl(DUHI(irhs,n-1),DULO(irhs,n-1)), &
                              !make_dbldbl(XIHI(irhs,n),XILO(irhs,n)))),make_dbldbl(DHI(irhs,n-1),DLO(irhs,n-1))))   
tempxilo = get_dbldbl_lo(div_dbldbl(sub_dbldbl(make_dbldbl(XIHI(irhs,n-1),XILO(irhs,n-1)),mul_dbldbl(make_dbldbl(DUHI(irhs,n-1),DULO(irhs,n-1)), &
                               make_dbldbl(XIHI(irhs,n),XILO(irhs,n)))),make_dbldbl(DHI(irhs,n-1),DLO(irhs,n-1))))
        XRHI(irhs,n-1) = tempxrhi
        XRLO(irhs,n-1) = tempxrlo
        XIHI(irhs,n-1) = tempxihi
        XILO(irhs,n-1) = tempxilo
        x (n-1,irhs) = DCMPLX( tempxrhi, tempxihi )

        !if(j==1.and.k==1) print *,"x n-1,1,1 =",REAL(x(n-1,irhs)),AIMAG(x(n-1,irhs)), DHI(irhs,n-1), DLO(irhs,n-1), DUHI(irhs,n-1), DULO(irhs,n-1)

        do i=n-2,1,-1
tempxrhi = get_dbldbl_hi( div_dbldbl( sub_dbldbl(make_dbldbl(XRHI(irhs,i),XRLO(irhs,i)), add_dbldbl(mul_dbldbl(make_dbldbl(DUHI(irhs,i),DULO(irhs,i)),make_dbldbl(XRHI(irhs,i+1),XRLO(irhs,i+1))), &
                              mul_dbldbl(make_dbldbl(DLHI(irhs,i+1),DLLO(irhs,i+1)),make_dbldbl(XRHI(irhs,i+2),XRLO(irhs,i+2))))), make_dbldbl(DHI(irhs,i),DLO(irhs,i)) ) )

!(XRHI(irhs,i)-DUHI(irhs,i)*XRHI(irhs,i+1)-DLHI(irhs,i+1)*XRHI(irhs,i+2))/DHI(irhs,i)
!get_dbldbl_hi( div_dbldbl( sub_dbldbl(make_dbldbl(XRHI(irhs,i),XRLO(irhs,i)), add_dbldbl(mul_dbldbl(make_dbldbl(DUHI(irhs,i),DULO(irhs,i)),make_dbldbl(XRHI(irhs,i+1),XRLO(irhs,i+1))), &
                              !mul_dbldbl(make_dbldbl(DLHI(irhs,i+1),DLLO(irhs,i+1)),make_dbldbl(XRHI(irhs,i+2),XRLO(irhs,i+2))))), make_dbldbl(DHI(irhs,i),DLO(irhs,i)) ) )
tempxrlo = get_dbldbl_lo( div_dbldbl( sub_dbldbl(make_dbldbl(XRHI(irhs,i),XRLO(irhs,i)), add_dbldbl(mul_dbldbl(make_dbldbl(DUHI(irhs,i),DULO(irhs,i)),make_dbldbl(XRHI(irhs,i+1),XRLO(irhs,i+1))), &
                              mul_dbldbl(make_dbldbl(DLHI(irhs,i+1),DLLO(irhs,i+1)),make_dbldbl(XRHI(irhs,i+2),XRLO(irhs,i+2))))), make_dbldbl(DHI(irhs,i),DLO(irhs,i)) ) )
tempxihi = get_dbldbl_hi( div_dbldbl( sub_dbldbl(make_dbldbl(XIHI(irhs,i),XILO(irhs,i)), add_dbldbl(mul_dbldbl(make_dbldbl(DUHI(irhs,i),DULO(irhs,i)),make_dbldbl(XIHI(irhs,i+1),XILO(irhs,i+1))), & 
                              mul_dbldbl(make_dbldbl(DLHI(irhs,i+1),DLLO(irhs,i+1)),make_dbldbl(XIHI(irhs,i+2),XILO(irhs,i+2))))), make_dbldbl(DHI(irhs,i),DLO(irhs,i)) ) ) 
!(XIHI(irhs,i)-DUHI(irhs,i)*XIHI(irhs,i+1)-DLHI(irhs,i+1)*XIHI(irhs,i+2))/DHI(irhs,i) 
!get_dbldbl_hi( div_dbldbl( sub_dbldbl(make_dbldbl(XIHI(irhs,i),XILO(irhs,i)), add_dbldbl(mul_dbldbl(make_dbldbl(DUHI(irhs,i),DULO(irhs,i)),make_dbldbl(XIHI(irhs,i+1),XILO(irhs,i+1))), &
                              !mul_dbldbl(make_dbldbl(DLHI(irhs,i+1),DLLO(irhs,i+1)),make_dbldbl(XIHI(irhs,i+2),XILO(irhs,i+2))))), make_dbldbl(DHI(irhs,i),DLO(irhs,i)) ) ) 
tempxilo = get_dbldbl_lo( div_dbldbl( sub_dbldbl(make_dbldbl(XIHI(irhs,i),XILO(irhs,i)), add_dbldbl(mul_dbldbl(make_dbldbl(DUHI(irhs,i),DULO(irhs,i)),make_dbldbl(XIHI(irhs,i+1),XILO(irhs,i+1))), &
                              mul_dbldbl(make_dbldbl(DLHI(irhs,i+1),DLLO(irhs,i+1)),make_dbldbl(XIHI(irhs,i+2),XILO(irhs,i+2))))), make_dbldbl(DHI(irhs,i),DLO(irhs,i)) ) )
        XRHI(irhs,i) = tempxrhi
        XRLO(irhs,i) = tempxrlo
        XIHI(irhs,i) = tempxihi
        XILO(irhs,i) = tempxilo
        x (i,irhs) = DCMPLX( tempxrhi, tempxihi )

        enddo
    end if
  end subroutine epZgtsv_pressure_piv_kernel_dd


  subroutine epZgtsv_pressure_piv_dd(n, ny, nz, amphk, acphk, apphk, ak1, ak2, x, ldx, tmp, joff, koff)
    implicit none
    integer, intent(in) :: n,ny,nz,ldx,joff,koff
    real(fp_kind), dimension(:), device, intent(IN) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(n,ny*nz), device :: x
    real(fp_kind), dimension(n,ny*nz), device :: tmp
    real(fp_kind), dimension(:,:), allocatable, device :: DLHI,DLLO,DHI,DLO,DUHI,DULO,XRHI,XRLO,XIHI,XILO
    type(dim3) :: tBlock, grid

    !allocate(DLHI(n,ny*nz),DLLO(n,ny*nz),DHI(n,ny*nz),DLO(n,ny*nz),DUHI(n,ny*nz),DULO(n,ny*nz),XRHI(n,ny*nz),XRLO(n,ny*nz),XIHI(n,ny*nz),XILO(n,ny*nz))
    allocate(DLHI(ny*nz,n),DLLO(ny*nz,n),DHI(ny*nz,n),DLO(ny*nz,n),DUHI(ny*nz,n),DULO(ny*nz,n),XRHI(ny*nz,n),XRLO(ny*nz,n),XIHI(ny*nz,n),XILO(ny*nz,n))

    tBlock = dim3(128,1,1)
    grid = dim3(ceiling(real(ny*nz)/tBlock%x), 1, 1)

!    tBlock = dim3(16,16,1)
!    grid = dim3(ceiling(real(ny)/tBlock%x), ceiling(real(nz)/tBlock%y), 1)

!    tBlock = dim3(16,16,1)
!    grid = dim3(1,1,1)

    call epZgtsv_pressure_piv_kernel_dd<<<grid,tBlock>>>(n,ny*nz,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,DLHI,DLLO,DHI,DLO,DUHI,DULO,XRHI,XRLO,XIHI,XILO,joff,koff)

    deallocate(DLHI,DLLO,DHI,DLO,DUHI,DULO,XRHI,XRLO,XIHI,XILO)

  end subroutine epZgtsv_pressure_piv_dd

#if 0

  attributes(global) subroutine epZgtsv_pressure_piv_kernel(n,nrhs,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,DL,DHI,DLO,DU)
    implicit none
    integer, value, intent(in) :: n,nrhs,ny,nz,ldx
    real(fp_kind), dimension(:), device, intent(in) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(ldx,nrhs), device :: x
    real(fp_kind), dimension(ldx,nrhs), device :: DL,DU,DHI,DLO
    !type(dbldbl), dimension(ldx,nrhs), device :: D
    real(fp_kind) :: acphT_b,dui,dlip,duip,fact
    complex(fp_kind) :: temp
    integer :: irhs,i,j,k
    j = threadIdx%x + (blockIdx%x-1)*blockDim%x
    k = threadIdx%y + (blockIdx%y-1)*blockDim%y
    irhs = (k-1)*ny + j
!    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x
!    j = mod(irhs-1,ny) + 1
!    k = (irhs+ny-1)/ny
    if ( irhs <= nrhs .and. j<=ny .and. k<=nz) then
        !D(1,irhs)=make_dbldbl( 1.0d0, 0.0d0 )
        DHI(1,irhs) = 1.0d0
        DLO(1,irhs) = 0.0d0
        acphT_b = 1.0d0/(acphk(1)-ak2(j)-ak1(k))
        !if(j==1.and.k==1) print *,"orig acphT_b at 1,1,1=",acphT_b
        DU(1,irhs) = get_dbldbl_hi(mul_dbldbl(make_dbldbl(apphk(1),0.0d0),make_dbldbl(acphT_b,0.0d0)))
        x(1,irhs) = acphT_b*x(1,irhs)
        !if(j==1.and.k==1) print *,"orig x 1,1,1
        !=",REAL(x(1,irhs)),AIMAG(x(1,irhs))
        do i=1,n-1
          acphT_b = 1.0d0/(acphk(i+1)-ak2(j)-ak1(k))
          duip = apphk(i+1)*acphT_b
          dlip = amphk(i+1)*acphT_b
!          if( abs( D(i,irhs) ) >= abs( dlip ) ) then
            fact   = get_dbldbl_hi( div_dbldbl( make_dbldbl(dlip,0.0d0), make_dbldbl( DHI(i,irhs),DLO(i,irhs) ) ) )
            DHI(i+1,irhs) =  get_dbldbl_hi( make_dbldbl( 1.0d0 - fact*DU(i,irhs) , 0.0d0 ) )
            DLO(i+1,irhs) =  get_dbldbl_lo( make_dbldbl( 1.0d0 - fact*DU(i,irhs) , 0.0d0 ) )
            !if(j==1.and.k==1) print *,"ogDF",i,fact,DU(i,irhs)
            !if(j==1.and.k==1.and.i==2) print *,"ogDF2",fact,DU(i,irhs)
            x(i+1,irhs) =  x(i+1,irhs)*acphT_b - fact*x(i,irhs)
            dlip = 0.d0
         DL(i+1,irhs) = dlip
         DU(i+1,irhs) = duip
        enddo
        x(n,irhs)%re = get_dbldbl_hi( div_dbldbl( make_dbldbl(x(n,irhs)%re,0.0d0), make_dbldbl( DHI(n,irhs), DLO(n,irhs) ) ) )
        x(n,irhs)%im = get_dbldbl_hi( div_dbldbl( make_dbldbl(x(n,irhs)%im,0.0d0), make_dbldbl( DHI(n,irhs), DLO(n,irhs) ) ) )
        !if(j==1.and.k==1) print *,"orig x n,1,1
        !=",REAL(x(n,irhs)),AIMAG(x(n,irhs)),D(n,irhs)
        x(n-1,irhs)%re = get_dbldbl_hi( div_dbldbl( make_dbldbl( (x(n-1,irhs)%re-DU(n-1,irhs)*x(n,irhs)%re), 0.d0 ), make_dbldbl( DHI(n-1,irhs), DLO(n-1,irhs) ) ) )
        x(n-1,irhs)%im = get_dbldbl_hi( div_dbldbl( make_dbldbl( (x(n-1,irhs)%im-DU(n-1,irhs)*x(n,irhs)%im), 0.d0 ), make_dbldbl( DHI(n-1,irhs), DLO(n-1,irhs) ) ) )
        !if(j==1.and.k==1) print *,"orig x n-1,1,1
        !=",REAL(x(n-1,irhs)),AIMAG(x(n-1,irhs)),D(n-1,irhs),DU(n-1,irhs)
        do i=n-2,1,-1
          x(i,irhs)%re = get_dbldbl_hi( div_dbldbl( make_dbldbl( (x(i,irhs)%re-DU(i,irhs)*x(i+1,irhs)%re-DL(i+1,irhs)*x(i+2,irhs)%re), 0.d0) , make_dbldbl( DHI(i,irhs), DLO(i,irhs) ) ) )
          x(i,irhs)%im = get_dbldbl_hi( div_dbldbl( make_dbldbl( (x(i,irhs)%im-DU(i,irhs)*x(i+1,irhs)%im-DL(i+1,irhs)*x(i+2,irhs)%im), 0.d0) , make_dbldbl( DHI(i,irhs), DLO(i,irhs) ) ) )
        enddo
    end if
  end subroutine epZgtsv_pressure_piv_kernel

  subroutine epZgtsv_pressure_piv(n, ny, nz, amphk, acphk, apphk, ak1, ak2, x, ldx, tmp)
    implicit none
    integer, intent(in) :: n,ny,nz,ldx
    real(fp_kind), dimension(:), device, intent(IN) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(n,ny*nz), device :: x
    real(fp_kind), dimension(n,ny*nz), device :: tmp
    real(fp_kind), dimension(:,:), allocatable, device :: DL,DU,DHI,DLO
    !type(dbldbl), dimension(:,:), allocatable, device :: D
    type(dim3) :: tBlock, grid

     allocate(DL(n,ny*nz),DHI(n,ny*nz),DLO(n,ny*nz),DU(n,ny*nz))

!    tBlock = dim3(128,1,1)
!    grid = dim3(ceiling(real(ny*nz)/tBlock%x), 1, 1)

    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(ny)/tBlock%x), ceiling(real(nz)/tBlock%y), 1)

!    tBlock = dim3(16,16,1)
!    grid = dim3(1,1,1)

    call epZgtsv_pressure_piv_kernel<<<grid,tBlock>>>(n,ny*nz,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,DL,DHI,DLO,DU)

    deallocate(DL,DHI,DLO,DU)

  end subroutine epZgtsv_pressure_piv

#endif

#if 1

  attributes(global) subroutine epZgtsv_pressure_piv_kernel(n,nrhs,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,DL,D,DU,joff,koff)
    implicit none
    integer, value, intent(in) :: n,nrhs,ny,nz,ldx,joff,koff
    real(fp_kind), dimension(:), device, intent(in) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(ldx,nrhs), device :: x
    real(fp_kind), dimension(ldx,nrhs), device :: DL,D,DU
    real(fp_kind) :: acphT_b,dui,dlip,duip,fact
    complex(fp_kind) :: temp
    integer :: irhs,i,j,k
    j = threadIdx%x + (blockIdx%x-1)*blockDim%x + joff - 1
    k = threadIdx%y + (blockIdx%y-1)*blockDim%y + koff - 1
    irhs = (k-1)*ny + j
!    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x
!    j = mod(irhs-1,ny) + 1
!    k = (irhs+ny-1)/ny
    if ( irhs <= nrhs .and. j<=ny .and. k<=nz) then
        D(1,irhs)=1.0d0
        acphT_b = 1.0d0/(acphk(1)-ak2(j)-ak1(k))
        !if(j==1.and.k==1) print *,"orig acphT_b at 1,1,1=",acphT_b
        DU(1,irhs) = apphk(1)*acphT_b
        x(1,irhs) = acphT_b*x(1,irhs)
        !if(j==1.and.k==1) print *,"orig x 1,1,1 =",REAL(x(1,irhs)),AIMAG(x(1,irhs))
        do i=1,n-1
          acphT_b = 1.0d0/(acphk(i+1)-ak2(j)-ak1(k))
          duip = apphk(i+1)*acphT_b
          dlip = amphk(i+1)*acphT_b
          if( abs( D(i,irhs) ) >= abs( dlip ) ) then
            fact   = dlip / D(i,irhs)
            D(i+1,irhs) =  1.0d0 - fact*DU(i,irhs)
!            if(j==1.and.k==1) print *,"ogDF",i,fact,DU(i,irhs)
            !if(j==1.and.k==1.and.i==2) print *,"ogDF2",fact,DU(i,irhs)
            x(i+1,irhs) =  x(i+1,irhs)*acphT_b - fact*x(i,irhs)
            dlip = 0.d0
#if 1
          else
!            print *,"pivot in row: ",i
            fact = D(i,irhs) / dlip
            D(i,irhs) = dlip
            D(i+1,irhs) = DU(i,irhs) - fact
            dlip = duip
            duip = -fact*duip
            DU(i,irhs) = 1.0d0
            temp = x(i,irhs)
            x(i,irhs) = x(i+1,irhs)*acphT_b
            x(i+1,irhs) = temp - fact*x(i+1,irhs)*acphT_b
         endif
#endif
         DL(i+1,irhs) = dlip
         DU(i+1,irhs) = duip
        enddo
        x(n,irhs) = x(n,irhs)/D(n,irhs)
        !if(j==1.and.k==1) print *,"orig x n,1,1 =",REAL(x(n,irhs)),AIMAG(x(n,irhs)),D(n,irhs)
        x(n-1,irhs) = (x(n-1,irhs)-DU(n-1,irhs)*x(n,irhs))/D(n-1,irhs)
        !if(j==1.and.k==1) print *,"orig x n-1,1,1 =",REAL(x(n-1,irhs)),AIMAG(x(n-1,irhs)),D(n-1,irhs),DU(n-1,irhs)
        do i=n-2,1,-1
          x(i,irhs) = (x(i,irhs)-DU(i,irhs)*x(i+1,irhs)-DL(i+1,irhs)*x(i+2,irhs))/D(i,irhs)
        enddo
    end if
  end subroutine epZgtsv_pressure_piv_kernel

  subroutine epZgtsv_pressure_piv(n, ny, nz, amphk, acphk, apphk, ak1, ak2, x, ldx, tmp, joff, koff)
    implicit none
    integer, intent(in) :: n,ny,nz,ldx,joff,koff
    real(fp_kind), dimension(:), device, intent(IN) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(n,ny*nz), device :: x
    real(fp_kind), dimension(n,ny*nz), device :: tmp
    real(fp_kind), dimension(:,:), allocatable, device :: DL,D,DU
    type(dim3) :: tBlock, grid

     allocate(DL(n,ny*nz),D(n,ny*nz),DU(n,ny*nz))

!    tBlock = dim3(128,1,1)
!    grid = dim3(ceiling(real(ny*nz)/tBlock%x), 1, 1)

    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(ny)/tBlock%x), ceiling(real(nz)/tBlock%y), 1)

!    tBlock = dim3(16,16,1)
!    grid = dim3(1,1,1)

    call epZgtsv_pressure_piv_kernel<<<grid,tBlock>>>(n,ny*nz,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,DL,D,DU,joff,koff)

    deallocate(DL,D,DU)

  end subroutine epZgtsv_pressure_piv

#endif

  attributes(global) subroutine epZgtsv_pressure_kernel(n,nrhs,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ny,nz,ldx
    real(fp_kind), dimension(:), device, intent(in) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(ldx,nrhs), device :: x
    real(fp_kind), dimension(ldx,nrhs), device :: tmp

    real(fp_kind) :: m,acphT_b,a,c
    integer :: irhs,i,j,k

    j = threadIdx%x + (blockIdx%x-1)*blockDim%x
    k = threadIdx%y + (blockIdx%y-1)*blockDim%y
    irhs = (k-1)*ny + j
!    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x
!    j = mod(irhs-1,ny) + 1
!    k = (irhs+ny-1)/ny

    if ( irhs <= nrhs .and. j<=ny .and. k<=nz) then

    acphT_b = 1.0d0/(acphk(1)-ak2(j)-ak1(k))
    tmp(1,irhs) = apphk(1)*acphT_b
      x(1,irhs) = x(1,irhs)*acphT_b

    do i=2,N
      acphT_b = 1.0d0/(acphk(i)-ak2(j)-ak1(k))
      a = amphk(i)*acphT_b
      c = apphk(i)*acphT_b
      m = 1.0d0 - tmp(i-1,irhs)*a
      tmp(i,irhs) = c/m
        x(i,irhs) = (x(i,irhs)*acphT_b - x(i-1,irhs)*a)/m
    end do

    do i=n-1,1,-1
      x(i,irhs) = x(i,irhs) - tmp(i,irhs)*x(i+1,irhs)
!      if(i==1.and.j==1.and.k==1) then
!        print *,"x(1,1,1)=",x(i,irhs)
!      endif
    end do

    end if

  end subroutine epZgtsv_pressure_kernel

  subroutine epZgtsv_pressure(n, ny, nz, amphk, acphk, apphk, ak1, ak2, x, ldx, tmp)
    implicit none
    integer, intent(in) :: n,ny,nz,ldx
    real(fp_kind), dimension(:), device, intent(IN) :: amphk,acphk,apphk,ak1,ak2
    complex(fp_kind), dimension(n,ny*nz), device :: x
    real(fp_kind), dimension(n,ny*nz), device :: tmp
    type(dim3) :: tBlock, grid

!    tBlock = dim3(128,1,1)
!    grid = dim3(ceiling(real(ny*nz)/tBlock%x), 1, 1)

    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(ny)/tBlock%x), ceiling(real(nz)/tBlock%y), 1)

!    tBlock = dim3(16,16,1)
!    grid = dim3(1,1,1)

    call epZgtsv_pressure_kernel<<<grid,tBlock>>>(n,ny*nz,ny,nz,amphk,acphk,apphk,ak1,ak2,x,ldx,tmp)

  end subroutine epZgtsv_pressure

end module ep_solve
#endif

#ifdef USE_CUDA
module ep_trans
  use precision
  use cudafor
contains
  attributes(global) subroutine trans_xy_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    real(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    real(fp_kind), dimension(n2,n1,n3), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n2 .and. j<=n1) then
      odata(i,j,k) = idata(j,i,k)
    end if
  end subroutine trans_xy_kernel
  subroutine trans_xy(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    real(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    real(fp_kind), dimension(n2,n1,n3), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n2)/tBlock%x), ceiling(real(n1)/tBlock%y), n3)
    call trans_xy_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_xy
  attributes(global) subroutine trans_yx_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    complex(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n2,n1,n3), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n2 .and. j<=n1) then
      odata(i,j,k) = idata(j,i,k)
    end if
  end subroutine trans_yx_kernel
  subroutine trans_yx(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    complex(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n2,n1,n3), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n2)/tBlock%x), ceiling(real(n1)/tBlock%y), n3)
    call trans_yx_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_yx
  attributes(global) subroutine trans_yz_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    complex(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n3,n1,n2), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n3 .and. j<=n1) then
      odata(i,j,k) = idata(j,k,i)
    end if
  end subroutine trans_yz_kernel
  subroutine trans_yz(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    complex(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n3,n1,n2), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n3)/tBlock%x), ceiling(real(n1)/tBlock%y), n2)
    call trans_yz_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_yz
  attributes(global) subroutine trans_zy_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    complex(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n2,n3,n1), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n2 .and. j<=n3) then
      odata(i,j,k) = idata(k,i,j)
    end if
  end subroutine trans_zy_kernel
  subroutine trans_zy(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    complex(fp_kind), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(fp_kind), dimension(n2,n3,n1), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n2)/tBlock%x), ceiling(real(n3)/tBlock%y), n1)
    call trans_zy_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_zy

end module ep_trans
#endif


