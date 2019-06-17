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

! This is the main 2D pencil decomposition module

module decomp_2d
  use nvtx
#ifdef USE_CUDA
  use cudafor
  use mpiDeviceUtil 
#endif
  implicit none

  include "mpif.h"

#ifdef GLOBAL_ARRAYS
#include "mafdecls.fh"
#include "global.fh"
#endif

  private        ! Make everything private unless declared public

  integer, parameter, public :: mytype = KIND(0.0D0)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef GLOBAL_ARRAYS
  integer, parameter, public :: ga_real_type = MT_F_DBL
  integer, parameter, public :: ga_complex_type = MT_F_DCPL
#endif

  integer, save, public :: mytype_bytes

  ! some key global variables
  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: nrank  ! local MPI rank 
  integer, save, public :: nproc  ! total number of processors

  ! parameters for 2D Cartesian topology 
  integer, save, dimension(2) :: dims, coord
  logical, save, dimension(2) :: periodic
  integer, save, public :: DECOMP_2D_COMM_CART_X, &
       DECOMP_2D_COMM_CART_Y, DECOMP_2D_COMM_CART_Z 
  integer, save :: DECOMP_2D_COMM_ROW, DECOMP_2D_COMM_COL
  integer, save :: row_rank, col_rank
  ! define neighboring blocks (to be used in halo-cell support)
  !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, save, dimension(1,6) :: neighbour 

  ! flags for periodic condition in three dimensions
  logical, save :: periodic_x, periodic_y, periodic_z

#ifdef USE_CUDA
    type(cudaEvent), allocatable, dimension(:) :: a2a_event
    integer(kind=cuda_stream_kind), public, save :: a2a_comp, a2a_h2d, a2a_d2h
    integer, allocatable, dimension(:) :: a2a_requests
    integer, allocatable, dimension(:,:) :: a2a_status
#endif

#ifdef SHM
  ! derived type to store shared-memory info
  TYPE, public :: SMP_INFO
     integer MPI_COMM          ! SMP associated with this communicator
     integer NODE_ME           ! rank in this communicator
     integer NCPU              ! size of this communicator
     integer SMP_COMM          ! communicator for SMP-node masters
     integer CORE_COMM         ! communicator for cores on SMP-node
     integer SMP_ME            ! SMP-node id starting from 1 ... NSMP
     integer NSMP              ! number of SMP-nodes in this communicator
     integer CORE_ME           ! core id starting from 1 ... NCORE
     integer NCORE             ! number of cores on this SMP-node
     integer MAXCORE           ! maximum no. cores on any SMP-node
     integer N_SND             ! size of SMP shared memory buffer
     integer N_RCV             ! size of SMP shared memory buffer
     integer(8) SND_P          ! SNDBUF address (cray pointer), for real 
     integer(8) RCV_P          ! RCVBUF address (cray pointer), for real
     integer(8) SND_P_c        ! for complex
     integer(8) RCV_P_c        ! for complex
  END TYPE SMP_INFO
#endif

  ! derived type to store decomposition info for a given global data size
  TYPE, public :: DECOMP_INFO
     ! staring/ending index and size of data held by current processor
     integer, dimension(3) :: xst, xen, xsz  ! x-pencil
     integer, dimension(3) :: yst, yen, ysz  ! y-pencil
     integer, dimension(3) :: zst, zen, zsz  ! z-pencil
     integer, dimension(3) :: wst, wen, wsz  ! x-pencil

     ! in addition to local information, processors also need to know 
     ! some global information for global communications to work 

     ! how each dimension is distributed along pencils
     integer, allocatable, dimension(:) :: &
          x1dist, y1dist, y2dist, z2dist, x2dist, z1dist, &
          x1st, y1st, y2st, z2st, x2st, z1st, &
          x1en, y1en, y2en, z2en, x2en, z1en

     ! send/receive buffer counts and displacements for MPI_ALLTOALLV
     integer, allocatable, dimension(:) :: &
          x1cnts, y1cnts, y2cnts, z2cnts, x2cnts, z1cnts, w1cnts, w2cnts
     integer, allocatable, dimension(:) :: &
          x1disp, y1disp, y2disp, z2disp, x2disp, z1disp, w1disp, w2disp

     integer, allocatable, dimension(:) :: &
          x1idx, y1idx, y2idx, z2idx, x2idx, z1idx, w1idx, w2idx

     ! buffer counts for MPI_ALLTOALL: either for evenly distributed data
     ! or for padded-alltoall
     integer :: x1count, y1count, y2count, z2count, x2count, z1count

     ! buffer counts, displacements and types for MPI_Alltoallw to transform
     ! directly between x- and z-pencils
     ! This is only for the complex datatype
     integer,dimension(:),allocatable::zcnts_xz,ztypes_xz
     integer,dimension(:),allocatable::xcnts_xz,xtypes_xz
     integer,dimension(:),allocatable::zdispls_xz,xdispls_xz

     ! evenly distributed data
     logical :: even

#ifdef SHM
     ! For shared-memory implementation

     ! one instance of this derived type for each communicator
     ! shared moemory info, such as which MPI rank belongs to which node
     TYPE(SMP_INFO) :: ROW_INFO, COL_INFO

     ! shared send/recv buffers for ALLTOALLV
     integer, allocatable, dimension(:) :: x1cnts_s, y1cnts_s, &
          y2cnts_s, z2cnts_s
     integer, allocatable, dimension(:) :: x1disp_s, y1disp_s, &
          y2disp_s, z2disp_s
     ! A copy of original buffer displacement (will be overwriten)
     integer, allocatable, dimension(:) :: x1disp_o, y1disp_o, &
          y2disp_o, z2disp_o
#endif
  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), save :: decomp_main
  TYPE(DECOMP_INFO), save :: decomp_ph,decomp_sp

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

  ! Additional start and end variables for CPU/GPU hybrid operation
  ! Set to full pencil dimensions in non hybrid cases
  integer, save, dimension(3), public :: xstart_cpu, xend_cpu, xsize_cpu
  integer, save, dimension(3), public :: xstart_gpu, xend_gpu, xsize_gpu

  ! These are the buffers used by MPI_ALLTOALL(V) calls
  integer, save :: decomp_buf_size = 0
#ifdef USE_CUDA
  real(mytype),    allocatable, dimension(:), public, pinned :: work1_r, work2_r
  complex(mytype), allocatable, dimension(:), public, pinned :: work1_c, work2_c
#else
  real(mytype),    allocatable, dimension(:), public :: work1_r, work2_r
  complex(mytype), allocatable, dimension(:) :: work1_c, work2_c
#endif

#ifdef USE_CUDA
  type(C_devptr) :: cptr1,cptr2
  real(mytype),    device, allocatable, dimension(:), public :: work1_r_d, work2_r_d
  complex(mytype), device, target, allocatable, dimension(:), public :: work1_c_d, work2_c_d
  integer,save :: alloc_work_r, alloc_work_c, first_time_alloc_work_r, first_time_alloc_work_c
#endif

  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       transpose_z_to_x, transpose_x_to_z, &
       decomp_info_init, decomp_info_finalize, partition, &
#ifdef GLOBAL_ARRAYS
       get_global_array, &
#endif
       alloc_x, alloc_y, alloc_z, &
       update_halo, decomp_2d_abort, &
       get_decomp_info



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are routines to perform global data transpositions
  ! 
  !   Four combinations are available, enough to cover all situations
  !    - transpose_x_to_y (X-pencil --> Y-pencil)
  !    - transpose_y_to_z (Y-pencil --> Z-pencil)
  !    - transpose_z_to_y (Z-pencil --> Y-pencil)
  !    - transpose_y_to_x (Y-pencil --> X-pencil)
  !
  !   Generic interface provided here to support multiple data types
  !    - real and complex types supported through generic interface
  !    - single/double precision supported through pre-processing
  !       * see 'mytype' variable at the beginning
  !    - an optional argument can be supplied to transpose data whose 
  !      global size is not the default nx*ny*nz 
  !       * as the case in fft r2c/c2r interface 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface transpose_x_to_y
     module procedure transpose_x_to_y_real
     module procedure transpose_x_to_y_complex
#ifdef USE_CUDA
     module procedure transpose_x_to_y_real_d
     module procedure transpose_x_to_y_complex_d
#endif
  end interface transpose_x_to_y
  
  interface transpose_y_to_z
     module procedure transpose_y_to_z_real
     module procedure transpose_y_to_z_complex
#ifdef USE_CUDA
     module procedure transpose_y_to_z_complex_d
#endif
  end interface transpose_y_to_z
  
  interface transpose_z_to_y
     module procedure transpose_z_to_y_real
     module procedure transpose_z_to_y_complex
#ifdef USE_CUDA
     module procedure transpose_z_to_y_complex_d
#endif
  end interface transpose_z_to_y

  interface transpose_y_to_x
     module procedure transpose_y_to_x_real
     module procedure transpose_y_to_x_complex
#ifdef USE_CUDA
     module procedure transpose_y_to_x_real_d
     module procedure transpose_y_to_x_complex_d
#endif
  end interface transpose_y_to_x

  interface transpose_z_to_x
! not available
!     module procedure transpose_z_to_x_real
     module procedure transpose_z_to_x_complex
#ifdef USE_CUDA
     module procedure transpose_z_to_x_complex_d
#endif

  end interface transpose_z_to_x

  interface transpose_x_to_z
! not available
!     module procedure transpose_x_to_z_real
     module procedure transpose_x_to_z_complex
#ifdef USE_CUDA
     module procedure transpose_x_to_z_complex_d
#endif

  end interface transpose_x_to_z

  interface update_halo
         module procedure update_halo_real_cpu
#ifdef USE_CUDA
         module procedure update_halo_real_gpu
#endif
  end interface update_halo

  interface alloc_x
     module procedure alloc_x_real
     module procedure alloc_x_complex
  end interface alloc_x

  interface alloc_y
     module procedure alloc_y_real
     module procedure alloc_y_complex
  end interface alloc_y

  interface alloc_z
     module procedure alloc_z_real
     module procedure alloc_z_complex
  end interface alloc_z

contains

#ifdef SHM_DEBUG
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For debugging, print the shared-memory structure
  subroutine print_smp_info(s)
    TYPE(SMP_INFO) :: s
    write(10,*) 'size of current communicator:', s%NCPU
    write(10,*) 'rank in current communicator:', s%NODE_ME
    write(10,*) 'number of SMP-nodes in this communicator:', s%NSMP
    write(10,*) 'SMP-node id (1 ~ NSMP):', s%SMP_ME
    write(10,*) 'NCORE - number of cores on this SMP-node', s%NCORE
    write(10,*) 'core id (1 ~ NCORE):', s%CORE_ME
    write(10,*) 'maximum no. cores on any SMP-node:', s%MAXCORE
    write(10,*) 'size of SMP shared memory SND buffer:', s%N_SND
    write(10,*) 'size of SMP shared memory RCV buffer:', s%N_RCV
  end subroutine print_smp_info
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc)
#ifdef USE_HYBRID
    use param, only: gc_split_ratio
#endif

    implicit none

    integer, intent(IN) :: nx,ny,nz,p_row,p_col
    logical, dimension(3), intent(IN), optional :: periodic_bc
    
    integer :: errorcode, ierror, row, col, deviceID, maxrc, m, istat
    
#ifdef SHM_DEBUG
    character(len=80) fname
#endif

    nx_global = nx
    ny_global = ny
    nz_global = nz

    if (present(periodic_bc)) then
       periodic_x = periodic_bc(1)
       periodic_y = periodic_bc(2)
       periodic_z = periodic_bc(3)
    else
       periodic_x = .false.
       periodic_y = .false.
       periodic_z = .false.
    end if

    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

#ifdef USE_CUDA
     first_time_alloc_work_r = 1
     first_time_alloc_work_c = 1
     alloc_work_r = 1
     alloc_work_c = 1
     call assignDevice(deviceID)
!     print *,"Rank: :",nrank," Assigning GPU :",deviceID

    !maxrc = row
    !if(col>row) maxrc = col
    maxrc = nproc

    allocate( a2a_event(0:maxrc*2) )
    allocate( a2a_requests(0:maxrc*2) )
    allocate( a2a_status(MPI_STATUS_SIZE,0:maxrc*2) )
    do m=0,maxrc*2
      istat = cudaEventCreate( a2a_event(m) )
    enddo
    istat = cudaStreamCreate( a2a_comp )
    istat = cudaStreamCreate( a2a_d2h  )
    istat = cudaStreamCreate( a2a_h2d  )

#endif

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, row, col)
    else
       if (nproc /= p_row*p_col) then
          errorcode = 1
          call decomp_2d_abort(errorcode, &
               'Invalid 2D processor grid - nproc /= p_row*p_col')
       else
          row = p_row
          col = p_col
          if (nrank==0) then
             write(*,*) 'the used processor grid is ', row,' by ', col
          end if

       end if
    end if
    
    ! Create 2D Catersian topology
    ! Note that in order to support periodic B.C. in the halo-cell code,
    ! need to create multiple topology objects: DECOMP_2D_COMM_CART_?,
    ! corresponding to three pencil orientations. They contain almost
    ! identical topological information but allow different combinations
    ! of periodic conditions.
    dims(1) = row
    dims(2) = col
    periodic(1) = periodic_y
    periodic(2) = periodic_z
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., &  ! do not reorder rank
         DECOMP_2D_COMM_CART_X, ierror)
    periodic(1) = periodic_x
    periodic(2) = periodic_z
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Y, ierror)
    periodic(1) = periodic_x
    periodic(2) = periodic_y
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Z, ierror)

    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    
    ! derive communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
         DECOMP_2D_COMM_COL,ierror)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
         DECOMP_2D_COMM_ROW,ierror)

    call MPI_COMM_RANK(DECOMP_2D_COMM_COL,col_rank,ierror)
    call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,row_rank,ierror)

    !print *,"rank: ",nrank,"col: ",col_rank,"row: ",row_rank


    ! gather information for halo-cell support code
    call init_neighbour
    
    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,decomp_main)
    
    ! make a copy of the decomposition information associated with the
    ! default global size in these global variables so applications can
    ! use them to create data structures 
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz

    ! If using GPU-CPU hybrid version, compute subpencil splits
    xstart_gpu = xstart; xend_gpu = xend; xsize_gpu = xsize
    xstart_cpu = xstart; xend_cpu = xend; xsize_cpu = xsize

#ifdef USE_HYBRID
    xend_gpu(3) = xstart(3) + floor(gc_split_ratio * xsize(3)) - 1
    xstart_cpu(3) = xend_gpu(3) + 1
    xend_cpu(3) = xend(3)
    xsize_gpu(3) = xend_gpu(3) - xstart_gpu(3) + 1
    xsize_cpu(3) = xend_cpu(3) - xstart_cpu(3) + 1
#endif
    

    !decomp_ph=decomp_main

#ifdef SHM_DEBUG
    ! print out shared-memory information
    write(fname,99) nrank
99  format('log',I2.2)
    open(10,file=fname)
    write(10,*)'I am mpi rank ', nrank, 'Total ranks ', nproc
    write(10,*)' '
    write(10,*)'Global data size:'
    write(10,*)'nx*ny*nz', nx,ny,nz
    write(10,*)' '
    write(10,*)'2D processor grid:'
    write(10,*)'p_row*p_col:', dims(1), dims(2)
    write(10,*)' '
    write(10,*)'Portion of global data held locally:'
    write(10,*)'xsize:',xsize
    write(10,*)'ysize:',ysize
    write(10,*)'zsize:',zsize
    write(10,*)' '
    write(10,*)'How pensils are to be divided and sent in alltoallv:'
    write(10,*)'x1dist:',decomp_main%x1dist
    write(10,*)'y1dist:',decomp_main%y1dist
    write(10,*)'y2dist:',decomp_main%y2dist
    write(10,*)'z2dist:',decomp_main%z2dist
    write(10,*)' '
    write(10,*)'######Shared buffer set up after this point######'
    write(10,*)' '
    write(10,*) 'col communicator detais:'
    call print_smp_info(decomp_main%COL_INFO)
    write(10,*)' '
    write(10,*) 'row communicator detais:'
    call print_smp_info(decomp_main%ROW_INFO)
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of per-core buffers'
    write(10,*)'x1cnts:',decomp_main%x1cnts
    write(10,*)'y1cnts:',decomp_main%y1cnts
    write(10,*)'y2cnts:',decomp_main%y2cnts
    write(10,*)'z2cnts:',decomp_main%z2cnts
    write(10,*)'x1disp:',decomp_main%x1disp
    write(10,*)'y1disp:',decomp_main%y1disp
    write(10,*)'y2disp:',decomp_main%y2disp
    write(10,*)'z2disp:',decomp_main%z2disp
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of shared buffers'
    write(10,*)'x1cnts:',decomp_main%x1cnts_s
    write(10,*)'y1cnts:',decomp_main%y1cnts_s
    write(10,*)'y2cnts:',decomp_main%y2cnts_s
    write(10,*)'z2cnts:',decomp_main%z2cnts_s
    write(10,*)'x1disp:',decomp_main%x1disp_s
    write(10,*)'y1disp:',decomp_main%y1disp_s
    write(10,*)'y2disp:',decomp_main%y2disp_s
    write(10,*)'z2disp:',decomp_main%z2disp_s
    write(10,*)' '
    close(10)
#endif

    ! determine the number of bytes per float number
    ! do not use 'mytype' which is compiler dependent
    ! also possible to use inquire(iolength=...) 
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)

#ifdef EVEN
    if (nrank==0) write(*,*) 'Padded ALLTOALL optimisation on'
#endif 

    return
  end subroutine decomp_2d_init
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to clean things up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_finalize
    
    implicit none
    integer ierr

    call decomp_info_finalize(decomp_main)

    decomp_buf_size = 0
    deallocate(work1_r, work2_r, work1_c, work2_c)

    call MPI_Finalize(ierr)
    
    return
  end subroutine decomp_2d_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the default decomposition object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_decomp_info(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(OUT) :: decomp

    !decomp = decomp_main

    return
  end subroutine get_decomp_info
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init(nx,ny,nz,decomp)
    use iso_c_binding
    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: buf_size, status, errorcode, i

    ! verify the global size can actually be distributed as pencils
    if (nx<dims(1) .or. ny<dims(1) .or. ny<dims(2) .or. nz<dims(2)) then
       errorcode = 6
       call decomp_2d_abort(errorcode, &
            'Invalid 2D processor grid. ' // &
            'Make sure that min(nx,ny) >= p_row and ' // &
            'min(ny,nz) >= p_col')
    end if
    
    if (mod(nx,dims(1))==0 .and. mod(ny,dims(1))==0 .and. &
         mod(ny,dims(2))==0 .and. mod(nz,dims(2))==0) then
       decomp%even = .true.
    else
       decomp%even = .false.
    end if

    ! distribute mesh points
    allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
             decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1), &
             decomp%x2dist(0:dims(2)-1),decomp%z1dist(0:dims(1)-1))

    allocate(decomp%x1st(0:dims(1)-1),decomp%x1en(0:dims(1)-1), &
             decomp%y1st(0:dims(1)-1),decomp%y1en(0:dims(1)-1), &
             decomp%z1st(0:dims(1)-1),decomp%z1en(0:dims(1)-1))

    allocate(decomp%x2st(0:dims(2)-1),decomp%x2en(0:dims(2)-1), &
             decomp%y2st(0:dims(2)-1),decomp%y2en(0:dims(2)-1), &
             decomp%z2st(0:dims(2)-1),decomp%z2en(0:dims(2)-1))

    call get_dist(nx,ny,nz,decomp)

    allocate(decomp%x1idx(0:dims(1)-1),decomp%y1idx(0:dims(1)-1), &
             decomp%y2idx(0:dims(2)-1),decomp%z2idx(0:dims(2)-1), &
             decomp%x2idx(0:dims(2)-1),decomp%z1idx(0:dims(1)-1))

    decomp%x1idx(0) = 1
    decomp%y1idx(0) = 1
    decomp%z1idx(0) = 1
    do i=1,dims(1)-1
      decomp%x1idx(i) = decomp%x1idx(i-1) + decomp%x1dist(i-1)
      decomp%y1idx(i) = decomp%y1idx(i-1) + decomp%y1dist(i-1)
      decomp%z1idx(i) = decomp%z1idx(i-1) + decomp%z1dist(i-1)
    end do

    decomp%x2idx(0) = 1
    decomp%y2idx(0) = 1
    decomp%z2idx(0) = 1
    do i=1,dims(2)-1
      decomp%x2idx(i) = decomp%x2idx(i-1) + decomp%x2dist(i-1)
      decomp%y2idx(i) = decomp%y2idx(i-1) + decomp%y2dist(i-1)
      decomp%z2idx(i) = decomp%z2idx(i-1) + decomp%z2dist(i-1)
    end do

    ! generate partition information - starting/ending index etc.
    call partition(nx, ny, nz, (/ 1,2,3 /), &
         decomp%xst, decomp%xen, decomp%xsz)
    call partition(nx, ny, nz, (/ 2,1,3 /), &
         decomp%yst, decomp%yen, decomp%ysz)
    call partition(nx, ny, nz, (/ 2,3,1 /), &
         decomp%zst, decomp%zen, decomp%zsz)
    call partition(nx, ny, nz, (/ 1,3,2 /), &
         decomp%wst, decomp%wen, decomp%wsz)

    ! prepare send/receive buffer displacement and count for ALLTOALL(V)
    allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
             decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1), &
             decomp%z1cnts(0:dims(1)-1),decomp%x2cnts(0:dims(2)-1), &
             decomp%w1cnts(0:dims(1)-1),decomp%w2cnts(0:dims(2)-1))

    allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
             decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1), &
             decomp%x2disp(0:dims(2)-1),decomp%z1disp(0:dims(1)-1), &
             decomp%w2disp(0:dims(2)-1),decomp%w1disp(0:dims(1)-1))

#ifndef USE_CUDA
    ! allocate arrays for MPI_ALLTOALLW calls
    allocate(decomp%xcnts_xz(nproc),decomp%zcnts_xz(nproc))
    allocate(decomp%xtypes_xz(nproc),decomp%ztypes_xz(nproc))
    allocate(decomp%xdispls_xz(nproc))
    allocate(decomp%zdispls_xz(nproc))
#endif
    call prepare_buffer(decomp)

    ! allocate memory for the MPI_ALLTOALL(V) buffers
    ! define the buffers globally for performance reason
    
    buf_size = max(decomp%wsz(1)*(decomp%wsz(2)+2)*decomp%wsz(3), &
               max(decomp%xsz(1)*(decomp%xsz(2)+2)*decomp%xsz(3), &
               max(decomp%ysz(1)*(decomp%ysz(2)+2)*decomp%ysz(3), &
                   decomp%zsz(1)*(decomp%zsz(2)+2)*decomp%zsz(3))) )


    ! check if additional memory is required
    ! *** TODO: consider how to share the real/complex buffers 
#ifdef USE_CUDA
    if (buf_size > decomp_buf_size) then
       decomp_buf_size = buf_size
    end if
    if (alloc_work_r+alloc_work_c>0 .and. first_time_alloc_work_r+first_time_alloc_work_c>0) then
      buf_size = max(buf_size,decomp_buf_size)
#else
    if (buf_size > decomp_buf_size) then
       decomp_buf_size = buf_size
#endif

#ifdef USE_CUDA
    if(alloc_work_r == 1) then
      !  print *,"allocating r: ",buf_size
        first_time_alloc_work_r = 0
#endif
       if (allocated(work1_r)) deallocate(work1_r)
       if (allocated(work2_r)) deallocate(work2_r)
#ifdef USE_CUDA
    end if
    if(alloc_work_c == 1) then
       !print *,"allocating c: ",buf_size
       first_time_alloc_work_c = 0   
#endif
       if (allocated(work1_c)) deallocate(work1_c)
       if (allocated(work2_c)) deallocate(work2_c)
#ifdef USE_CUDA
    end if
    if(alloc_work_r == 1) then
#endif
       allocate(work1_r(buf_size), STAT=status)
       allocate(work2_r(buf_size), STAT=status)
#ifdef USE_CUDA
    end if
    if(alloc_work_c == 1) then
#endif
       allocate(work1_c(buf_size), STAT=status)
       allocate(work2_c(buf_size), STAT=status)
#ifdef USE_CUDA
    end if
    if(alloc_work_r == 1) then
       if (allocated(work1_r_d)) deallocate(work1_r_d)
       if (allocated(work2_r_d)) deallocate(work2_r_d)
    end if
    if(alloc_work_c == 1) then
       if (allocated(work1_c_d)) deallocate(work1_c_d)
       !if (allocated(work2_c_d)) deallocate(work2_c_d)
    end if
    if(alloc_work_c == 1) then
       allocate(work1_c_d(buf_size), STAT=status)
       !allocate(work2_c_d(buf_size), STAT=status)
       work1_c_d = (0.d0,0.d0)
    end if
    if(alloc_work_r == 1) then
       cptr1 = C_devloc(work1_c_d(buf_size/2 + 1))
       call C_F_pointer(cptr1,work1_r_d,[buf_size])
       cptr2 = C_devloc(work1_c_d)
       call C_F_pointer(cptr2,work2_r_d,[buf_size])
       !allocate(work1_r_d(buf_size), STAT=status)
       !allocate(work2_r_d(buf_size), STAT=status)
       work1_r_d = 0.d0
       work2_r_d = 0.d0
    end if
#endif
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
    end if

    return
  end subroutine decomp_info_init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory associated with a DECOMP_INFO object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_finalize(decomp)

    implicit none

    integer :: i
    integer :: ierror
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    deallocate(decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist)
    deallocate(decomp%z1dist,decomp%x2dist)
    deallocate(decomp%x1st,decomp%y1st,decomp%y2st,decomp%z2st)
    deallocate(decomp%z1st,decomp%x2st)
    deallocate(decomp%x1en,decomp%y1en,decomp%y2en,decomp%z2en)
    deallocate(decomp%z1en,decomp%x2en)
    deallocate(decomp%x1cnts,decomp%y1cnts,decomp%y2cnts,decomp%z2cnts)
    deallocate(decomp%z1cnts,decomp%x2cnts,decomp%w1cnts,decomp%w2cnts)
    deallocate(decomp%x1disp,decomp%y1disp,decomp%y2disp,decomp%z2disp)
    deallocate(decomp%z1disp,decomp%x2disp,decomp%w1disp,decomp%w2disp)
    deallocate(decomp%x1idx,decomp%y1idx,decomp%z1idx)
    deallocate(decomp%x2idx,decomp%y2idx,decomp%z2idx)

#ifndef USE_CUDA
    do i=1,nproc
      if (decomp%ztypes_xz(i).ne.MPI_DATATYPE_NULL) then
        call MPI_Type_free(decomp%ztypes_xz(i),ierror)
      endif
      if (decomp%xtypes_xz(i).ne.MPI_DATATYPE_NULL) then
        call MPI_Type_free(decomp%xtypes_xz(i),ierror)
      endif
    enddo    
    deallocate(decomp%xcnts_xz,decomp%zcnts_xz)
    deallocate(decomp%xtypes_xz,decomp%ztypes_xz)
    deallocate(decomp%xdispls_xz)
    deallocate(decomp%zdispls_xz)
#endif

    return
  end subroutine decomp_info_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find sub-domain information held by current processor
  !   INPUT: 
  !     nx, ny, nz - global data dimension
  !     pdim(3)    - number of processor grid in each dimension, 
  !                  valid values: 1 - distibute locally; 
  !                                2 - distribute across p_row; 
  !                                3 - distribute across p_col
  !   OUTPUT:
  !     lstart(3)  - starting index
  !     lend(3)    - ending index
  !     lsize(3)   - size of the sub-block (redundant) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    integer, dimension(3), intent(IN) :: pdim
    integer, dimension(3), intent(OUT) :: lstart, lend, lsize

    integer, allocatable, dimension(:) :: st,en,sz
    integer :: i, gsize

    do i = 1, 3
 
      if (i==1) then
        gsize = nx
      else if (i==2) then
        gsize = ny
      else if (i==3) then
        gsize = nz
      end if

      if (pdim(i) == 1) then        ! all local
        lstart(i) = 1
        lend(i)   = gsize
        lsize(i)  = gsize
      elseif (pdim(i) == 2) then    ! distribute across dims(1)
        allocate(st(0:dims(1)-1))
        allocate(en(0:dims(1)-1))
        allocate(sz(0:dims(1)-1))
        call distribute(gsize,dims(1),st,en,sz)
        lstart(i) = st(coord(1))
        lend(i)   = en(coord(1))
        lsize(i)  = sz(coord(1))
        deallocate(st,en,sz)
      elseif (pdim(i) == 3) then    ! distribute across dims(2)
        allocate(st(0:dims(2)-1))
        allocate(en(0:dims(2)-1))
        allocate(sz(0:dims(2)-1))
        call distribute(gsize,dims(2),st,en,sz)
        lstart(i) = st(coord(2))
        lend(i)   = en(coord(2))
        lsize(i)  = sz(coord(2))
        deallocate(st,en,sz)
      end if    

    end do
    return   

  end subroutine partition

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine distribute(data1,proc,st,en,sz)
  
    implicit none
    ! data1 -- data size in any dimension to be partitioned
    ! proc  -- number of processors in that dimension
    ! st    -- array of starting index
    ! en    -- array of ending index
    ! sz    -- array of local size  (redundent)
    integer data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size1,nl,nu
  
    size1=data1/proc
    nu = data1 - size1 * proc
    nl = proc - nu
    st(0) = 1
    sz(0) = size1
    en(0) = size1
    do i=1,nl-1
      st(i) = st(i-1) + size1
      sz(i) = size1
      en(i) = en(i-1) + size1
    end do
    size1 = size1 + 1
    do i=nl,proc-1
      st(i) = en(i-1) + 1
      sz(i) = size1
      en(i) = en(i-1) + size1
    end do
    en(proc-1)= data1 
    sz(proc-1)= data1-st(proc-1)+1
  
    return
  end subroutine distribute

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Define how each dimension is distributed across processors
  !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
  !    such global information is required locally at MPI_ALLTOALLV time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist(nx,ny,nz,decomp)

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    call distribute(nx,dims(1),decomp%x1st,decomp%x1en,decomp%x1dist)
    call distribute(ny,dims(1),decomp%y1st,decomp%y1en,decomp%y1dist)
    call distribute(nz,dims(1),decomp%z1st,decomp%z1en,decomp%z1dist)

    call distribute(nx,dims(2),decomp%x2st,decomp%x2en,decomp%x2dist)
    call distribute(ny,dims(2),decomp%y2st,decomp%y2en,decomp%y2dist)
    call distribute(nz,dims(2),decomp%z2st,decomp%z2en,decomp%z2dist)

    return
  end subroutine get_dist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_buffer(decomp)
    
    implicit none
    
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: i, k
    integer :: rank_x, rank_z
    integer :: subsize_y, offset_y
    integer :: ierror

    ! MPI_ALLTOALLV buffer information

    do i=0, dims(1)-1
       decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
       decomp%z1cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z1dist(i)
       decomp%w1cnts(i) = decomp%x1dist(i)*decomp%wsz(2)*decomp%wsz(3)

       if (i==0) then
          decomp%x1disp(i) = 0  ! displacement is 0-based index
          decomp%y1disp(i) = 0
          decomp%z1disp(i) = 0
          decomp%w1disp(i) = 0
       else
          decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
          decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
          decomp%z1disp(i) = decomp%z1disp(i-1) + decomp%z1cnts(i-1)
          decomp%w1disp(i) = decomp%w1disp(i-1) + decomp%w1cnts(i-1)
       end if
    end do
    
    do i=0, dims(2)-1
       decomp%x2cnts(i) = decomp%x2dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
       decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
       decomp%w2cnts(i) = decomp%x2dist(i)*decomp%wsz(2)*decomp%wsz(3)

       if (i==0) then
          decomp%x2disp(i) = 0  ! displacement is 0-based index
          decomp%y2disp(i) = 0  ! displacement is 0-based index
          decomp%z2disp(i) = 0
          decomp%w2disp(i) = 0
       else
          decomp%x2disp(i) = decomp%x2disp(i-1) + decomp%x2cnts(i-1)
          decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
          decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
          decomp%w2disp(i) = decomp%w2disp(i-1) + decomp%w2cnts(i-1)
       end if
    end do

    ! MPI_ALLTOALL buffer information

    ! For evenly distributed data, following is an easier implementation.
    ! But it should be covered by the more general formulation below.
    !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
    !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1) 
    !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
    !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

    ! For unevenly distributed data, pad smaller messages. Note the 
    ! last blocks along pencils always get assigned more mesh points
    ! for X <=> Y transposes
    decomp%x1count = decomp%x1dist(dims(1)-1) * &
         decomp%y1dist(dims(1)-1) * decomp%xsz(3)
    decomp%y1count = decomp%x1count
    ! for Y <=> Z transposes
    decomp%y2count = decomp%y2dist(dims(2)-1) * &
         decomp%z2dist(dims(2)-1) * decomp%zsz(1)
    decomp%z2count = decomp%y2count

#ifndef USE_CUDA
    ! Information for MPI_Alltoallw for complex X <=> Z transposes
    decomp%xdispls_xz(:)=0
    decomp%zdispls_xz(:)=0
    decomp%xcnts_xz(:)=0
    decomp%zcnts_xz(:)=0
    decomp%xtypes_xz(:)=MPI_DATATYPE_NULL
    decomp%ztypes_xz(:)=MPI_DATATYPE_NULL
    do k=0,dims(1)-1
    do i=0,dims(2)-1
! Actually, rank_x and rank_z are the same..
           call MPI_Cart_rank(DECOMP_2D_COMM_CART_X,(/k,i/),rank_x,ierror)
! this call fails, since DECOMP_2D_COMM_CART_Z is not yet created
!           call MPI_Cart_rank(DECOMP_2D_COMM_CART_Z,(/k,i/),rank_z,ierror)
           rank_z=rank_x
!JD no checks on x- or z-dimension, since we transform from z- into x-pencils, so these always overlap.
           if (decomp%zst(2).le.decomp%y1en(k) .and. &
               decomp%zen(2).ge.decomp%y1st(k)) then
             decomp%zcnts_xz(rank_z+1)=1
             subsize_y=min(decomp%zen(2),decomp%y1en(k))-max(decomp%zst(2),decomp%y1st(k))+1
             offset_y =max(decomp%zst(2),decomp%y1st(k))-decomp%zst(2)
             call MPI_Type_create_subarray(3,decomp%zsz, &
               (/decomp%zsz(1),subsize_y,decomp%z2dist(i)/), &
               (/0,offset_y,decomp%z2st(i)-decomp%zst(3)/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%ztypes_xz(rank_z+1),ierror)
             call MPI_Type_commit(decomp%ztypes_xz(rank_z+1),ierror)
!JD send to process with x-pencil defined by (k,i)
!JD x-bounds are taken from the z-pencils
!             send: decomp%zst(1):decomp%zen(1)
!JD y-bounds are the overlapping region of both pencils.
!                   max(decomp%zst(2),decomp%y1st(k)):min(decomp%zen(2),decomp%y1en(k))
!JD z-bounds are taken from the x-pencils.
!                   decomp%z2st(i):decomp%z2en(i)
           endif
!JD no checks on x- or z-dimension, since we transform from z- into x-pencils, so these always overlap.
           if (decomp%xst(2).le.decomp%y2en(i) .and. &
               decomp%xen(2).ge.decomp%y2st(i)) then
             decomp%xcnts_xz(rank_x+1)=1
             subsize_y=min(decomp%xen(2),decomp%y2en(i))-max(decomp%xst(2),decomp%y2st(i))+1
             offset_y =max(decomp%xst(2),decomp%y2st(i))-decomp%xst(2)
             call MPI_Type_create_subarray(3,decomp%xsz, &
               (/decomp%x1dist(k),subsize_y,decomp%xsz(3)/), &
               (/decomp%x1st(k)-decomp%xst(1),offset_y,0/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%xtypes_xz(rank_x+1),ierror)
             call MPI_Type_commit(decomp%xtypes_xz(rank_x+1),ierror)
!JD recv from process with z-pencil defined by (k,i)
!JD x-bounds are taken from the z-pencils
!             send: decomp%x1st(k):decomp%x1en(k)
!JD y-bounds are the overlapping region of both pencils.
!                   max(decomp%xst(2),decomp%y2st(i)):min(decomp%xen(2),decomp%y2en(i))
!JD z-bounds are taken from the x-pencils.
!                   decomp%xst(3):decomp%xen(3)
           endif
         enddo
         enddo
#endif
    return
  end subroutine prepare_buffer  

#ifdef GLOBAL_ARRAYS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create global arrays that mapped to pencil decompisitions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_global_array(ga, ipencil, data_type, decomp)
    
    implicit none

    integer, intent(OUT) :: ga
    integer, intent(IN) :: ipencil ! 1=X-pencil; 2=Y-pencil; 3=Z-pencil
    integer, intent(IN) :: data_type
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer, dimension(3) :: nblock
    integer, allocatable, dimension(:) :: map
    integer :: offset, i, errorcode
    logical :: success

    ga = ga_create_handle()
    call ga_set_data(ga, 3, &
         (/decomp%xsz(1),decomp%ysz(2),decomp%zsz(3)/), data_type)
    allocate(map(1+dims(1)+dims(2)))

    ! generate the GA irreg distribution parameters using 
    ! 2DECOMP's decomposition information
    if (ipencil==1) then  ! X-pencil
       nblock(1) = 1
       nblock(2) = dims(1)
       nblock(3) = dims(2)
       map(1) = 1
       offset = nblock(1)+1
       do i=0, dims(1)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%y1dist(i-1)
          end if
       end do
       offset = nblock(1) + nblock(2) + 1
       do i=0, dims(2)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%z2dist(i-1)
          end if
       end do
    else if (ipencil==2) then  ! Y-pencil
       nblock(1) = dims(1)
       nblock(2) = 1
       nblock(3) = dims(2)
       offset = 1
       do i=0, dims(1)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%x1dist(i-1)
          end if
       end do
       map(nblock(1)+1) = 1
       offset = nblock(1) + nblock(2) + 1
       do i=0, dims(2)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%z2dist(i-1)
          end if
       end do
    else if (ipencil==3) then  ! Z-pencil
       nblock(1) = dims(1)
       nblock(2) = dims(2)
       nblock(3) = 1
       offset = 1
       do i=0, dims(1)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%x1dist(i-1)
          end if
       end do
       offset = nblock(1)+1
       do i=0, dims(2)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%y2dist(i-1)
          end if
       end do
       map(nblock(1)+nblock(2)+1) = 1
    end if

    call ga_set_irreg_distr(ga, map, nblock)
    success = ga_allocate(ga)
    if (.not.success) then
       errorcode = 7
       call decomp_2d_abort(errorcode, &
            'Failed to create global arrays')
    end if

    deallocate(map)

    return
  end subroutine get_global_array

#endif


#ifdef OCC
  ! For non-blocking communication code, progress the comminication stack
  subroutine transpose_test(handle)

    implicit none

    integer :: handle, ierror

    call NBC_TEST(handle,ierror)

    return
  end subroutine transpose_test
#endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transposition routines 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define EPA2A
#include "transpose_x_to_y.F90"
#include "transpose_y_to_z.F90"
#include "transpose_z_to_y.F90"
#include "transpose_y_to_x.F90"
#include "transpose_x_to_z.F90"
#include "transpose_z_to_x.F90"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Auto-tuning algorithm to select the best 2D processor grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)
    use param, only: lvlhalo
    use nvtx
    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, allocatable, dimension(:) :: factors
    double precision :: t1, t2, best_time
    integer :: nfact, i, row, col, ierror, errorcode

#ifdef USE_CUDA
    real(mytype),    device, allocatable, dimension(:,:,:) :: rx, ry, hx
    complex(mytype), device, allocatable, dimension(:,:,:) :: cx, cy, cz 
    TYPE(DECOMP_INFO) :: decomp_r, decomp_c
    integer :: istat
#else
    real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3
    TYPE(DECOMP_INFO) :: decomp
#endif

call nvtxStartRange("best_2d_grid",1)

    if (nrank==0) write(*,*) 'In auto-tuning mode......'

    best_time = huge(t1)
    best_p_row = -1
    best_p_col = -1
    
    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

#ifdef USE_CUDA
    alloc_work_r = 0
    alloc_work_c = 0

    do i=1, nfact
       row = factors(i)
       col = iproc / row

       ! enforce the limitation of 2D decomposition
       if (min(nx_global,ny_global/2+1)>=row .and. &
            min(ny_global/2+1,nz_global)>=col) then

          ! 2D Catersian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .true.
          periodic(2) = .true.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,DECOMP_2D_COMM_CART_X, ierror)
          if(ierror) print *,nrank,"ERROR cart create"
          call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
          if(ierror) print *,nrank,"ERROR cart coords"
          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
               DECOMP_2D_COMM_COL,ierror)
          if(ierror) print *,nrank,"ERROR cart sub"
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
               DECOMP_2D_COMM_ROW,ierror)
          if(ierror) print *,nrank,"ERROR cart sub"

          call MPI_COMM_RANK(DECOMP_2D_COMM_COL,col_rank,ierror); if(ierror) print *,nrank,"ERROR commrank"
          call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,row_rank,ierror); if(ierror) print *,nrank,"ERROR commrank"

          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global    ,nz_global,decomp_r)
          call decomp_info_init(nx_global,ny_global/2+1,nz_global,decomp_c)
          call decomp_info_finalize(decomp_r)
          call decomp_info_finalize(decomp_c)
          call MPI_Comm_free(DECOMP_2D_COMM_COL,ierror) ; if(ierror) print *,nrank,"ERROR commfree"
          call MPI_Comm_free(DECOMP_2D_COMM_ROW,ierror) ; if(ierror) print *,nrank,"ERROR commfree"
          call MPI_Comm_free(DECOMP_2D_COMM_CART_X, ierror) ; if(ierror) print *,nrank,"ERROR commfree"
      end if
    end do
    alloc_work_r = 1
    alloc_work_c = 1
#endif

!    Make initial communication to un-bias results

    row = factors(1)
    col = iproc / row

       if (min(nx_global,ny_global/2+1)>=row .and. &
            min(ny_global/2+1,nz_global)>=col) then

          ! 2D Catersian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .true.
          periodic(2) = .true.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,DECOMP_2D_COMM_CART_X, ierror)
          if(ierror) print *,nrank,"ERROR cart create"
          call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
          if(ierror) print *,nrank,"ERROR cart coords"
          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
               DECOMP_2D_COMM_COL,ierror)
          if(ierror) print *,nrank,"ERROR cart sub"
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
               DECOMP_2D_COMM_ROW,ierror)
          if(ierror) print *,nrank,"ERROR cart sub"


#ifdef USE_CUDA
    call MPI_COMM_RANK(DECOMP_2D_COMM_COL,col_rank,ierror); if(ierror) print *,nrank,"ERROR commrank"
    call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,row_rank,ierror); if(ierror) print *,nrank,"ERROR commrank"
    call init_neighbour
          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global    ,nz_global,decomp_r)
          call decomp_info_init(nx_global,ny_global/2+1,nz_global,decomp_c)

          ! arrays for X,Y and Z-pencils
          allocate(rx(decomp_r%xst(1):decomp_r%xen(1),decomp_r%xst(2):decomp_r%xen(2),decomp_r%xst(3):decomp_r%xen(3)))
          allocate(cx(decomp_c%wst(1):decomp_c%wen(1),decomp_c%wst(2):decomp_c%wen(2),decomp_c%wst(3):decomp_c%wen(3)))
          allocate(ry(decomp_r%yst(1):decomp_r%yen(1),decomp_r%yst(2):decomp_r%yen(2),decomp_r%yst(3):decomp_r%yen(3)))
          allocate(cy(decomp_c%yst(1):decomp_c%yen(1),decomp_c%yst(2):decomp_c%yen(2),decomp_c%yst(3):decomp_c%yen(3)))
          allocate(cz(decomp_c%zst(1):decomp_c%zen(1),decomp_c%zst(2):decomp_c%zen(2),decomp_c%zst(3):decomp_c%zen(3)))
          allocate(hx(decomp_r%xst(1):decomp_r%xen(1),decomp_r%xst(2)-lvlhalo:decomp_r%xen(2)+lvlhalo,decomp_r%xst(3)-lvlhalo:decomp_r%xen(3)+lvlhalo))

          rx=0.d0; cx=(0.d0,0.d0); ry=0.d0; cy=(0.d0,0.d0); cz=(0.d0,0.d0); hx=0.0d0

          istat=cudaDeviceSynchronize()
          call MPI_BARRIER(MPI_COMM_WORLD,ierror)

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(rx,ry,decomp_r)
          call transpose_y_to_z(cy,cz,decomp_c)
          call transpose_z_to_x(cz,cx,decomp_c)
          call transpose_x_to_z(cx,cz,decomp_c)
          call transpose_z_to_y(cz,cy,decomp_c)
          call transpose_y_to_x(ry,rx,decomp_r)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          t2 = MPI_WTIME() - t1

          !deallocate(rx,ry)
          deallocate(rx,cx,ry,cy,cz,hx)
          call decomp_info_finalize(decomp_r)
          call decomp_info_finalize(decomp_c)
#else          
          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global,nz_global,decomp)

          ! arrays for X,Y and Z-pencils
          allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
          allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(u1,u2,decomp)
          call transpose_y_to_x(u2,u1,decomp)
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2)
          call decomp_info_finalize(decomp)
#endif
          call MPI_Comm_free(DECOMP_2D_COMM_COL,ierror) ; if(ierror) print *,nrank,"ERROR commfree"
          call MPI_Comm_free(DECOMP_2D_COMM_ROW,ierror) ; if(ierror) print *,nrank,"ERROR commfree"
          call MPI_Comm_free(DECOMP_2D_COMM_CART_X, ierror) ; if(ierror) print *,nrank,"ERROR commfree"

       end if

    do i=1, nfact

       row = factors(i)
       col = iproc / row

       ! enforce the limitation of 2D decomposition
       if (min(nx_global,ny_global/2+1)>=row .and. &
            min(ny_global/2+1,nz_global)>=col) then

          ! 2D Catersian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .true.
          periodic(2) = .true.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,DECOMP_2D_COMM_CART_X, ierror)
          if(ierror) print *,nrank,"ERROR cart create"
          call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
          if(ierror) print *,nrank,"ERROR cart coords"
          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
               DECOMP_2D_COMM_COL,ierror)
          if(ierror) print *,nrank,"ERROR cart sub"
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
               DECOMP_2D_COMM_ROW,ierror)
          if(ierror) print *,nrank,"ERROR cart sub"
          
#ifdef USE_CUDA
    call MPI_COMM_RANK(DECOMP_2D_COMM_COL,col_rank,ierror); if(ierror) print *,nrank,"ERROR commrank"
    call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,row_rank,ierror); if(ierror) print *,nrank,"ERROR commrank"
    call init_neighbour

          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global    ,nz_global,decomp_r)
          call decomp_info_init(nx_global,ny_global/2+1,nz_global,decomp_c)

          ! arrays for X,Y and Z-pencils
          allocate(rx(decomp_r%xst(1):decomp_r%xen(1),decomp_r%xst(2):decomp_r%xen(2),decomp_r%xst(3):decomp_r%xen(3)))
          allocate(cx(decomp_c%wst(1):decomp_c%wen(1),decomp_c%wst(2):decomp_c%wen(2),decomp_c%wst(3):decomp_c%wen(3)))
          allocate(ry(decomp_r%yst(1):decomp_r%yen(1),decomp_r%yst(2):decomp_r%yen(2),decomp_r%yst(3):decomp_r%yen(3)))
          allocate(cy(decomp_c%yst(1):decomp_c%yen(1),decomp_c%yst(2):decomp_c%yen(2),decomp_c%yst(3):decomp_c%yen(3)))
          allocate(cz(decomp_c%zst(1):decomp_c%zen(1),decomp_c%zst(2):decomp_c%zen(2),decomp_c%zst(3):decomp_c%zen(3)))
          allocate(hx(decomp_r%xst(1):decomp_r%xen(1),decomp_r%xst(2)-lvlhalo:decomp_r%xen(2)+lvlhalo,decomp_r%xst(3)-lvlhalo:decomp_r%xen(3)+lvlhalo))

          rx=0.d0; cx=(0.d0,0.d0); ry=0.d0; cy=(0.d0,0.d0); cz=(0.d0,0.d0); hx=0.0d0

          istat=cudaDeviceSynchronize()
          call MPI_BARRIER(MPI_COMM_WORLD,ierror)

call nvtxStartRange("try",i)

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(rx,ry,decomp_r)
          call transpose_y_to_z(cy,cz,decomp_c)
          call transpose_z_to_x(cz,cx,decomp_c)
          call transpose_x_to_z(cx,cz,decomp_c)
          call transpose_z_to_y(cz,cy,decomp_c)
          call transpose_y_to_x(ry,rx,decomp_r)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          call update_halo(hx,lvlhalo)
          t2 = MPI_WTIME() - t1

call nvtxEndRange

          !deallocate(rx,ry)
          deallocate(rx,cx,ry,cy,cz,hx)
          call decomp_info_finalize(decomp_r)
          call decomp_info_finalize(decomp_c)
#else
          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global,nz_global,decomp)

          ! arrays for X,Y and Z-pencils
          allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
          allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          allocate(u3(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(u1,u2,decomp)
          call transpose_y_to_z(u2,u3,decomp)
          call transpose_z_to_y(u3,u2,decomp)
          call transpose_y_to_x(u2,u1,decomp)
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2,u3)
          call decomp_info_finalize(decomp)
#endif
          !call MPI_Comm_free(DECOMP_2D_COMM_COL,ierror)
          !call MPI_Comm_free(DECOMP_2D_COMM_ROW,ierror)
          !call MPI_Comm_free(DECOMP_2D_COMM_CART_X, ierror)
          call MPI_Comm_free(DECOMP_2D_COMM_COL,ierror) ; if(ierror) print *,nrank,"ERROR commfree"
          call MPI_Comm_free(DECOMP_2D_COMM_ROW,ierror) ; if(ierror) print *,nrank,"ERROR commfree"
          call MPI_Comm_free(DECOMP_2D_COMM_CART_X, ierror) ; if(ierror) print *,nrank,"ERROR commfree"


          call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                   MPI_COMM_WORLD,ierror)
          if(ierror) print *,nrank,"ERROR allreduce"
          t1 = t1 / dble(nproc)

          if (nrank==0) then
             write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
          end if

          if (best_time > t1) then
             best_time = t1
             best_p_row = row
             best_p_col = col
          end if

       end if
       
    end do ! loop through processer grid

    deallocate(factors)

    if (best_p_row/=-1) then
       if (nrank==0) then
          write(*,*) 'the best processor grid is probably ', &
               best_p_row, ' by ', best_p_col
       end if
    else
       errorcode = 9
       call decomp_2d_abort(errorcode, &
            'The processor-grid auto-tuning code failed. ' // &
            'The number of processes requested is probably too large.')
    end if

call nvtxEndRange

    return
  end subroutine best_2d_grid

#include "factor.F90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.F90"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error handling
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_abort(errorcode, msg)

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    integer :: ierror
    
    if (nrank==0) then
       write(*,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
    end if
    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)

    return
  end subroutine decomp_2d_abort


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.F90"
    
  
end module decomp_2d

