! cuDecomp adapter, included in decomp_2d's CONTAINS section.
! Keep the solver's physical (x,y,z) memory order for every pencil.
  subroutine cu_check(status, operation)
    integer, intent(in) :: status
    character(len=*), intent(in) :: operation
    if (status /= CUDECOMP_RESULT_SUCCESS) then
      call decomp_2d_abort(status, 'cuDecomp: '//operation)
    end if
  end subroutine cu_check

  ! One selector controls transposes and the matching halo transport.
  ! Called after MPI initialization and before any descriptors are created.
  subroutine cu_select_backend(backend)
    character(len=*), intent(in), optional :: backend
    character(len=:), allocatable :: name
    character(len=16) :: label
    integer :: i, code

    name = 'mpi-p2p'
    if (present(backend)) name = trim(adjustl(backend))
    do i=1,len(name)
      code = iachar(name(i:i))
      if (code>=iachar('A') .and. code<=iachar('Z')) name(i:i) = achar(code+32)
    end do
    cu_halo_backend = CUDECOMP_HALO_COMM_MPI
    select case (name)
    case ('1','mpi-p2p')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
      label = 'mpi-p2p'
    case ('2','mpi-p2p-pl')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL
      label = 'mpi-p2p-pl'
    case ('3','mpi-a2a')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_MPI_A2A
      label = 'mpi-a2a'
    case ('4','nccl')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
      cu_halo_backend = CUDECOMP_HALO_COMM_NCCL
      label = 'nccl'
    case ('5','nccl-pl')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_NCCL_PL
      cu_halo_backend = CUDECOMP_HALO_COMM_NCCL
      label = 'nccl-pl'
    case ('6','nvshmem')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM
      cu_halo_backend = CUDECOMP_HALO_COMM_NVSHMEM
      label = 'nvshmem'
    case ('7','nvshmem-pl')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
      cu_halo_backend = CUDECOMP_HALO_COMM_NVSHMEM
      label = 'nvshmem-pl'
    case ('8','nvshmem-sm')
      cu_transpose_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_SM
      cu_halo_backend = CUDECOMP_HALO_COMM_NVSHMEM
      label = 'nvshmem-sm'
    case default
      call decomp_2d_abort(6, 'Unknown cuDecomp backend: '//name// &
        '. Use 1:mpi-p2p, 2:mpi-p2p-pl, 3:mpi-a2a, 4:nccl, 5:nccl-pl, '// &
        '6:nvshmem, 7:nvshmem-pl, or 8:nvshmem-sm')
      return
    end select
    if (nrank==0) write(*,'(A,A,A,I0)') 'cuDecomp backend: ', trim(label), &
      '; halo backend ID: ', cu_halo_backend
  end subroutine cu_select_backend

  subroutine cu_create_grid(gdims, grid, swapped)
    integer, intent(in) :: gdims(3)
    type(cudecompGridDesc), intent(out) :: grid
    logical, intent(in) :: swapped
    type(cudecompGridDescConfig) :: config
    integer :: axis

    call cu_check(cudecompGridDescConfigSetDefaults(config), 'configuration defaults')
    config%gdims = gdims
    config%pdims = dims
    config%rank_order = CUDECOMP_RANK_ORDER_ROW_MAJOR
    config%transpose_comm_backend = cu_transpose_backend
    config%halo_comm_backend = cu_halo_backend
    config%transpose_axis_contiguous = .false.
    if (swapped) then
      ! AFiD's spectral W pencil distributes z over rows and y over columns.
      ! In coordinates (x,z,y), W is X and the physical Z pencil is Y.
      ! The memory order below retains physical (x,y,z) storage in both.
      do axis=1,3
        config%transpose_mem_order(:,axis) = [1,3,2]
      end do
    end if
    call cu_check(cudecompGridDescCreate(cu_handle, grid, config), 'create grid')
  end subroutine cu_create_grid


  subroutine cu_decomp_init(nx, ny, nz, decomp)
    use iso_fortran_env, only: int64
    integer, intent(in) :: nx, ny, nz
    type(DECOMP_INFO), intent(inout) :: decomp
    type(cudecompPencilInfo) :: pencil
    integer(int64) :: nwork
    integer :: axis

    call cu_create_grid([nx,ny,nz], decomp%cu_grid, .false.)
    call cu_create_grid([nx,nz,ny], decomp%cu_grid_xzy, .true.)
    ! Validate starts as well as sizes: an offset mismatch corrupts uneven grids.
    do axis=1,3
      call cu_check(cudecompGetPencilInfo(cu_handle, decomp%cu_grid, pencil, axis), 'pencil geometry')
      select case(axis)
      case(1)
        call cu_validate_pencil(pencil, decomp%xst, decomp%xsz)
      case(2)
        call cu_validate_pencil(pencil, decomp%yst, decomp%ysz)
      case(3)
        call cu_validate_pencil(pencil, decomp%zst, decomp%zsz)
      end select
    end do
    call cu_check(cudecompGetPencilInfo(cu_handle, decomp%cu_grid_xzy, pencil, 1), 'W pencil geometry')
    call cu_validate_pencil(pencil, decomp%wst, decomp%wsz)
    call cu_check(cudecompGetPencilInfo(cu_handle, decomp%cu_grid_xzy, pencil, 2), 'Z pencil geometry')
    call cu_validate_pencil(pencil, decomp%zst, decomp%zsz)
    call cu_check(cudecompGetTransposeWorkspaceSize(cu_handle, decomp%cu_grid, nwork), 'workspace size')
    ! Complex allocation also covers real transposes; sizes are in elements.
    call cu_check(cudecompMalloc(cu_handle, decomp%cu_grid, decomp%cu_work, max(1_int64,nwork)), 'workspace')
    call cu_check(cudecompGetTransposeWorkspaceSize(cu_handle, decomp%cu_grid_xzy, nwork), 'xzy workspace size')
    call cu_check(cudecompMalloc(cu_handle, decomp%cu_grid_xzy, decomp%cu_work_xzy, max(1_int64,nwork)), &
                  'xzy workspace')
  end subroutine cu_decomp_init

  subroutine cu_validate_pencil(pencil, first, extent)
    type(cudecompPencilInfo), intent(in) :: pencil
    integer, intent(in) :: first(3), extent(3)
    if (any(pencil%lo /= first) .or. any(pencil%shape /= extent)) then
      call decomp_2d_abort(6, 'cuDecomp pencil geometry does not match solver decomposition')
    end if
  end subroutine cu_validate_pencil

  subroutine cu_halo_finalize()
    if (cu_halo_level < 0) return
    call cu_check(cudecompFree(cu_handle, cu_halo_grid, cu_halo_work), 'free halo workspace')
    call cu_check(cudecompGridDescDestroy(cu_handle, cu_halo_grid), 'destroy halo grid')
    cu_halo_level = -1
    cu_halo_dims = 0
  end subroutine cu_halo_finalize

  subroutine cu_update_halo(field, level)
    use iso_fortran_env, only: int64
    real(mytype), device, contiguous, intent(inout) :: field(:,:,:)
    integer, intent(in) :: level
    integer :: halos(3), halo_dims(3)
    integer(int64) :: nwork
    type(cudecompPencilInfo) :: pencil
    logical :: periods(3)

    if (level == 0) return
    if (level < 0) call decomp_2d_abort(6, 'Negative halo extent')
    halos = [0,level,level]
    periods = [periodic_x,periodic_y,periodic_z]
    ! X is local and includes the extra boundary plane for solver fields.
    ! Pressure corrections and tuning buffers instead have nx_global entries.
    halo_dims = [size(field,1),dims]
    ! Also used during process-grid tuning, before decomp_main exists.
    if (any(cu_halo_dims /= halo_dims) .or. cu_halo_level /= level) then
      call cu_halo_finalize()
      call cu_create_grid([size(field,1),ny_global,nz_global], cu_halo_grid, .false.)
      call cu_check(cudecompGetHaloWorkspaceSize(cu_handle, cu_halo_grid, 1, halos, nwork), 'halo workspace size')
      call cu_check(cudecompMalloc(cu_handle, cu_halo_grid, cu_halo_work, max(1_int64,nwork)), 'halo workspace')
      cu_halo_dims = halo_dims
      cu_halo_level = level
    end if
    call cu_check(cudecompGetPencilInfo(cu_handle, cu_halo_grid, pencil, 1, halos), 'halo geometry')
    if (any(shape(field) /= pencil%shape)) call decomp_2d_abort(6, 'cuDecomp halo array shape mismatch')
    ! Complete producer kernels on any solver stream before communication.
    call cu_check(cudaDeviceSynchronize(), 'synchronize before halo exchange')
    ! Y before Z propagates edges and corners, matching the existing stencil.
    call cu_check(cudecompUpdateHalosX(cu_handle, cu_halo_grid, field, cu_halo_work, &
                  CUDECOMP_DOUBLE, halos, periods, 2), 'Y halo exchange')
    call cu_check(cudecompUpdateHalosX(cu_handle, cu_halo_grid, field, cu_halo_work, &
                  CUDECOMP_DOUBLE, halos, periods, 3), 'Z halo exchange')
    call cu_check(cudaDeviceSynchronize(), 'complete halo exchange')
  end subroutine cu_update_halo

  subroutine cu_transpose_real(src, dst, decomp, source_axis, dest_axis)
    real(mytype), device, contiguous, intent(in) :: src(:,:,:)
    real(mytype), device, contiguous, intent(out) :: dst(:,:,:)
    type(DECOMP_INFO), intent(in) :: decomp
    integer, intent(in) :: source_axis, dest_axis
    real(mytype), device, allocatable :: middle(:,:,:)
    call cu_check(cudaDeviceSynchronize(), 'synchronize before transpose')
    select case(10*source_axis + dest_axis)
    case(12)
      call cu_check(cudecompTransposeXToY(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose XToY')
    case(21)
      call cu_check(cudecompTransposeYToX(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose YToX')
    case(23)
      call cu_check(cudecompTransposeYToZ(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose YToZ')
    case(32)
      call cu_check(cudecompTransposeZToY(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose ZToY')
    case(13)
      allocate(middle(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
      call cu_check(cudecompTransposeXToY(cu_handle, decomp%cu_grid, src, middle, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose XToY')
      call cu_check(cudecompTransposeYToZ(cu_handle, decomp%cu_grid, middle, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose YToZ')
      call cu_check(cudaDeviceSynchronize(), 'complete intermediate transpose')
      deallocate(middle)
    case(31)
      allocate(middle(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
      call cu_check(cudecompTransposeZToY(cu_handle, decomp%cu_grid, src, middle, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose ZToY')
      call cu_check(cudecompTransposeYToX(cu_handle, decomp%cu_grid, middle, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE), 'transpose YToX')
      call cu_check(cudaDeviceSynchronize(), 'complete intermediate transpose')
      deallocate(middle)
    case default
      call decomp_2d_abort(6, 'Unsupported cuDecomp transpose')
    end select
    call cu_check(cudaDeviceSynchronize(), 'complete transpose')
  end subroutine cu_transpose_real

  subroutine cu_transpose_complex(src, dst, decomp, source_axis, dest_axis)
    complex(mytype), device, contiguous, intent(in) :: src(:,:,:)
    complex(mytype), device, contiguous, intent(out) :: dst(:,:,:)
    type(DECOMP_INFO), intent(in) :: decomp
    integer, intent(in) :: source_axis, dest_axis
    call cu_check(cudaDeviceSynchronize(), 'synchronize before transpose')
    select case(10*source_axis + dest_axis)
    case(12)
      call cu_check(cudecompTransposeXToY(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE_COMPLEX), 'transpose XToY')
    case(21)
      call cu_check(cudecompTransposeYToX(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE_COMPLEX), 'transpose YToX')
    case(23)
      call cu_check(cudecompTransposeYToZ(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE_COMPLEX), 'transpose YToZ')
    case(32)
      call cu_check(cudecompTransposeZToY(cu_handle, decomp%cu_grid, src, dst, &
                    decomp%cu_work, CUDECOMP_DOUBLE_COMPLEX), 'transpose ZToY')
    case(13)
      call cu_check(cudecompTransposeXToY(cu_handle, decomp%cu_grid_xzy, src, dst, &
                    decomp%cu_work_xzy, CUDECOMP_DOUBLE_COMPLEX), 'spectral transpose XToY')
    case(31)
      call cu_check(cudecompTransposeYToX(cu_handle, decomp%cu_grid_xzy, src, dst, &
                    decomp%cu_work_xzy, CUDECOMP_DOUBLE_COMPLEX), 'spectral transpose YToX')
    case default
      call decomp_2d_abort(6, 'Unsupported cuDecomp transpose')
    end select
    call cu_check(cudaDeviceSynchronize(), 'complete transpose')
  end subroutine cu_transpose_complex
