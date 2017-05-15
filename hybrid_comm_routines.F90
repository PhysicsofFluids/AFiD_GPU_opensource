#ifdef USE_HYBRID
module hybrid_comm_routines
  contains
  subroutine update_interface(in_h, in_d)
    use cudafor
    use param, only: fp_kind, nx, nym, lvlhalo, istat
    use decomp_2d, only: xstart, xend, xstart_cpu, xend_gpu, a2a_d2h, a2a_h2d
    implicit none
    real(fp_kind), dimension(1:nx, xstart(2)-lvlhalo:xend(2)+lvlhalo, xstart(3)-lvlhalo:xend(3)+lvlhalo) :: in_h
    real(fp_kind), device, dimension(1:nx, xstart(2)-lvlhalo:xend(2)+lvlhalo, xstart(3)-lvlhalo:xend_gpu(3)+lvlhalo) :: in_d

    !JR GPU->CPU
    istat = cudaMemcpyAsync(in_h(1, xstart(2)-lvlhalo, xstart_cpu(3)-1), in_d(1, xstart(2)-lvlhalo, xend_gpu(3)), nx * (nym + 2), stream = a2a_d2h)
    !JR CPU->GPU
    istat = cudaMemcpyAsync(in_d(1, xstart(2)-lvlhalo, xend_gpu(3)+1), in_h(1, xstart(2)-lvlhalo, xstart_cpu(3)), nx * (nym + 2), stream = a2a_h2d)

  end subroutine update_interface
  subroutine transfer_halo_d2h(in_h, in_d)
    use cudafor
    use param, only: fp_kind, nx, nym, lvlhalo, istat
    use decomp_2d, only: xstart, xend, xend_gpu, xsize_gpu, a2a_d2h
    implicit none
    real(fp_kind), dimension(1:nx, xstart(2)-lvlhalo:xend(2)+lvlhalo, xstart(3)-lvlhalo:xend(3)+lvlhalo) :: in_h
    real(fp_kind), device, dimension(1:nx, xstart(2)-lvlhalo:xend(2)+lvlhalo, xstart(3)-lvlhalo:xend_gpu(3)+lvlhalo), intent(IN) :: in_d

    !JR Before CPU halo exchanges, copy boundary planes for halo communication from GPU to CPU
    istat = cudaMemcpyAsync(in_h(1, xstart(2)-lvlhalo, xstart(3)), in_d(1, xstart(2)-lvlhalo, xstart(3)), nx * (nym + 2), stream = a2a_d2h)

    istat = cudaMemcpy2DAsync(in_h(1, xstart(2), xstart(3)), nx * (nym + 2), in_d(1, xstart(2), xstart(3)), &
                           nx * (nym + 2), nx, xsize_gpu(3), stream = a2a_d2h)
    istat = cudaMemcpy2DAsync(in_h(1, xend(2), xstart(3)), nx * (nym + 2), in_d(1, xend(2), xstart(3)), &
                         nx * (nym + 2), nx, xsize_gpu(3), stream = a2a_d2h)

  end subroutine transfer_halo_d2h

  subroutine transfer_halo_h2d(in_d, in_h)
    use cudafor
    use param, only: fp_kind, nx, nym, lvlhalo, istat
    use decomp_2d, only: xstart, xend, xend_gpu, xsize_gpu, a2a_d2h, a2a_h2d
    implicit none
    real(fp_kind), dimension(1:nx, xstart(2)-lvlhalo:xend(2)+lvlhalo, xstart(3)-lvlhalo:xend(3)+lvlhalo), intent(IN) :: in_h
    real(fp_kind), device, dimension(1:nx, xstart(2)-lvlhalo:xend(2)+lvlhalo, xstart(3)-lvlhalo:xend_gpu(3)+lvlhalo) :: in_d

    !JR After CPU halo exchange, copy updated halo planes from CPU to GPU
    istat = cudaMemcpyAsync(in_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), in_h(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2), stream = a2a_h2d)

    istat = cudaMemcpy2DAsync(in_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2), in_h(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), &
                         nx * (nym + 2), nx, xsize_gpu(3)+2, stream = a2a_h2d)
    istat = cudaMemcpy2DAsync(in_d(1, xend(2)+lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2), in_h(1, xend(2)+lvlhalo, xstart(3)-lvlhalo), &
                         nx * (nym + 2), nx, xsize_gpu(3)+2, stream = a2a_h2d)
  end subroutine transfer_halo_h2d
end module hybrid_comm_routines
#endif

