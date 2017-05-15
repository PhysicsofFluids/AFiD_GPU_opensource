#ifdef USE_CUDA
module dph_routines
  contains

  subroutine copy_dph_to_dphhalo(dph,dphhalo)
    use nvtx
    use param, only: fp_kind, nxm, lvlhalo
    use decomp_2d, only: xstart, xend, update_halo

    implicit none
    !args
    real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),intent(IN) :: dph 
    real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(OUT) :: dphhalo
#ifdef USE_CUDA
    attributes(device) :: dph, dphhalo
#endif
    !locals
    integer :: i,j,k

    call nvtxStartRangeAsync("COPYdph ",9)
    !$cuf kernel do(3) <<<*,*>>>
    do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
        do k=1,nxm
          dphhalo(k,j,i) = dph(k,j,i)
        enddo
      enddo
    enddo
    call nvtxEndRangeAsync

  end subroutine copy_dph_to_dphhalo

  subroutine update_dphhalo(dph,dphhalo)
    use nvtx
    use param, only: fp_kind, nxm, lvlhalo
    use decomp_2d, only: xstart, xend, update_halo

    implicit none
    !args
    real(fp_kind), dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),intent(IN) :: dph 
    real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo),intent(OUT) :: dphhalo
#ifdef USE_CUDA
    attributes(device) :: dph, dphhalo
#endif
    call nvtxStartRangeAsync("update_halo ",13)
    call update_halo(dphhalo,lvlhalo)
    call nvtxEndRangeAsync
  end subroutine update_dphhalo
    
#ifdef USE_HYBRID
    subroutine transfer_dphhalo_bulk_d2h(dphhalo_h, dph_d)
      use cudafor
      use param, only: fp_kind, nxm, nym, lvlhalo, istat
      use decomp_2d, only: xstart, xend, xstart_cpu, xsize_cpu, a2a_d2h
      implicit none
      real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: dphhalo_h
      real(fp_kind), device, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(IN) :: dph_d
      integer :: k

      ! Copy bulk data plane by plane from device to host due to halo points
      do k = xstart_cpu(3)-1, xend(3)
        istat = cudaMemcpyAsync(dphhalo_h(1, xstart(2), k), dph_d(1, xstart(2), k), nxm * nym, stream = a2a_d2h)
      end do

    end subroutine transfer_dphhalo_bulk_d2h

    subroutine transfer_dphhalo_halo_d2h(in_h, in_d)
      use cudafor
      use param, only: fp_kind, nxm, nym, lvlhalo, istat
      use decomp_2d, only: xstart, xend, xstart_cpu, xsize_cpu, a2a_d2h
      implicit none
      real(fp_kind), dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: in_h
      real(fp_kind), device, dimension(1:nxm,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo), intent(IN) :: in_d

      istat = cudaMemcpyAsync(in_h(1, xstart(2)-lvlhalo, xend(3)+lvlhalo), in_d(1, xstart(2)-lvlhalo, xend(3)+lvlhalo), nxm * (nym + 2), stream = a2a_d2h)

      istat = cudaMemcpy2DAsync(in_h(1, xstart(2)-lvlhalo, xstart_cpu(3)-lvlhalo), nxm * (nym + 2), in_d(1, xstart(2)-lvlhalo, &
                                xstart_cpu(3)-lvlhalo), nxm * (nym + 2), nxm, xsize_cpu(3)+2, stream = a2a_d2h)

      istat = cudaMemcpy2DAsync(in_h(1, xend(2)+lvlhalo, xstart_cpu(3)-lvlhalo), nxm * (nym + 2), in_d(1, xend(2)+lvlhalo, &
                                xstart_cpu(3)-lvlhalo), nxm * (nym + 2), nxm, xsize_cpu(3)+2, stream = a2a_d2h)
    end subroutine transfer_dphhalo_halo_d2h
#endif
end module dph_routines
#endif
