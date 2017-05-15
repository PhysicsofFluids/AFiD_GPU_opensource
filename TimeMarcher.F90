!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         ! 
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TimeMarcher
  use nvtx

#if defined(USE_CUDA) && defined(USE_HYBRID)
  use cudafor
  use dph_routines
  use hybrid_comm_routines
  use param, only : fp_kind,nsst,alm,gam,rom,al,ga,ro,beta,dt,ren,lvlhalo,nxm,nx,nym,ap3ck_h=>ap3ck,ap3ck=>ap3ck_d,ac3ck_h=>ac3ck,ac3ck=>ac3ck_d,am3ck_h=>am3ck,am3ck=>am3ck_d,     &
                    ap3sk_h=>ap3sk,ap3sk=>ap3sk_d,ac3sk_h=>ac3sk,ac3sk=>ac3sk_d,am3ssk_h=>am3ssk,am3ssk=>am3ssk_d,ac3ssk_h=>ac3ssk,ac3ssk=>ac3ssk_d,ap3ssk_h=>ap3ssk,ap3ssk=>ap3ssk_d, &
                    am3sk_h=>am3sk,am3sk=>am3sk_d,amphk_h=>amphk,amphk=>amphk_d,acphk_h=>acphk,acphk=>acphk_d,apphk_h=>apphk,apphk=>apphk_d, &
                    tempbp_h=>tempbp,tempbp=>tempbp_d,temptp_h=>temptp,temptp=>temptp_d,udx3c_h=>udx3c,udx3c=>udx3c_d,udx3m_h=>udx3m,udx3m=>udx3m_d,kmv_h=>kmv, &
                    kmv=>kmv_d,kpv_h=>kpv,kpv=>kpv_d, istat
  use local_arrays, only : vx_h=>vx,vx=>vx_d,vy_h=>vy,vy=>vy_d,vz_h=>vz,vz=>vz_d,temp_h=>temp,temp=>temp_d,pr_h=>pr, pr=>pr_d,qcap_h=>qcap,dq_h=>dq, &
                           hro_h=>hro,hro=>hro_d,dph_h=>dph,dph=>dph_d,rux_h=>rux,rux=>rux_d,ruy_h=>ruy,ruy=>ruy_d,ruz_h=>ruz,ruz=>ruz_d,rutemp_h=>rutemp,rutemp=>rutemp_d, &
                           rhs_h=>rhs,dphhalo_h=>dphhalo
  use decomp_2d, only: xstart, xend, xstart_cpu, xend_cpu, xstart_gpu, xend_gpu, a2a_h2d, a2a_d2h, a2a_comp, rhs=>work1_r_d, dphhalo=>work2_r_d
  use fftw_params, only: qcap=>rbuff1_d, dq=>rbuff2_d

#elif defined(USE_CUDA)
  use cudafor
  use dph_routines
  use param, only : fp_kind,nsst,alm,gam,rom,al,ga,ro,beta,dt,ren,lvlhalo,nxm,nx,ap3ck=>ap3ck_d,ac3ck=>ac3ck_d,am3ck=>am3ck_d,     &
                    ap3sk=>ap3sk_d,ac3sk=>ac3sk_d,am3ssk=>am3ssk_d,ac3ssk=>ac3ssk_d,ap3ssk=>ap3ssk_d,am3sk=>am3sk_d,amphk=>amphk_d,acphk=>acphk_d,apphk=>apphk_d, &
                    tempbp=>tempbp_d,temptp=>temptp_d,udx3c=>udx3c_d,udx3m=>udx3m_d,kmv=>kmv_d,kpv=>kpv_d, istat
  use local_arrays, only: vx=>vx_d,vy=>vy_d,vz=>vz_d,temp=>temp_d,pr=>pr_d,hro=>hro_d, &
                          dph=>dph_d,rux=>rux_d,ruy=>ruy_d,ruz=>ruz_d,rutemp=>rutemp_d
  use decomp_2d, only: rhs=>work1_r_d, dphhalo=>work2_r_d
  use fftw_params, only: qcap=>rbuff1_d, dq=>rbuff2_d

#else
  use param, only : nsst,alm,gam,rom,fp_kind,al,ga,ro,beta,dt,ren,lvlhalo,nxm,ap3ck,ac3ck,am3ck, &
                    ap3sk,ac3sk,am3ssk,ac3ssk,ap3ssk,am3sk,amphk,acphk,apphk,tempbp,temptp,udx3c,udx3m,kmv,kpv
  use local_arrays, only: vx,vy,vz,temp,pr,rhs,qcap,dq,hro,dph,dphhalo,rux,ruy,ruz,rutemp
  use fftw_params, only: ry1, cy1, cz1, dphc1
#endif

  use mpih
  use decomp_2d
  use solver_interfaces
  implicit none
  integer :: ns
  integer :: j,k,i
!EHP
  character(len=4)::itcount

  beta=dt/ren*real(0.5,fp_kind)

  do ns=1,nsst                                                 
!EHP
    write(itcount,'(i4)') ns
    call nvtxStartRange("Step "//itcount,ns)


!RO     Coefficients for time marching integration (alpha, gamma, rho)

    al=alm(ns)
    ga=gam(ns)
    ro=rom(ns)

    call nvtxStartRange("ExplicitTermsVX "//itcount,1)
    call ExplicitTermsVX(vx,vy,vz,temp,qcap,udx3c)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,1)
    call ExplicitTermsVX(vx_h,vy_h,vz_h,temp_h,qcap_h,udx3c_h)
    call nvtxEndRangeAsync
#endif
    call nvtxEndRange

    call nvtxStartRange("ExplicitTermsVY "//itcount,2)
    call ExplicitTermsVY(vx,vy,vz,dph,udx3m,kmv,kpv)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,2)
    call ExplicitTermsVY(vx_h,vy_h,vz_h,dph_h,udx3m_h,kmv_h,kpv_h)
    call nvtxEndRangeAsync
#endif
    call nvtxEndRange

    call nvtxStartRange("ExplicitTermsVZ "//itcount,3)
    call ExplicitTermsVZ(vx,vy,vz,dq,udx3m,kmv,kpv)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,3)
    call ExplicitTermsVZ(vx_h,vy_h,vz_h,dq_h,udx3m_h,kmv_h,kpv_h)
    call nvtxEndRangeAsync
#endif
    call nvtxEndRange

    call nvtxStartRange("ExplicitTermsTemp "//itcount,4)
    call ExplicitTermsTemp(vx,vy,vz,temp,hro,udx3c)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,4)
    call ExplicitTermsTemp(vx_h,vy_h,vz_h,temp_h,hro_h,udx3c_h)
    call nvtxEndRangeAsync
#endif
    call nvtxEndRange


    call nvtxStartRange("ImplicitAndUpdateVY "//itcount,6)
    call ImplicitAndUpdateVY(vy,pr,rhs,ruy,dph,am3sk,ac3sk,ap3sk,kmv,kpv)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,6)
    call ImplicitAndUpdateVY(vy_h,pr_h,rhs_h,ruy_h,dph_h,am3sk_h,ac3sk_h,ap3sk_h,kmv_h,kpv_h)
    call nvtxEndRangeAsync
#endif
    call nvtxEndRange


    call nvtxStartRange("ImplicitAndUpdateVZ "//itcount,7)
    call ImplicitAndUpdateVZ(vz,pr,rhs,ruz,dq,am3sk,ac3sk,ap3sk,kmv,kpv)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,7)
    call ImplicitAndUpdateVZ(vz_h,pr_h,rhs_h,ruz_h,dq_h,am3sk_h,ac3sk_h,ap3sk_h,kmv_h,kpv_h)
    call nvtxEndRangeAsync
    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange


    call nvtxStartRange("ImplicitAndUpdateVX "//itcount,5)
    call ImplicitAndUpdateVX(vx,rhs,rux,qcap,pr,am3ck,ac3ck,ap3ck,udx3c,am3ssk,ap3ssk)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR Start asyncronous exchange of interpencil VY/VZ CPU/GPU interface and boundary data 
    call update_interface(vy_h, vy)
    call update_interface(vz_h, vz)

    call transfer_halo_d2h(vy_h, vy)
    call transfer_halo_d2h(vz_h, vz)

    call nvtxStartRangeAsync("CPU"//itcount,5)
    call ImplicitAndUpdateVX(vx_h,rhs_h,rux_h,qcap_h,pr_h,am3ck_h,ac3ck_h,ap3ck_h,udx3c_h,am3ssk_h,ap3ssk_h)
    call nvtxEndRangeAsync
    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange


    call nvtxStartRange("update_halovyx "//itcount,9)

#if !defined(USE_HYBRID)
    call update_halo(vy,lvlhalo)
    call update_halo(vz,lvlhalo)
#elif defined(USE_CUDA) && defined(USE_HYBRID)

    !JR Start asyncronous exchange of interpencil VY/VZ CPU/GPU interface data 
    call update_interface(vx_h, vx)

    !JR Update halo on CPU, transfer halos to GPU
    call nvtxStartRangeAsync("CPU"//itcount,9)
    call update_halo(vy_h,lvlhalo)
    call transfer_halo_h2d(vy, vy_h)
    call update_halo(vz_h,lvlhalo)
    call transfer_halo_h2d(vz, vz_h)
    call nvtxEndRangeAsync


    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange

    call nvtxStartRange("CalcLocalDivergence "//itcount,12)
    call CalcLocalDivergence(vx,vy,vz,dph,udx3m)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    call nvtxStartRangeAsync("CPU"//itcount,12)
    call CalcLocalDivergence(vx_h,vy_h,vz_h,dph_h,udx3m_h)
    call nvtxEndRangeAsync

    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange

    call nvtxStartRange("ImplicitAndUpdateTemp "//itcount,8)
    call ImplicitAndUpdateTemp(temp,hro,rhs,rutemp,tempbp,temptp,am3ck,ac3ck,ap3ck,am3ssk,ac3ssk,ap3ssk)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR Before poisson solve, asynchronsouly transfer dph data from CPU subdomain to GPU
    istat = cudaMemcpyAsync(dph(1, xstart(2), xstart_cpu(3)), dph_h(1, xstart(2), xstart_cpu(3)), nxm * nym * xsize_cpu(3), stream = a2a_h2d)

    call nvtxStartRangeAsync("CPU"//itcount,8)
    call ImplicitAndUpdateTemp(temp_h,hro_h,rhs_h,rutemp_h,tempbp_h,temptp_h,am3ck_h,ac3ck_h,ap3ck_h,am3ssk_h,ac3ssk_h,ap3ssk_h)
    call nvtxEndRangeAsync
    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange

#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR Before poisson solve, set default stream back to true default stream. Seems to be needed for synchronization.
    istat = cudaforSetDefaultStream(0)
#endif

    call nvtxStartRange("SolvePressureCorrection "//itcount,11)
#ifdef USE_CUDA
    !JR A note on arguments, here we are reusing buffers for FFTs to save memory
    call SolvePressureCorrection(dq, qcap, dq, qcap)
#else
    call SolvePressureCorrection(ry1, cy1, cz1, dphc1)
#endif

    call nvtxEndRange

#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR After poisson solve, set default stream back to asynchronous compute stream
    istat = cudaforSetDefaultStream(a2a_comp)
    istat = cudaDeviceSynchronize
#endif

!EP this copy can be avoided by changing transpose_x_to_y_real and
!transpose_y_to_x_real so these routines can handle arrays with
!halo. This copy is a defacto array temporary. Using inferred size
!arrays in the transpose calls results in 5 more of these, and more
!memory usage.  Time spent on this copy is 0.1% for 65^3 grid.
    call nvtxStartRange("dph_to_dphhalo "//itcount,9)
#ifdef USE_CUDA
    call copy_dph_to_dphhalo(dph, dphhalo)

#ifdef USE_HYBRID
!JR Transfer bulk dphhalo data to host from dph on device during dph to dphhalo copy operation
    call transfer_dphhalo_bulk_d2h(dphhalo_h, dph)
#endif

    call update_dphhalo(dph, dphhalo)

#ifdef USE_HYBRID
    istat = cudaDeviceSynchronize
#endif

#else
    call nvtxStartRange("COPYdph "//itcount,9)
    do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
        do k=1,nxm
          dphhalo(k,j,i) = dph(k,j,i)
        enddo
      enddo
    enddo
    call nvtxEndRange

    call nvtxStartRange("update_halo "//itcount,13)
    call update_halo(dphhalo,lvlhalo)
    call nvtxEndRange
#endif
    call nvtxEndRange


    call nvtxStartRange("CorrectVelocity "//itcount,14)
    call CorrectVelocity(vx,vy,vz,dphhalo,udx3c,kmv)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR Complete transfer of dphhalo data (halos only) to CPU
    call transfer_dphhalo_halo_d2h(dphhalo_h, dphhalo)
    istat = cudaStreamSynchronize(a2a_d2h)

    !JR Transfer temperature interface data
    call update_interface(temp_h, temp)
    call transfer_halo_d2h(temp_h, temp)

    call nvtxStartRangeAsync("CPU"//itcount,14)
    call CorrectVelocity(vx_h,vy_h,vz_h,dphhalo_h,udx3c_h,kmv_h)
    call nvtxEndRangeAsync
    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange


    call nvtxStartRange("CorrectPressure "//itcount,15)
    call CorrectPressure(pr,dphhalo,amphk,acphk,apphk,kmv,kpv)
#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR Start asyncronous exchange of interpencil VX/VY/VZ CPU/GPU interface and boundary data 
    call update_interface(vx_h, vx)
    call update_interface(vy_h, vy)
    call update_interface(vz_h, vz)

    call transfer_halo_d2h(vx_h, vx)
    call transfer_halo_d2h(vy_h, vy)
    call transfer_halo_d2h(vz_h, vz)

    call nvtxStartRangeAsync("CPU"//itcount,15)
    call CorrectPressure(pr_h,dphhalo_h,amphk_h,acphk_h,apphk_h,kmv_h,kpv_h)
    call nvtxEndRangeAsync
    istat = cudaDeviceSynchronize
#endif
    call nvtxEndRange

#if defined(USE_CUDA) && defined(USE_HYBRID)
    !JR Start asyncronous exchange of interpencil Pr CPU/GPU interface data 
    call update_interface(pr_h, pr)

    istat = cudaDeviceSynchronize
#endif


    call nvtxStartRange("update_allhalos "//itcount,16)
#if !defined(USE_HYBRID)
    call update_halo(vx,lvlhalo)
    call update_halo(vy,lvlhalo)
    call update_halo(vz,lvlhalo)
    call update_halo(temp,lvlhalo)
    call update_halo(pr,lvlhalo)

#elif defined(USE_CUDA) && defined(USE_HYBRID)
    istat = cudaStreamSynchronize(a2a_comp)

    !JR Start asyncronous exchange of interpencil Pr CPU/GPU boundary data 
    call transfer_halo_d2h(pr_h, pr)

    call nvtxStartRangeAsync("CPU"//itcount,16)

    !JR Update halos on CPU, transfer updated halos to GPU
    call update_halo(vx_h,lvlhalo)
    call transfer_halo_h2d(vx, vx_h)

    call update_halo(vy_h,lvlhalo)
    call transfer_halo_h2d(vy, vy_h)

    call update_halo(vz_h,lvlhalo)
    call transfer_halo_h2d(vz, vz_h)

    call update_halo(temp_h,lvlhalo)

    istat = cudaDeviceSynchronize ! JR Wait for D2H transfer of Pr boundary data to complete
    call transfer_halo_h2d(temp, temp_h)

    call update_halo(pr_h,lvlhalo)
    call transfer_halo_h2d(pr, pr_h)
    call nvtxEndRangeAsync

    istat = cudaDeviceSynchronize
#endif

    call nvtxEndRange

    call nvtxEndRange

    enddo
    
  return                                                            
end                                                               
