!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection,        !
!     CreateFFTTmpArrays, DestroyFFTTmpArrays             ! 
!                                                         ! 
!    PURPOSE: Compute the pressure correction by solving  !
!     a Poisson equation                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define BATCH_RND (32)
!#define NO_BATCH

      subroutine SolvePressureCorrection(ry, cy, cz, dphc)
      use, intrinsic :: iso_c_binding
      use param
      use fftw_params
      use local_arrays, only: dph
#ifdef USE_CUDA
      use local_arrays, only: dph_d
      use ep_trans
      use ep_solve
#endif
      use decomp_2d
      use decomp_2d_fft
      use mpih
      implicit none
      real(fp_kind), dimension(ph%yst(1):ph%yen(1), ph%yst(2):ph%yen(2), ph%yst(3):ph%yen(3)) :: ry
      complex(fp_kind), dimension(sp%yst(1):sp%yen(1), sp%yst(2):sp%yen(2), sp%yst(3):sp%yen(3)) :: cy
      complex(fp_kind), dimension(sp%zst(1):sp%zen(1), sp%zst(2):sp%zen(2), sp%zst(3):sp%zen(3)) :: cz
#ifndef USE_CUDA
      complex(fp_kind), dimension(sp%xst(1):sp%xen(1), sp%xst(2):sp%xen(2), sp%xst(3):sp%xen(3)) :: dphc
#else
      complex(fp_kind), dimension(sp%wst(1):sp%wen(1), sp%wst(2):sp%wen(2), sp%wst(3):sp%wen(3)) :: dphc
      attributes(device) :: ry, cy, cz, dphc
#endif

      integer :: i,j,k,info
      complex(fp_kind) :: acphT_b
      complex(fp_kind) :: appph(nxm-2)
      complex(fp_kind), dimension(nxm) :: acphT,apph,amph
      integer :: phpiv(nxm)
      integer :: nymh

#ifdef DEBUG
       print *,"-----In SolvePressure",minval(dph),maxval(dph)
#endif

      nymh=nym/2+1

#ifdef USE_CUDA
      call transpose_x_to_y(dph_d,ry,ph)
#else
      call transpose_x_to_y(dph,ry,ph)
#endif

#ifdef USE_CUDA
      call trans_xy(ry, work1_r_d, SIZE(ry,1), SIZE(ry,2), SIZE(ry,3))
#ifdef NO_BATCH
      do i=1,(ph%yen(1)-ph%yst(1)+1)*(ph%yen(3)-ph%yst(3)+1)
        istat = cufftExecD2Z(cufft_plan_fwd_y,work1_r_d((i-1)*nym+1),work1_c_d((i-1)*nymh+1))
      end do
#else
      istat = cufftExecD2Z(cufft_plan_fwd_y,work1_r_d,work1_c_d)
#endif
      call trans_yx(work1_c_d, cy, SIZE(cy,2), SIZE(cy,1), SIZE(cy,3))
#else
      call dfftw_execute_dft_r2c(fwd_guruplan_y,ry,cy)
#endif

#ifdef USE_CUDA
      call transpose_y_to_z(cy,cz,sp)
#else
      call transpose_y_to_z(cy,cz,sp)
#endif

#ifdef USE_CUDA
      call trans_yz(cz,work1_c_d,SIZE(cz,1),SIZE(cz,2),SIZE(cz,3))
#ifdef NO_BATCH
      do i=1,((sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1))
        istat = cufftExecZ2Z(cufft_plan_fwd_z,work1_c_d((i-1)*nzm+1),work1_c_d((i-1)*nzm+1),CUFFT_FORWARD)
      end do
#else
      istat = cufftExecZ2Z(cufft_plan_fwd_z,work1_c_d,work1_c_d,CUFFT_FORWARD)
#endif
      call trans_zy(work1_c_d,cz,SIZE(cz,3),SIZE(cz,1),SIZE(cz,2))
#else
      call dfftw_execute_dft(fwd_guruplan_z,cz,cz)
#endif

!EP   Normalize. FFT does not do this
#ifdef USE_CUDA
!$cuf kernel do(3) <<<*,*>>>
      do k=sp%zst(3),sp%zen(3)
        do j=sp%zst(2),sp%zen(2)
          do i=sp%zst(1),sp%zen(1)
            cz(i,j,k) = cz(i,j,k)/(nzm*nym)
          end do
        end do
      end do
#else
      cz = cz / (nzm*nym)
#endif

#ifdef USE_CUDA
      call transpose_z_to_x(cz,dphc,sp)      
#else
!!      call transpose_z_to_x(cz,dphc,sp)
      call transpose_z_to_y(cz,cy,sp)
      call transpose_y_to_x(cy,dphc,sp)
#endif

!RO   Solve the tridiagonal matrix with complex coefficients
#ifdef USE_CUDA
   !call epZgtsv_pressure_piv_dd( nxm, SIZE(dphc_d,2), SIZE(dphc_d,3), amphk_d, acphk_d, apphk_d, ak1_d, ak2_d, dphc_d, SIZE(dphc_d,1), scratch_d, sp%xst(2), sp%xst(3))
   !call epZgtsv_pressure_piv( nxm, SIZE(dphc_d,2), SIZE(dphc_d,3), amphk, acphk, apphk, ak1, ak2, dphc_d, SIZE(dphc_d,1), scratch, sp%xst(2), sp%xst(3))
   !call epZgtsv_pressure( nxm, SIZE(dphc_d,2), SIZE(dphc_d,3), amphk, acphk, apphk, ak1, ak2, dphc_d, SIZE(dphc_d,1), scratch)
   call tepZgtsv_pressure( nxm, SIZE(dphc,2), SIZE(dphc,3), amphk_d, acphk_d, apphk_d, ak1_d, ak2_d, dphc, SIZE(dphc,1), sp%wst(2), sp%wst(3),work1_c_d)
#else
!$OMP  PARALLEL DO COLLAPSE(2)                                          &
!$OMP   DEFAULT(none)                                                   &
!$OMP   SHARED(sp,nxm)                                                  &
!$OMP   SHARED(acphk,ak2,ak1,dphc,apphk,amphk)                          &
!$OMP   PRIVATE(apph,amph,acphT,acphT_b)                                &
!$OMP   PRIVATE(phpiv,info,appph)
      do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
         do k = 1,nxm
          acphT_b=real(1.0,fp_kind)/(acphk(k)-ak2(j)-ak1(i))
          dphc(k,j,i)=dphc(k,j,i)*acphT_b
          apph(k)=apphk(k)*acphT_b
          amph(k)=amphk(k)*acphT_b
          acphT(k)=real(1.0,fp_kind)
         enddo

         call zgttrf(nxm, amph(2), acphT, apph(1), appph, phpiv, info)

         call zgttrs('N',nxm,1,amph(2),acphT,apph(1),appph,phpiv,      &
                       dphc(1,j,i), nxm, info)
        enddo
      enddo
!$OMP END PARALLEL DO
#endif

#ifdef USE_CUDA
      call transpose_x_to_z(dphc,cz,sp)
#else
!      call transpose_x_to_z(dphc,cz,sp)
      call transpose_x_to_y(dphc,cy,sp)
      call transpose_y_to_z(cy,cz,sp)
#endif

#ifdef USE_CUDA
      call trans_yz(cz,work1_c_d,SIZE(cz,1),SIZE(cz,2),SIZE(cz,3))
#ifdef NO_BATCH
      do i=1,((sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1))
        istat = cufftExecZ2Z(cufft_plan_bwd_z,work1_c_d((i-1)*nzm+1),work1_c_d((i-1)*nzm+1),CUFFT_INVERSE)
      end do
#else
      istat = cufftExecZ2Z(cufft_plan_bwd_z,work1_c_d,work1_c_d,CUFFT_INVERSE)
#endif
      call trans_zy(work1_c_d,cz,SIZE(cz,3),SIZE(cz,1),SIZE(cz,2))
#else
      call dfftw_execute_dft(bwd_guruplan_z,cz,cz)
#endif

#ifdef USE_CUDA
      call transpose_z_to_y(cz,cy,sp)
#else
      call transpose_z_to_y(cz,cy,sp)
#endif

#ifdef USE_CUDA
      call trans_yx(cy, work1_c_d, SIZE(cy,1), SIZE(cy,2),SIZE(cy,3))

#ifdef NO_BATCH
      do i=1,(sp%yen(1)-sp%yst(1)+1)*(sp%yen(3)-sp%yst(3)+1)
        istat = cufftExecZ2D(cufft_plan_bwd_y,work1_c_d((i-1)*nymh+1),work1_r_d((i-1)*nym+1))
      end do
#else
      istat = cufftExecZ2D(cufft_plan_bwd_y,work1_c_d,work1_r_d)
#endif

      call trans_xy(work1_r_d, ry, SIZE(ry,2), SIZE(ry,1), SIZE(ry,3))
#else
      call dfftw_execute_dft_c2r(bwd_guruplan_y,cy,ry)
#endif

#ifdef USE_CUDA
      call transpose_y_to_x(ry,dph_d,ph)
#else
      call transpose_y_to_x(ry,dph,ph)
#endif

#ifdef DEBUG
       print *,"In SolvePressure",minval(dph),maxval(dph)
#endif

      return
      end

!======================================================================

      subroutine CreateFFTTmpArrays
      use fftw_params
      use decomp_2d
      use decomp_2d_fft
      use param, only: nx
      implicit none

      integer :: nd_exp, nd_ry1, nd_cy1, nd_cz1, nd_dphc1

      allocate(ry1(ph%yst(1):ph%yen(1),                                 &
     &             ph%yst(2):ph%yen(2),                                 &
     &             ph%yst(3):ph%yen(3)))
!      allocate(rz1(ph%zst(1):ph%zen(1),                                 &
!     &             ph%zst(2):ph%zen(2),                                 &
!     &             ph%zst(3):ph%zen(3)))
      allocate(cy1(sp%yst(1):sp%yen(1),                                 &
     &             sp%yst(2):sp%yen(2),                                 &
     &             sp%yst(3):sp%yen(3)))
      allocate(cz1(sp%zst(1):sp%zen(1),                                 &
     &             sp%zst(2):sp%zen(2),                                 &
     &             sp%zst(3):sp%zen(3)))
      allocate(dphc1(sp%xst(1):sp%xen(1),                                &
     &             sp%xst(2):sp%xen(2),                                 &
     &             sp%xst(3):sp%xen(3)))

#ifdef USE_CUDA
     ! JR Replacing old FFT temporary buffers. Leaving this in place for future use if necessary.
     ! allocate(ry1_d(ph%yst(1):ph%yen(1),                                 &
     !&               ph%yst(2):ph%yen(2),                                 &
     !&               ph%yst(3):ph%yen(3)))
     ! allocate(cy1_d(sp%yst(1):sp%yen(1),                                 &
     !&               sp%yst(2):sp%yen(2),                                 &
     !&               sp%yst(3):sp%yen(3)))
     ! allocate(cz1_d(sp%zst(1):sp%zen(1),                                 &
     !&             sp%zst(2):sp%zen(2),                                 &
     !&             sp%zst(3):sp%zen(3)))
     ! allocate(dphc1_d(sp%wst(1):sp%wen(1),                                &
     !&             sp%wst(2):sp%wen(2),                                 &
     !&             sp%wst(3):sp%wen(3)))

     ! JR Allocating buffers for fft ops. Size to be reused for explicit terms. 
     nd_exp = nx * xsize(2) * xsize(3)  ! required number of doubles for explicit terms
     nd_ry1 = ph%ysz(1) * ph%ysz(2) * ph%ysz(3)   ! required number of doubles for ry1
     nd_cy1 = 2*sp%ysz(1) * sp%ysz(2) * sp%ysz(3)   ! required number of doubles for ry1
     nd_cz1 = 2*sp%zsz(1) * sp%zsz(2) * sp%zsz(3)   ! required number of doubles for rz1
     nd_dphc1 = 2*sp%wsz(1) * sp%wsz(2) * sp%wsz(3)   ! required number of doubles for rz1

     allocate(rbuff1_d(max(nd_exp, max(nd_cy1, nd_dphc1))))  !sized to max size of exp term, cy1, or dphc1
     allocate(rbuff2_d(max(nd_exp, max(nd_cz1, nd_ry1))))  !sized to max size of exp term, cz1, or ry1


#endif

      return
      end

!======================================================================

      subroutine DestroyFFTTmpArrays
      use fftw_params
      implicit none

      if(allocated(dphc1)) deallocate(dphc1)
!      if(allocated(rz1)) deallocate(rz1)
      if(allocated(cz1)) deallocate(cz1)
      if(allocated(ry1)) deallocate(ry1)
      if(allocated(cy1)) deallocate(cy1)

      return
      end

!======================================================================

      subroutine InitPressureSolverPlans
      use param
      use fftw_params
      use decomp_2d_fft
      use mpih
#ifdef USE_CUDA
      integer(int_ptr_kind()) :: worksize,max_worksize

      max_worksize = 0

#ifdef NO_BATCH
        batch = 1
#else
        batch = (sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
#endif
!        istat = cufftPlan1D(cufft_plan_fwd_z,nzm,CUFFT_Z2Z,batch)
        istat = cufftCreate( cufft_plan_fwd_z )
        istat = cufftSetAutoAllocation( cufft_plan_fwd_z, 0 )
        istat = cufftMakePlanMany(cufft_plan_fwd_z,1,nzm,0,1,nzm,0,1,nzm,CUFFT_Z2Z,batch,worksize)
        max_worksize = max(worksize,max_worksize)

!        istat = cufftPlan1D(cufft_plan_bwd_z,nzm,CUFFT_Z2Z,batch)        
        istat = cufftCreate( cufft_plan_bwd_z )
        istat = cufftSetAutoAllocation( cufft_plan_bwd_z, 0 )
        istat = cufftMakePlanMany(cufft_plan_bwd_z,1,nzm,0,1,nzm,0,1,nzm,CUFFT_Z2Z,batch,worksize)
        max_worksize = max(worksize,max_worksize)

#ifdef NO_BATCH
        batch = 1
#else
        batch = (ph%yen(1)-ph%yst(1)+1)*(ph%yen(3)-ph%yst(3)+1)
        batch = ceiling(real(batch)/BATCH_RND)*BATCH_RND
#endif
        !istat = cufftPlan1D(cufft_plan_fwd_y,nym,CUFFT_D2Z,batch)
        istat = cufftCreate( cufft_plan_fwd_y )
        istat = cufftSetAutoAllocation( cufft_plan_fwd_y, 0 )
        istat = cufftMakePlanMany(cufft_plan_fwd_y,1,nym,0,1,nym,0,1,nym/2+1,CUFFT_D2Z,batch,worksize)
        max_worksize = max(worksize,max_worksize)

#ifdef NO_BATCH
        batch = 1
#else
        batch = (sp%yen(1)-sp%yst(1)+1)*(sp%yen(3)-sp%yst(3)+1)
        batch = ceiling(real(batch)/BATCH_RND)*BATCH_RND
#endif
        !istat = cufftPlan1D(cufft_plan_bwd_y,nym,CUFFT_Z2D,batch)
        istat = cufftCreate( cufft_plan_bwd_y )
        istat = cufftSetAutoAllocation( cufft_plan_bwd_y, 0 )
        istat = cufftMakePlanMany(cufft_plan_bwd_y,1,nym,0,1,nym/2+1,0,1,nym,CUFFT_Z2D,batch,worksize)
        max_worksize = max(worksize,max_worksize)

        allocate(cufft_workspace(max_worksize/16))

        istat = cufftSetWorkArea( cufft_plan_fwd_z, cufft_workspace )
        istat = cufftSetWorkArea( cufft_plan_bwd_z, cufft_workspace )
        istat = cufftSetWorkArea( cufft_plan_fwd_y, cufft_workspace )
        istat = cufftSetWorkArea( cufft_plan_bwd_y, cufft_workspace )

#else
        iodim(1)%n=nzm
        iodim(1)%is=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
        iodim(1)%os=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
        iodim_howmany(1)%n=(sp%zen(1)-sp%zst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(sp%zen(2)-sp%zst(2)+1)
        iodim_howmany(2)%is=(sp%zen(1)-sp%zst(1)+1)
        iodim_howmany(2)%os=(sp%zen(1)-sp%zst(1)+1)
        fwd_guruplan_z=fftw_plan_guru_dft(1,iodim,                      &
     &    2,iodim_howmany,cz1,cz1,                                      &
     &    FFTW_FORWARD,FFTW_ESTIMATE)
        iodim(1)%n=nzm
        bwd_guruplan_z=fftw_plan_guru_dft(1,iodim,                      &
     &    2,iodim_howmany,cz1,cz1,                                      &
     &    FFTW_BACKWARD,FFTW_ESTIMATE)

        if (.not.c_associated(bwd_guruplan_z)) then
          if (ismaster) print*,'Failed to create guru plan. You should'
          if (ismaster) print*,'link with FFTW3 before MKL'
          if (ismaster) print*,'Please check linking order.'
          call MPI_Abort(MPI_COMM_WORLD,1,info)
        endif

        iodim(1)%n=nym
        iodim(1)%is=ph%yen(1)-ph%yst(1)+1
        iodim(1)%os=sp%yen(1)-sp%yst(1)+1
        iodim_howmany(1)%n=(ph%yen(1)-ph%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(ph%yen(3)-ph%yst(3)+1)
        iodim_howmany(2)%is=(ph%yen(1)-ph%yst(1)+1)                     &
     &    *(ph%yen(2)-ph%yst(2)+1)
        iodim_howmany(2)%os=(sp%yen(1)-sp%yst(1)+1)                     &
     &    *(sp%yen(2)-sp%yst(2)+1)
        fwd_guruplan_y=fftw_plan_guru_dft_r2c(1,iodim,                  &
     &    2,iodim_howmany,ry1,cy1,                                      &
     &    FFTW_ESTIMATE)

        iodim(1)%n=nym
        iodim(1)%is=sp%yen(1)-sp%yst(1)+1
        iodim(1)%os=ph%yen(1)-ph%yst(1)+1
        iodim_howmany(1)%n=(sp%yen(1)-sp%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(sp%yen(3)-sp%yst(3)+1)
        iodim_howmany(2)%is=(sp%yen(1)-sp%yst(1)+1)                     &
     &    *(sp%yen(2)-sp%yst(2)+1)
        iodim_howmany(2)%os=(ph%yen(1)-ph%yst(1)+1)                     &
     &    *(ph%yen(2)-ph%yst(2)+1)
        bwd_guruplan_y=fftw_plan_guru_dft_c2r(1,iodim,                  &
     &    2,iodim_howmany,cy1,ry1,                                      &
     &    FFTW_ESTIMATE)
#endif
        return
      end

