!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         ! 
!    PURPOSE: Calculate maximum velocity and temperature. !
!     volume averaged Nusselt number and Reynolds number  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GlobalQuantities
#if defined(USE_CUDA) && defined(USE_HYBRID)
  use param, g3rc_h=>g3rc, g3rc=>g3rc_d
  use local_arrays, only: vy_h=>vy, vy=>vy_d, vx_h=>vx, vx=>vx_d, vz_h=>vz, vz=>vz_d, temp_h=>temp, temp=>temp_d
#elif defined (USE_CUDA)
  use param, g3rc=>g3rc_d
  use local_arrays, only: vy=>vy_d, vx=>vx_d, vz=>vz_d, temp=>temp_d
#else
  use param
  use local_arrays, only: vy,vx,vz,temp
#endif
#ifdef USE_HYBRID
  use decomp_2d, only: xstart=>xstart_gpu, xend=>xend_gpu, xstart_cpu, xend_cpu
#else
  use decomp_2d, only: xstart,xend
#endif
  use mpih
  use nvtx
  implicit none
  integer :: jc,kc,kp,ic
  real(fp_kind) :: anusin,vol,vxcen,fac2,tempcen
  real(fp_kind) :: vy_rms_vol,vz_rms_vol
  real(fp_kind) :: vx_rms_vol,vzvyvx_rms_vol,rradpr
  real(fp_kind) :: vmax1, vmax2, vmax3

#ifdef USE_HYBRID
  real(fp_kind) :: anusin_h, tempmax_h, tempmin_h, tempm_h
  real(fp_kind) :: vy_rms_vol_h,vz_rms_vol_h
  real(fp_kind) :: vx_rms_vol_h,vzvyvx_rms_vol_h
  real(fp_kind) :: vmax1_h, vmax2_h, vmax3_h

  real(fp_kind), device :: anusin_d, tempmax_d, tempmin_d, tempm_d
  real(fp_kind), device :: vy_rms_vol_d,vz_rms_vol_d
  real(fp_kind), device :: vx_rms_vol_d,vzvyvx_rms_vol_d
  real(fp_kind), device :: vmax1_d, vmax2_d, vmax3_d

  tempmax_d=-huge(real(0.,fp_kind))
  tempmin_d= huge(real(0.,fp_kind))
  tempm_d       =real(0.,fp_kind)
  anusin_d      =real(0.,fp_kind) 
  vx_rms_vol_d = real(0.,fp_kind)
  vy_rms_vol_d = real(0.,fp_kind)
  vz_rms_vol_d = real(0.,fp_kind)
  vzvyvx_rms_vol_d = real(0.,fp_kind)
  vmax1_d = real(0.,fp_kind)
  vmax2_d = real(0.,fp_kind)
  vmax3_d = real(0.,fp_kind)  
#endif
  
  !EP   Initialize
  tempmax=-huge(real(0.,fp_kind))
  tempmin= huge(real(0.,fp_kind))
  tempm       =real(0.,fp_kind)
  anusin      =real(0.,fp_kind) 
  vx_rms_vol = real(0.,fp_kind)
  vy_rms_vol = real(0.,fp_kind)
  vz_rms_vol = real(0.,fp_kind)
  vzvyvx_rms_vol = real(0.,fp_kind)
  vmax1 = real(0.,fp_kind)
  vmax2 = real(0.,fp_kind)
  vmax3 = real(0.,fp_kind)  
  vol  = real(1.0,fp_kind)/(alx3*dx*real(nzm,fp_kind)*real(nym,fp_kind))

#ifdef USE_CUDA
  !$cuf kernel do (3) <<<*,*>>>
#endif
  !EP   Loop over volume
  do ic=xstart(3),xend(3)
     do jc=xstart(2),xend(2)
        do kc=1,nxm
           kp = kc + 1
           fac2 = g3rc(kc)
#ifdef USE_HYBRID
! JR Reduce into device resident variables to avoid implicit sync after CUF kernel
           vmax1_d = max(vmax1_d,abs(vz(kc,jc,ic)))
           vmax2_d = max(vmax2_d,abs(vy(kc,jc,ic)))
           vmax3_d = max(vmax3_d,abs(vx(kc,jc,ic)))
           tempmax_d = max(tempmax_d,temp(kc,jc,ic))
           tempmin_d = min(tempmin_d,temp(kc,jc,ic))
#else
           vmax1 = max(vmax1,abs(vz(kc,jc,ic)))
           vmax2 = max(vmax2,abs(vy(kc,jc,ic)))
           vmax3 = max(vmax3,abs(vx(kc,jc,ic)))
           tempmax = max(tempmax,temp(kc,jc,ic))
           tempmin = min(tempmin,temp(kc,jc,ic))
#endif
           vxcen = (vx(kc,jc,ic)+vx(kp,jc,ic))*real(0.5,fp_kind)
           tempcen = (temp(kc,jc,ic)+temp(kp,jc,ic))*real(0.5,fp_kind)

#ifdef USE_HYBRID
           anusin_d=anusin_d+tempcen*vxcen*fac2
           tempm_d=tempm_d+tempcen*fac2
           vx_rms_vol_d = vx_rms_vol_d + fac2*vx(kc,jc,ic)**2
           vy_rms_vol_d = vy_rms_vol_d + fac2*vy(kc,jc,ic)**2
           vz_rms_vol_d = vz_rms_vol_d + fac2*vz(kc,jc,ic)**2
           vzvyvx_rms_vol_d = vzvyvx_rms_vol_d + fac2* &
                &    (vz(kc,jc,ic)**2+vy(kc,jc,ic)**2+vx(kc,jc,ic)**2)  
#else
           anusin=anusin+tempcen*vxcen*fac2
           tempm=tempm+tempcen*fac2
           vx_rms_vol = vx_rms_vol + fac2*vx(kc,jc,ic)**2
           vy_rms_vol = vy_rms_vol + fac2*vy(kc,jc,ic)**2
           vz_rms_vol = vz_rms_vol + fac2*vz(kc,jc,ic)**2
           vzvyvx_rms_vol = vzvyvx_rms_vol + fac2* &
                &    (vz(kc,jc,ic)**2+vy(kc,jc,ic)**2+vx(kc,jc,ic)**2)  
#endif
        enddo
     enddo
  enddo

#ifdef USE_HYBRID
  !EP   Loop over volume
  call nvtxStartRangeAsync("CPU", 4)
  !$OMP PARALLEL DO DEFAULT(none) &
  !$OMP PRIVATE(ic, jc, kc, kp, fac2, vxcen, tempcen) &
  !$OMP SHARED(xstart_cpu, xend_cpu, vx_h, vy_h, vz_h, temp_h, g3rc_h, nxm) &
  !$OMP REDUCTION(max:vmax1, vmax2, vmax3, tempmax) &
  !$OMP REDUCTION(+:anusin, tempm, vx_rms_vol, vy_rms_vol, vz_rms_vol, vzvyvx_rms_vol) &
  !$OMP REDUCTION(min:tempmin)
  do ic=xstart_cpu(3),xend_cpu(3)
     do jc=xstart_cpu(2),xend_cpu(2)
        do kc=1,nxm
           kp = kc + 1
           fac2 = g3rc_h(kc)
           vmax1 = max(vmax1,abs(vz_h(kc,jc,ic)))
           vmax2 = max(vmax2,abs(vy_h(kc,jc,ic)))
           vmax3 = max(vmax3,abs(vx_h(kc,jc,ic)))
           tempmax = max(tempmax,temp_h(kc,jc,ic))
           tempmin = min(tempmin,temp_h(kc,jc,ic))
           vxcen = (vx_h(kc,jc,ic)+vx_h(kp,jc,ic))*real(0.5,fp_kind)
           tempcen = (temp_h(kc,jc,ic)+temp_h(kp,jc,ic))*real(0.5,fp_kind)
           anusin=anusin+tempcen*vxcen*fac2
           tempm=tempm+tempcen*fac2
           vx_rms_vol = vx_rms_vol + fac2*vx_h(kc,jc,ic)**2
           vy_rms_vol = vy_rms_vol + fac2*vy_h(kc,jc,ic)**2
           vz_rms_vol = vz_rms_vol + fac2*vz_h(kc,jc,ic)**2
           vzvyvx_rms_vol = vzvyvx_rms_vol + fac2* &
                &    (vz_h(kc,jc,ic)**2+vy_h(kc,jc,ic)**2+vx_h(kc,jc,ic)**2)  
        enddo
     enddo
  enddo
  call nvtxEndRangeAsync

!JR Reduce host and device variables
  vmax1_h = vmax1_d
  vmax1 = max(vmax1_h, vmax1) 
  vmax2_h = vmax2_d
  vmax2 = max(vmax2_h, vmax2) 
  vmax3_h = vmax3_d
  vmax3 = max(vmax3_h, vmax3) 
  tempmax_h = tempmax
  tempmax = max(tempmax_h, tempmax)
  tempmin_h = tempmin
  tempmin = min(tempmin_h, tempmin)

  anusin_h = anusin_d 
  anusin = anusin + anusin_h
  tempm_h = tempm_d 
  tempm = tempm + tempm_h
  vx_rms_vol_h = vx_rms_vol_d 
  vx_rms_vol = vx_rms_vol + vx_rms_vol_h
  vy_rms_vol_h = vy_rms_vol_d 
  vy_rms_vol = vy_rms_vol + vy_rms_vol_h
  vz_rms_vol_h = vz_rms_vol_d 
  vz_rms_vol = vz_rms_vol + vz_rms_vol_h
  vzvyvx_rms_vol_h = vzvyvx_rms_vol_d 
  vzvyvx_rms_vol = vzvyvx_rms_vol + vzvyvx_rms_vol_h
#endif


  vmax(1) = vmax1
  vmax(2) = vmax2
  vmax(3) = vmax3
  
  !EP   Reduce
  
  call MpiSumRealScalar(tempm)
  call MpiSumRealScalar(anusin)
  call MpiSumRealScalar(vx_rms_vol)
  call MpiSumRealScalar(vy_rms_vol)
  call MpiSumRealScalar(vz_rms_vol)
  call MpiSumRealScalar(vzvyvx_rms_vol)
  call MpiMinRealScalar(tempmin)
  call MpiMaxRealScalar(tempmax)
  call MpiMaxRealScalar(vmax(1))
  call MpiMaxRealScalar(vmax(2))
  call MpiMaxRealScalar(vmax(3))
  
  !EP   Write logs
  if(ismaster) then
     
     anusin=real(1.0,fp_kind) + sqrt(pra*ray)*anusin*vol
     
     open(95,file='nu_vol.out',status='unknown',access='sequential', &
          position='append')
     write(95,*) time, anusin
     close(95)
     
     rradpr=sqrt(ray/pra)
     tempm=tempm*vol
     vx_rms_vol=sqrt(vx_rms_vol*vol)*rradpr
     vy_rms_vol=sqrt(vy_rms_vol*vol)*rradpr
     vz_rms_vol=sqrt(vz_rms_vol*vol)*rradpr
     vzvyvx_rms_vol=sqrt(vzvyvx_rms_vol*vol)*rradpr
     
     open(94,file='rms_vel.out',status='unknown',position='append', &
          access='sequential')
     write(94,*) time,vz_rms_vol,vy_rms_vol,vx_rms_vol, &
          & vzvyvx_rms_vol
     close(94)
     
  endif
  
  return   
end subroutine GlobalQuantities
