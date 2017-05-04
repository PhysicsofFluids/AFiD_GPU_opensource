!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateVariables.F90                        !
!    CONTAINS: subroutine DeallocateVariables             !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DeallocateVariables
      use param
      use local_arrays
      use stat_arrays
      use AuxiliaryRoutines
      implicit none
      
      call DestroyReal1DArray(zc)
      call DestroyReal1DArray(zm)
      call DestroyReal1DArray(ak1)
      call DestroyReal1DArray(ao)

      call DestroyReal1DArray(yc)
      call DestroyReal1DArray(ym)
      call DestroyReal1DArray(ak2)
      call DestroyReal1DArray(ap)

      call DestroyReal1DArray(xc)
      call DestroyReal1DArray(xm)
      call DestroyReal1DArray(g3rc)
      call DestroyReal1DArray(g3rm)

      call DestroyReal1DArray(udx3c)
      call DestroyReal1DArray(udx3m)

      call DestroyReal1DArray(ap3ck)
      call DestroyReal1DArray(ac3ck)
      call DestroyReal1DArray(am3ck)

      call DestroyReal1DArray(ap3sk)
      call DestroyReal1DArray(ac3sk)
      call DestroyReal1DArray(am3sk)

      call DestroyReal1DArray(ap3ssk)
      call DestroyReal1DArray(ac3ssk)
      call DestroyReal1DArray(am3ssk)

      call DestroyReal1DArray(amphk)
      call DestroyReal1DArray(acphk)
      call DestroyReal1DArray(apphk)

      call DestroyReal1DArray(amkl)
      call DestroyReal1DArray(apkl)
      call DestroyReal1DArray(ackl)
      call DestroyReal1DArray(fkl)
 
      call DestroyInt1dArray(kmc)
      call DestroyInt1dArray(kpc)
      call DestroyInt1dArray(kmv)
      call DestroyInt1dArray(kpv)

      call DestroyReal2DArray(tempbp)
      call DestroyReal2DArray(temptp)

      call DestroyReal2DArray(vx_m1)
      call DestroyReal2DArray(vy_m1)
      call DestroyReal2DArray(vz_m1)
      call DestroyReal2DArray(tp_m1)

      call DestroyReal2DArray(vx_m2)
      call DestroyReal2DArray(vy_m2)
      call DestroyReal2DArray(vz_m2)
      call DestroyReal2DArray(tp_m2)

      call DestroyReal2DArray(vx_m3)
      call DestroyReal2DArray(vy_m3)
      call DestroyReal2DArray(vz_m3)
      call DestroyReal2DArray(tp_m3)

      call DestroyReal2DArray(vx_m4)
      call DestroyReal2DArray(vy_m4)
      call DestroyReal2DArray(vz_m4)
      call DestroyReal2DArray(tp_m4)

      call DestroyReal2DArray(tpvx_m1)

      call DestroyReal2DArray(disste)
      call DestroyReal2DArray(dissth)

      call DestroyReal3DArray(rhs)

      call DestroyReal3DArray(rux)
      call DestroyReal3DArray(ruy)
      call DestroyReal3DArray(ruz)
      call DestroyReal3DArray(rutemp)

!JR   Using basic deallocate for pinned data
      deallocate(vx, vy, vz, temp, pr, dph, dphhalo)

! Clean up GPU memory
#ifdef USE_CUDA
     deallocate(vx_m1_d,vx_m2_d,vx_m3_d,vx_m4_d)
     deallocate(vy_m1_d,vy_m2_d,vy_m3_d,vy_m4_d)
     deallocate(vz_m1_d,vz_m2_d,vz_m3_d,vz_m4_d)
     deallocate(tp_m1_d,tp_m2_d,tp_m3_d,tp_m4_d,tpvx_m1_d)
     deallocate(disste_d,dissth_d)

     deallocate(vx_d,vy_d,vz_d,pr_d,temp_d)
     deallocate(rux_d,ruy_d,ruz_d,hro_d,rutemp_d)
     deallocate(dph_d, amkl_d, apkl_d, ackl_d, fkl_d)
#endif

      return 
      end   

