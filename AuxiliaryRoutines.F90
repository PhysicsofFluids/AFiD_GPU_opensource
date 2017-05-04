!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: AuxiliaryRoutines.F90                          !
!    CONTAINS: subroutines Allocate*,Destroy*             !
!                                                         ! 
!    PURPOSE: Auxiliary routines used for memory allocs   !
!     and memory freeing                                  !
!    
!    Generic interface for host and device variables      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module AuxiliaryRoutines
      use precision

      interface AllocateReal1DArray
         module procedure AllocateReal1DArray_cpu
#ifdef USE_CUDA
         module procedure AllocateReal1DArray_gpu
#endif
      end interface AllocateReal1DArray

      interface AllocateInt1DArray
         module procedure AllocateInt1DArray_cpu
#ifdef USE_CUDA
         module procedure AllocateInt1DArray_gpu
#endif
      end interface AllocateInt1DArray

      interface AllocateReal2DArray
         module procedure AllocateReal2DArray_cpu
#ifdef USE_CUDA
         module procedure AllocateReal2DArray_gpu
#endif
      end interface AllocateReal2DArray 

      interface AllocateReal3DArray
         module procedure AllocateReal3DArray_cpu
#ifdef USE_CUDA
         module procedure AllocateReal3DArray_gpu
#endif
      end interface AllocateReal3DArray 

      contains 

        subroutine AllocateReal1DArray_cpu(var,st1,en1)
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1
        real(fp_kind),allocatable,dimension(:),intent(inout) :: var
        integer :: alloc_stat, errorcode
  
        if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
        if (alloc_stat /= 0) then
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if
  
        var =real(0.0,fp_kind)
  
        return
        end subroutine AllocateReal1DArray_cpu
      
!===========================================================================

        subroutine AllocateInt1DArray_cpu(var,st1,en1)
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1
        integer,allocatable,dimension(:),intent(inout) :: var
        integer :: alloc_stat, errorcode
  
        if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if
  
        var =0
  
        return
        end subroutine AllocateInt1DArray_cpu
  
!===========================================================================

        subroutine AllocateReal2DArray_cpu(var,st1,en1,st2,en2)
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1,st2,en2
        real(fp_kind),allocatable,dimension(:,:),intent(inout) :: var
        integer :: alloc_stat, errorcode
  
        if (.not. allocated(var)) allocate(var(st1:en1,st2:en2), stat=alloc_stat)
      
        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if
  
        var =real(0.0,fp_kind)
  
        return
        end subroutine AllocateReal2DArray_cpu

!===========================================================================

        subroutine AllocateReal3DArray_cpu(var,st1,en1,st2,en2,st3,en3)
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1,st2,en2,st3,en3
        real(fp_kind),allocatable,dimension(:,:,:),intent(inout) :: var
        integer :: alloc_stat, errorcode
  
        if (.not. allocated(var)) allocate(var(st1:en1,st2:en2,st3:en3), stat=alloc_stat)
      
        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if
  
        var =real(0.0,fp_kind)
  
        return
        end subroutine AllocateReal3DArray_cpu

#ifdef USE_CUDA

        subroutine AllocateReal1DArray_gpu(var,st1,en1)
        use cudafor
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1
        real(fp_kind),device,allocatable,dimension(:),intent(inout) :: var
        integer :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)

        if (alloc_stat /= 0) then
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if

        var =real(0.0,fp_kind)

        return
        end subroutine AllocateReal1DArray_gpu

        subroutine AllocateInt1DArray_gpu(var,st1,en1)
        use cudafor
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1
        integer,device,allocatable,dimension(:),intent(inout) :: var
        integer :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)

        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if

        var =0

        return
        end subroutine AllocateInt1DArray_gpu


        subroutine AllocateReal2DArray_gpu(var,st1,en1,st2,en2)
        use cudafor
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1,st2,en2
        real(fp_kind),device,allocatable,dimension(:,:),intent(inout) :: var
        integer :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate(var(st1:en1,st2:en2), stat=alloc_stat)

        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if

        var =real(0.0,fp_kind)

        return
        end subroutine AllocateReal2DArray_gpu

        subroutine AllocateReal3DArray_gpu(var,st1,en1,st2,en2,st3,en3)
        use cudafor
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1,st2,en2,st3,en3
        real(fp_kind),device,allocatable,dimension(:,:,:),intent(inout) :: var
        integer :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate(var(st1:en1,st2:en2,st3:en3),stat=alloc_stat)

        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if

        var =real(0.0,fp_kind)

        return
        end subroutine AllocateReal3DArray_gpu
#endif
!===========================================================================

        subroutine DestroyReal1DArray(var)
        use decomp_2d
        implicit none
        real(fp_kind),allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal1DArray

!===========================================================================

        subroutine DestroyInt1DArray(var)
        use decomp_2d
        implicit none
        integer,allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyInt1DArray

!===========================================================================

        subroutine DestroyReal2DArray(var)
        use decomp_2d
        implicit none
        real(fp_kind),allocatable,dimension(:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal2DArray

!===========================================================================

        subroutine DestroyReal3DArray(var)
        use decomp_2d
        implicit none
        real(fp_kind),allocatable,dimension(:,:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal3DArray
      
      end module AuxiliaryRoutines
