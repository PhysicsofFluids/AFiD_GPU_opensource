!***********************************************************      
! This module declares the  global variables used by the solver.
! It also contains the interface to the NVTX instrumentation library 
! and a precision module to control use of single or double precision.
! A fp_kind type is defined, real data should be declared as real(fp_kind) 

! ---------
! precision
! ---------

module precision
  integer, parameter :: singlePrecision = kind(0.0)
  integer, parameter :: doublePrecision = kind(0.0d0)

#ifdef SINGLE
  integer, parameter :: fp_kind = singlePrecision
#else
  integer, parameter :: fp_kind = doublePrecision
#endif
end module precision

! ----
! nvtx
! ----

module nvtx
  use iso_c_binding
#ifdef USE_CUDA
  use cudafor
#endif
  implicit none

#ifdef USE_NVTX
  integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00',Z'00ff00ff',Z'0000ffff', &
       Z'00ff0000', Z'00ffffff']
  character(len=256),private :: tempName

  type, bind(C):: nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
     type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes
  
  interface nvtxRangePush
     ! push range with custom label and standard color
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR,len=*) :: name
     end subroutine nvtxRangePushA
     
     ! push range with custom label and custom color
     subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import:: nvtxEventAttributes
       type(nvtxEventAttributes):: event
     end subroutine nvtxRangePushEx
  end interface nvtxRangePush
  
  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop
#endif
  
contains
  
  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef USE_NVTX
    type(nvtxEventAttributes):: event
#ifdef USE_CUDA
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
#endif
  end subroutine nvtxStartRange

  subroutine nvtxStartRangeAsync(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef USE_NVTX
    type(nvtxEventAttributes):: event

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
#endif
  end subroutine nvtxStartRangeAsync

  
  subroutine nvtxEndRange
#ifdef USE_NVTX
#ifdef USE_CUDA
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif
    call nvtxRangePop
#endif
  end subroutine nvtxEndRange

  subroutine nvtxEndRangeAsync
#ifdef USE_NVTX
    call nvtxRangePop
#endif
  end subroutine nvtxEndRangeAsync

  
end module nvtx




! Declaration of global variables
!***********************************************************      
module param
  use precision
#ifdef USE_CUDA
  use cudafor
!  use cusparse
#endif
  implicit none
  !==========================================================			
  !       read from input file bou.in
  !==========================================================
  integer   :: nx, ny, nz
  integer   :: nsst, nread, ntst, ireset
  real(fp_kind)      :: walltimemax,tout,tmax
  real(fp_kind)      :: alx3,str3
  integer   :: istr3
  real(fp_kind)      :: ylen,zlen
  real(fp_kind)      :: ray,pra,dt,resid
  integer   :: inslws,inslwn
  integer   :: starea,tsta
  real(fp_kind)      :: dtmin,dtmax,limitCFL,limitVel
  integer   :: nson,idtv
  real(fp_kind)   :: gc_split_ratio ! GPU-CPU hybrid split ratio
  !=================================================
  !       end of input file
  !=================================================
  real(fp_kind) :: time
  !******* Grid parameters**************************
  real(fp_kind) :: dx,dy,dz,dxq,dyq,dzq
#ifdef USE_CUDA
  real(fp_kind), device :: dx_d,dy_d,dz_d,dxq_d,dyq_d,dzq_d
#endif
  !        
  real(fp_kind), allocatable, dimension(:) :: xc,xm
  real(fp_kind), allocatable, dimension(:) :: yc,ym
  real(fp_kind), allocatable, dimension(:) :: zc,zm
  real(fp_kind), allocatable, dimension(:) :: g3rc,g3rm
#ifdef USE_CUDA
  real(fp_kind), device,  allocatable, dimension(:) :: xc_d,xm_d
  real(fp_kind), device, allocatable, dimension(:) :: g3rc_d, g3rm_d
#endif
  !====================================================
  !******* QUANTITIES FOR DERIVATIVES******************
  real(fp_kind), allocatable, dimension(:) :: udx3c,udx3m
#ifdef USE_CUDA
  real(fp_kind), device, allocatable, dimension(:) :: udx3c_d,udx3m_d
#endif  
  !==========================================================
  !******* Grid indices**************************************
  integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv
#ifdef USE_CUDA
  integer, device, allocatable, dimension(:) :: kmc_d,kpc_d,kmv_d,kpv_d
#endif
  !===========================================================
  !******* Metric coefficients *******************************
  real(fp_kind), allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
  real(fp_kind), allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
  real(fp_kind), allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   
#ifdef USE_CUDA
  real(fp_kind), device, allocatable, dimension(:) :: ap3ck_d,ac3ck_d,am3ck_d
  real(fp_kind), device, allocatable, dimension(:) :: ap3sk_d,ac3sk_d,am3sk_d
  real(fp_kind), device, allocatable, dimension(:) :: ap3ssk_d,ac3ssk_d,am3ssk_d   
#endif
  !============================================================
  !******* Variables for Implicit, FFTW and Poisson solver****************
  real(fp_kind), allocatable, dimension(:) :: ak2,ap
  real(fp_kind), allocatable, dimension(:) :: ak1,ao
  real(fp_kind), allocatable, dimension(:) :: amphk,acphk,apphk
  real(fp_kind), allocatable, dimension(:) :: amkl, apkl, ackl, fkl
#ifdef USE_CUDA
  real(fp_kind), device, allocatable, dimension(:) :: ak2_d,ap_d
  real(fp_kind), device, allocatable, dimension(:) :: ak1_d,ao_d
  real(fp_kind), device, allocatable, dimension(:) :: amphk_d, acphk_d, apphk_d
  real(fp_kind), device, allocatable, dimension(:) :: amkl_d, apkl_d, ackl_d, fkl_d
#endif
  !===========================================================
  !******* Other variables ***********************************
  integer  :: nxm, nym, nzm
  real(fp_kind) :: ren, pec
  real(fp_kind) :: pi
  real(fp_kind) :: al,ga,ro
  real(fp_kind) :: beta
  real(fp_kind) :: qqmax,qqtot
  real(fp_kind) :: re
  real(fp_kind) :: tempmax,tempmin,tempm
  integer :: ntime
  integer, parameter:: ndv=3
  real(fp_kind), dimension(1:ndv) :: vmax
  real(fp_kind), dimension(1:3) :: gam,rom,alm
  real(fp_kind), allocatable, dimension(:,:) :: tempbp,temptp
#ifdef USE_CUDA
  real(fp_kind),device,  allocatable, dimension(:,:) :: tempbp_d,temptp_d
#endif
  
  logical :: dumpslabs=.false.
  logical :: statcal=.false.
  logical :: disscal=.false.
  logical :: readflow=.false.
  logical :: readstats=.false.
  logical :: ismaster=.false.
  logical :: resetlogstime=.false.
  logical :: variabletstep=.true.
  
  integer :: lvlhalo=1
  
#ifdef USE_CUDA
  integer(int_ptr_kind()):: freeMem,totalMem
  real(fp_kind):: rFreeMem,rUsedMem,rTotalMem,minUsedMem,maxUsedMem,minFreeMem,maxFreeMem
  integer:: istat
#endif
end module param

!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******
module local_arrays
#ifdef USE_CUDA
  use cudafor
#endif
  use param
  implicit none
  real(fp_kind) :: nusslw,nussup
  real(fp_kind),allocatable,dimension(:,:,:) :: vx,vy,vz
  real(fp_kind),allocatable,dimension(:,:,:) :: pr,temp,rhs
  real(fp_kind),allocatable,dimension(:,:,:) :: rux,ruy,ruz,rutemp
  real(fp_kind),allocatable,dimension(:,:,:) :: dph,qcap,dq,hro,dphhalo
#ifdef USE_CUDA
  attributes(pinned) :: vx, vy, vz, pr, temp, dph, dphhalo

  real(fp_kind),device,allocatable,dimension(:,:,:) :: vx_d,vy_d,vz_d
  real(fp_kind),device,allocatable,dimension(:,:,:) :: pr_d,temp_d,rhs_d
  real(fp_kind),device,allocatable,dimension(:,:,:) :: rux_d,ruy_d,ruz_d,rutemp_d
  real(fp_kind),device,allocatable,dimension(:,:,:) :: dph_d,qcap_d,dq_d,hro_d,dphhalo_d
#endif
end module local_arrays

!===============================================================
module stat_arrays
  use precision
  implicit none
  integer :: stat_columns
  real(fp_kind),allocatable, dimension(:,:) :: vz_m1,vz_m2,vz_m3,vz_m4
  real(fp_kind),allocatable, dimension(:,:) :: vy_m1,vy_m2,vy_m3,vy_m4
  real(fp_kind),allocatable, dimension(:,:) :: vx_m1,vx_m2,vx_m3,vx_m4
  real(fp_kind),allocatable, dimension(:,:) :: tp_m1,tp_m2,tp_m3,tp_m4,tpvx_m1
  real(fp_kind),allocatable, dimension(:,:) :: disste,dissth
#ifdef USE_CUDA
  real(fp_kind), device, allocatable, dimension(:,:) :: vz_m1_d,vz_m2_d,vz_m3_d,vz_m4_d
  real(fp_kind), device, allocatable, dimension(:,:) :: vy_m1_d,vy_m2_d,vy_m3_d,vy_m4_d
  real(fp_kind), device, allocatable, dimension(:,:) :: vx_m1_d,vx_m2_d,vx_m3_d,vx_m4_d
  real(fp_kind), device, allocatable, dimension(:,:) :: tp_m1_d,tp_m2_d,tp_m3_d,tp_m4_d,tpvx_m1_d
  real(fp_kind), device, allocatable, dimension(:,:) :: disste_d,dissth_d
#endif

  real(fp_kind):: tstat=real(0.,fp_kind)
  integer :: nstatsamples
end module stat_arrays
!=====================================================       
module stat3_param
  use precision
  implicit none
  integer :: kslab(1:9)
  real(fp_kind)    :: xslab(1:9)
end module stat3_param
!=====================================================       
module mpih
  implicit none
  include 'mpif.h'
  integer :: ierr
  integer, parameter :: master=0
  integer :: MDP = MPI_DOUBLE_PRECISION
end module mpih
!====================================================
module fftw_params
  !        use param, only: m2m,m2mh,m1m
  use precision
  use iso_c_binding
#ifdef USE_CUDA
  use cufft
#endif
  
  type, bind(C) :: fftw_iodim
     integer(C_INT) n, is, os
  end type fftw_iodim
  
  interface
     type(C_PTR) function fftw_plan_guru_dft(rank,dims, &
          howmany_rank,howmany_dims,in,out,sign,flags) &
          bind(C, name='fftw_plan_guru_dft')
       import
       integer(C_INT), value :: rank
       type(fftw_iodim), dimension(*), intent(in) :: dims
       integer(C_INT), value :: howmany_rank
       type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
       integer(C_INT), value :: sign
       integer(C_INT), value :: flags
     end function fftw_plan_guru_dft
     
     type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims, &
          howmany_rank,howmany_dims,in,out,flags) &
          bind(C, name='fftw_plan_guru_dft_r2c')
       import
       integer(C_INT), value :: rank
       type(fftw_iodim), dimension(*), intent(in) :: dims
       integer(C_INT), value :: howmany_rank
       type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
       real(C_DOUBLE), dimension(*), intent(out) :: in
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
       integer(C_INT), value :: flags
     end function fftw_plan_guru_dft_r2c
     
     type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims, &
          howmany_rank,howmany_dims,in,out,flags)  &
          bind(C, name='fftw_plan_guru_dft_c2r')
       import
       integer(C_INT), value :: rank
       type(fftw_iodim), dimension(*), intent(in) :: dims
       integer(C_INT), value :: howmany_rank
       type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
       real(C_DOUBLE), dimension(*), intent(out) :: out
       integer(C_INT), value :: flags
     end function fftw_plan_guru_dft_c2r
  end interface
  
  integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACKWARD,FFTW_ESTIMATE
  parameter (FFTW_PATIENT=32)   
  parameter (FFTW_ESTIMATE=64)   
  parameter (FFTW_FORWARD=-1)   
  parameter (FFTW_BACKWARD=1)   
  type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y 
  type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
  logical :: planned=.false.
  
  real(fp_kind),allocatable,dimension(:,:,:) :: ry1,rz1
  complex(fp_kind),allocatable,dimension(:,:,:) :: cy1,cz1,dphc1

  type(fftw_iodim),dimension(1) :: iodim
  type(fftw_iodim),dimension(2) :: iodim_howmany

#ifdef USE_CUDA
  integer :: batch
  integer :: cufft_plan_fwd_y, cufft_plan_bwd_y
  integer :: cufft_plan_fwd_z, cufft_plan_bwd_z
  real(fp_kind),device,allocatable,dimension(:,:,:) :: ry1_d,rz1_d !unused
  complex(fp_kind),device,allocatable,dimension(:,:,:) ::cy1_d,cz1_d,dphc1_d !unused
  complex(fp_kind),device,allocatable,dimension(:) :: cufft_workspace
  real(fp_kind),device,allocatable,dimension(:) :: rbuff1_d, rbuff2_d ! buffers for fft ops
#endif
end module fftw_params
