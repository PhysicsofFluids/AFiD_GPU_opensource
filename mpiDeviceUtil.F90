#ifdef USE_CUDA
module mpiDeviceUtil
  interface
     subroutine quicksort(base, nmemb, elemsize, compar) &
          bind(C,name='qsort')
       use iso_c_binding
       implicit none
       !pgi$ ignore_tkr base,nmemb,elemsize,compar
       type(C_PTR), value :: base
       integer(C_SIZE_T), value :: nmemb, elemsize
       type(C_FUNPTR), value :: compar
     end subroutine quicksort

     integer function strcmp(a,b) bind(C,name='strcmp')
       use iso_c_binding
       implicit none
       !pgi$ ignore_tkr a,b
       type(C_PTR), value :: a, b
     end function strcmp
  end interface
contains
  subroutine assignDevice(dev)
    use mpih
    use cudafor
    implicit none
    integer :: dev
    character (len=MPI_MAX_PROCESSOR_NAME), allocatable :: hosts(:)
    character (len=MPI_MAX_PROCESSOR_NAME) :: hostname
    integer :: namelength, color, i, j
    integer :: nProcs, myrank, newComm, newRank

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    ! allocate array of hostnames
    allocate(hosts(0:nProcs-1))
  
    ! Every process collects the hostname of all the nodes
    call MPI_GET_PROCESSOR_NAME(hostname, namelength, ierr)
    hosts(myrank)=hostname(1:namelength)

    do i=0,nProcs-1
       call MPI_BCAST(hosts(i),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,i, &
            MPI_COMM_WORLD,ierr)
    end do
  
    ! sort the list of names
    call quicksort(hosts,nProcs,MPI_MAX_PROCESSOR_NAME,strcmp)

    ! assign the same color to the same node
    color=0
    do i=0,nProcs-1
       if (i > 0) then
          if ( lne(hosts(i-1),hosts(i)) ) color=color+1
       end if
       if ( leq(hostname,hosts(i)) ) exit
    end do
  
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,newComm,ierr)
    call MPI_COMM_RANK(newComm, newRank, ierr)

    dev = newRank
    ierr = cudaSetDevice(dev)
#if 0
    do i=0,nProcs-1
      if(myrank == i) then
          write(6,"(A8,I4,A8,A12,A12,I2)") "Rank: ",myrank,"Host: ",hostname(1:namelength),"Using GPU: ",dev
      endif
      do j=0,1000
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end do
    end do
#endif
    deallocate(hosts)
  end subroutine assignDevice

  ! lexical .eq.
  function leq(s1, s2) result(res)
    implicit none
    character (len=*) :: s1, s2
    logical :: res    
    res = .false.
    if (lle(s1,s2) .and. lge(s1,s2)) res = .true.
  end function leq

  ! lexical .ne.
  function lne(s1, s2) result(res)
    implicit none
    character (len=*) :: s1, s2
    logical :: res    
    res = .not. leq(s1, s2)
  end function lne
end module mpiDeviceUtil
#endif
