      program AFiD
      use nvtx
      use mpih
      use param
      use local_arrays, only: vx,vy,vz,temp,pr,dph,dphhalo,nusslw,nussup
#ifdef USE_CUDA
      use local_arrays, only: vx_d,vy_d,vz_d,temp_d,pr_d,dph_d,dphhalo_d,nusslw,nussup
#endif
      use hdf5
      use decomp_2d
      use decomp_2d_fft
      use stat_arrays, only: nstatsamples
#ifdef DEBUG
      use ep_debug
#endif
!$    use omp_lib
      implicit none
      integer :: errorcode, nthreads,version
      real(fp_kind)    :: instCFL,dmax
      real(fp_kind)    :: ti(2), tin(3), minwtdt
      real(fp_kind) :: ts
      integer :: prow=0,pcol=0
      character(100) :: arg

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
!
      call ReadInputFile

      if (command_argument_count().eq.2) then
        call get_command_argument(1,arg)
        read(arg,'(i10)') prow
        call get_command_argument(2,arg)
        read(arg,'(i10)') pcol
      endif

      call decomp_2d_init(nxm,nym,nzm,prow,pcol, &
     & (/ .false.,.true.,.true. /))

!MF  Assign GPU to MPI task is done in decomp_2d_init 

      ts=MPI_WTIME()
      tin(1) = MPI_WTIME()

      call MpiBarrier

      call HdfStart

      if (nrank.eq.master) ismaster = .true.

      if (ismaster) write(6,*) 'MPI tasks=', nproc

!$    if (ismaster) then 
!$OMP PARALLEL
!$OMP MASTER
!$        nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$        write(6,*) 'OMP threads=', nthreads
!$    end if

      if(ismaster) then
!m==========================================    
      call ResetLogs
!m====================================================
      write(6,112)ylen/alx3,zlen/alx3
  112 format(//,20x,'R A Y B E N ',//,10x, &
       '3D Cell with aspect-ratio:  D1/H = ',f5.2,' D2/H = ',f5.2)
      write(6,142) 
  142 format(//,8x,'Periodic lateral wall boundary condition')
      write(6,202) ray,pra
  202 format(/,5x,'Parameters: ',' Ra=',e10.3,' Pr= ',e10.3)
      if(variabletstep) then
         write(6,204) limitCFL
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
      else 
         write(6,205) dtmax,limitCFL
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
      endif
#ifdef USE_CUDA
#ifndef USE_HYBRID
       write(6,*) '********  CUDA version *******'
#else
       write(6,*) '********  Hybrid CUDA + CPU version *******'
       write(6,'(5X,A5,1X,F5.2,1X,A5,1X,F5.2)') 'GPU:' , GC_SPLIT_RATIO, 'CPU:', 1.d0-GC_SPLIT_RATIO
#endif
#endif
      endif

      call InitTimeMarchScheme

      call InitVariables

      call CreateGrid

      call WriteGridInfo

      if (dumpslabs) call InitializeSlabDump

!m===================================                                                      
!m===================================
!m===================================
      if(ismaster) then
      write(6,754)nx,ny,nz                                              
  754 format(/,5x,'grid resolution: ',' nx= ',i5,' ny= ',i5, &
      ' nz= ',i5/)                       
      write(6,755) 1.d0/dx,1.d0/dy,1.d0/dz,dt,ntst                  
  755 format(/,2x,' dx=',e10.3,' dy=',e10.3,' dz=',e10.3,' dt=' &
      ,e10.3,' ntst=',i7,/)
      endif

!m===================================
!m===================================     
      
      time=real(0.,fp_kind)
      if(statcal) nstatsamples = 0

      call InitPressureSolver
      call SetTempBCs

#ifdef USE_CUDA
      istat=CudaMemGetInfo(freeMem,totalMem)
      rTotalMem = totalMem/(1024.**2)
      rFreeMem =   freeMem/(1024.**2)
      rUsedMem = (totalMem-freeMem)/(1024**2)
      call MPI_ALLREDUCE(rUsedMem,maxUsedMem,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rUsedMem,minUsedMem,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rFreeMem,maxFreeMem,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rFreeMem,minFreeMem,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(nrank==0) then
        write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
        write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",maxFreeMem," - ",minFreeMem," / ",rTotalMem," MBytes"
        print *," "
      endif
#endif


      if(readflow) then

        if(ismaster) write(6,*) ' Reading initial condition from file'

        call ReadFlowField

      else

        if(ismaster) write(6,*) ' Creating initial condition'

        ntime=0                                                         
        time=real(0.,fp_kind)
        instCFL=real(0.,fp_kind)
        
        call CreateInitialConditions

      endif                                                             

!EP   Update all relevant halos
      call update_halo(vx,lvlhalo)
      call update_halo(vy,lvlhalo)
      call update_halo(vz,lvlhalo)
      call update_halo(temp,lvlhalo)
      call update_halo(pr,lvlhalo)

#if defined(USE_CUDA) && defined(USE_HYBRID)
      istat = cudaforSetDefaultStream(a2a_comp)

!JR   Copy initial GPU subdomain data from CPU to GPU
      istat = cudaMemcpy(vx_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), vx(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2) *(xsize_gpu(3)+2))
      istat = cudaMemcpy(vy_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), vy(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2) * (xsize_gpu(3)+2))
      istat = cudaMemcpy(vz_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), vz(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2) * (xsize_gpu(3)+2))
      istat = cudaMemcpy(temp_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), temp(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2) * (xsize_gpu(3)+2))
      istat = cudaMemcpy(pr_d(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), pr(1, xstart(2)-lvlhalo, xstart(3)-lvlhalo), nx * (nym + 2) * (xsize_gpu(3)+2))

#elif defined (USE_CUDA)
!  Copy the initial solution to GPU
      vx_d=vx
      vy_d=vy
      vz_d=vz
      temp_d=temp
      pr_d=pr
#endif


!     Compute the global quantities and PlateNu from the restarted field
      if(readflow) then
         call GlobalQuantities
         call CalcPlateNu
       end if

     

!EP   Check divergence. Should be reduced to machine precision after the first
!     phcalc. Here it can still be high.

         call CheckDivergence(dmax)

      if(ismaster) write(6,*)' Initial maximum divergence: ',dmax

      if(ismaster) then
        tin(2) = MPI_WTIME()
        write(6,'(a,f6.2,a)') '  Initialization Time = ', tin(2) -tin(1), ' sec.'
        print *," "
      endif


!EP   Write some values
      if(variabletstep) then
       if(ismaster) then
          write(6,"(A8,A8,A6,2x,A8,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12)") &
                "WallDt","CFL","ntime","time","vmax(1)","vmax(2)","vmax(3)","dmax","tempm","tempmax","tempmin","nuslw","nussup"
          write(6,"(F8.3,F8.3,I6,2x,F8.3,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5)") &
                minwtdt,instCFL*dt,ntime,time,vmax(1),vmax(2),vmax(3),dmax,tempm,tempmax,tempmin,nusslw,nussup
       endif 
       !if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin
      else
       if(ismaster) then
          write(6,"(A8,A8,A6,2x,A8,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12,1x,A12)") &
                "WallDt","CFL","ntime","time","vmax(1)","vmax(2)","vmax(3)","dmax","tempm","tempmax","tempmin"
          write(6,"(F8.3,F8.3,I6,2x,F8.3,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5)") &
                minwtdt,instCFL*dt,ntime,time,vmax(1),vmax(2),vmax(3),dmax,tempm,tempmax,tempmin,nusslw,nussup
       endif
       !if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin !RO Fix??
      end if

                                                                        
!  ********* starts the time dependent calculation ***
      errorcode = 0 !EP set errocode to 0 (OK)
      minwtdt = huge(real(0.0,fp_kind)) !EP initialize minimum time step walltime
      do ntime=0,ntst                                           
        call nvtxStartRange("TimeStep",mod(ntime,2)+4) !EHP

        ti(1) = MPI_WTIME()
!EP   Determine timestep size
        call nvtxStartRange("CalcMaxCFL",2) !EHP
           call CalcMaxCFL(instCFL)
        !if(mod(time,tout).lt.dt) then
        ! if(ismaster) then
        !   write(*,*) 'CalcMaxCFL=',instCFL
        ! endif
        !endif
        call nvtxEndRange !EHP

        if(variabletstep) then
          if(ntime.gt.1) then
            if(instCFL.lt.1.0d-8) then !EP prevent fp-overflow
              dt=dtmax
            else
              dt=limitCFL/instCFL
            endif
            if(dt.gt.dtmax) dt=dtmax
              !print *,"TCFL =",dt*instCFL
          else
            dt=dtmin
          endif
            if(dt.lt.dtmin) errorcode = 166
        else  
!RO    fixed time-step
          instCFL=instCFL*dt
          if(instCFL.gt.limitCFL) errorcode = 165
        endif

        call TimeMarcher

        !if(mod(time,tout).lt.dt) then
         !if(ismaster) then
          !write(6,*) ' ---------------------------------------- '
          !write(6,"(A6,F8.3,A10,I8,A6,F8.3)") ' T = ',time,' NTIME = ',ntime,' DT = ',dt
         !endif
        !endif

        time=time+dt


        if(ntime.eq.1.or.mod(time,tout).lt.dt) then

          call nvtxStartRange("GlobalQuantities",4) !EHP
          call GlobalQuantities
          call nvtxEndRange

          if(vmax(1).gt.limitVel.and.vmax(2).gt.limitVel) errorcode = 266

            call nvtxStartRange("CalcMaxCFL",3) !EHP
            call CalcMaxCFL(instCFL)
            call nvtxEndRange !EHP

            call nvtxStartRange("CheckDivergence",4) !EHP
            call CheckDivergence(dmax)
            call nvtxEndRange !EHP

            call nvtxStartRange("CalcPlateNu",5) !EHP
            call CalcPlateNu
            call nvtxEndRange !EHP



            if(time.gt.tsta) then

             call nvtxStartRange("Stats",6) !EHP
             if (statcal)  call CalcStats
             call nvtxEndRange !EHP
             call nvtxStartRange("Slab",7) !EHP
             if (dumpslabs) call SlabDumper
             call nvtxEndRange !EHP
             call nvtxStartRange("Dis",8) !EHP
             if (disscal.and.statcal) call CalcDissipationNu
             call nvtxEndRange !EHP

            endif

            if(.not.variabletstep) instCFL=instCFL*dt

            if(dmax.gt.resid) errorcode = 169

        endif

        if(time.gt.tmax) errorcode = 333

        ti(2) = MPI_WTIME()
        minwtdt = min(minwtdt,ti(2) - ti(1))
        if(mod(time,tout).lt.dt) then
          if(ismaster) then
          !write(6,*) 'Maximum divergence = ', dmax
          write(6,"(F8.3,F8.3,I6,2x,F8.3,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5,1x,ES12.5)") &
                minwtdt,instCFL*dt,ntime,time,vmax(1),vmax(2),vmax(3),dmax,tempm,tempmax,tempmin,nusslw,nussup
          !write(6,'(a,f8.3,a)') 'Minimum Iteration Time = ', minwtdt, &
          !        ' sec.'
          endif
          minwtdt = huge(0.0d0)
        endif

       if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334

      call MpiBcastInt(errorcode)

!EP   Conditional exits
      if(errorcode.ne.0) then

!EP    dt too small
        if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)

!EP   cfl too high    
        if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)
      
!EP   velocities diverged
        if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)
          
!EP   mass not conserved
        if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)

!EP   Physical time exceeded tmax, no error; normal quit
        if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)

!EP   walltime exceeded walltimemax, no error; normal quit
        if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)

        errorcode = 100 !EP already finalized
      
        exit

      endif

      call nvtxEndRange !EHP
      enddo !EP main loop

      call QuitRoutine(tin,.true.,errorcode)
      
      end                                                               

