# Makefile for Cartesius system
# 1) first login on node gcn1
# 2) Load the required modules:
#    module load mpi/mvapich2-gdr/2.1-cuda70 hdf5/mpich/pgi/1.8.16 mkl cuda/7.0.28 fftw3/intel/3.3.3-avx pgi/15.7
# 3) Compile:
#   MPICH_FC=pgf90 make -f Makefile_cartesius
#
all: afid_gpu afid_cpu afid_hyb

FC = pgf90
LINKER= mpif90

# Edit the following variables for target system___________________
#CUDA_HOME = /usr/local/cuda-7.5
HDF5_ROOT=$(SURFSARA_HDF5_ROOT)
MPI_ROOT=/hpc/sw/pgi-15.10/linux86-64/2015/mpi/mpich

FFTW3_LIBS = -L$(SURFSARA_FFTW3_LIB) -lfftw3
HDF5_LIBS = -L$(SURFSARA_HDF5_ROOT)/lib  -lhdf5_fortran -lhdf5  -lz -ldl -lm
BLAS_LIBS= -L$(SURFSARA_MKL_LIB) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread


LDFLAGS = -lpmi2 -Mcudalib=cufft -L$(CUDA_HOME)/lib64 -lnvToolsExt  $(FFTW3_LIBS) $(BLAS_LIBS) $(HDF5_LIBS)   -lirc

# _________________________________________________________________


INCLUDE=  -I$(HDF5_ROOT)/include
INCLUDE+=  -I$(MPI_ROOT)/include 

# Flags for CPU version
afid_cpu: FFLAGS = -O3 -mp -Kieee -DPASS_ARRAYS
#-DDEBUG

# Flags for CUDA  version
afid_gpu: FFLAGS = -O3 -DUSE_CUDA -DPASS_ARRAYS -Mcuda=cc35,cuda7.0,ptxinfo,madconst -Kieee -Mlarge_arrays -DUSE_NVTX
#-DDEBUG
#-Mcuda=nofma

# Flags for hybrid version
afid_hyb: FFLAGS = -O3 -DUSE_CUDA -DPASS_ARRAYS -Mcuda=cc35,cuda7.0,ptxinfo,madconst -Kieee -Mlarge_arrays -DUSE_HYBRID -mp -DUSE_NVTX
#-DDEBUG
#-Mcuda=nofma

OBJ= param.o \
     mpiDeviceUtil.o \
     decomp_2d.o \
     tridiag.o \
     decomp_2d_fft.o \
     dph_routines.o \
     hybrid_comm_routines.o \
     AuxiliaryRoutines.o \
     CalcDissipationNu.o  \
     CalcMaxCFL.o \
     CalcPlateNu.o \
     CheckDivergence.o  \
     CreateGrid.o  \
     CreateInitialConditions.o \
     DeallocateVariables.o \
     DebugRoutines.o      \
     GlobalQuantities.o \
     HdfRoutines.o \
     HdfReadContinua.o  \
     interp.o \
     InitVariables.o \
     InitTimeMarchScheme.o \
     LocateLargeDivergence.o \
     MpiAuxRoutines.o \
     QuitRoutine.o \
     ReadFlowField.o \
     ReadInputFile.o \
     ResetLogs.o \
     SetTempBCs.o \
     SlabDumpRoutines.o \
     SolverInterfaces.o \
     SolvePressureCorrection.o \
     TimeMarcher.o \
     WriteFlowField.o \
     WriteGridInfo.o \
     InitPressureSolver.o \
     StatReadReduceWrite.o \
     StatRoutines.o \
     main.o  

# Objects which require two compilation passes for hybrid version
OBJ2= CalcLocalDivergence.o \
      CorrectPressure.o \
      CorrectVelocity.o \
      ExplicitTermsVX.o \
      ExplicitTermsVY.o \
      ExplicitTermsVZ.o \
      ExplicitTermsTemp.o \
      ImplicitAndUpdateVX.o \
      ImplicitAndUpdateVY.o \
      ImplicitAndUpdateVZ.o \
      ImplicitAndUpdateTemp.o \
      SolveImpEqnUpdate_X.o \
      SolveImpEqnUpdate_YZ.o\
      SolveImpEqnUpdate_Temp.o

# Create the subdirectory for objects
OBJDIRCPU := $(shell mkdir -p Obj_cpu)
OBJDIRGPU := $(shell mkdir -p Obj_gpu)
OBJDIRHYB := $(shell mkdir -p Obj_hyb Obj_hyb/cpu Obj_hyb/gpu)

OBJ_CPU = $(foreach O, ${OBJ}, Obj_cpu/$(O))
OBJ_CPU += $(foreach O, ${OBJ2}, Obj_cpu/$(O))
OBJ_GPU = $(foreach O, ${OBJ}, Obj_gpu/$(O))
OBJ_GPU += $(foreach O, ${OBJ2}, Obj_gpu/$(O))
OBJ_HYB = $(foreach O, ${OBJ}, Obj_hyb/$(O))
OBJ_HYB += $(foreach O, ${OBJ2}, Obj_hyb/cpu/$(O))
OBJ_HYB += $(foreach O, ${OBJ2}, Obj_hyb/gpu/$(O))


Obj_gpu/%.o: %.F90  $(AUX_DEP)
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_gpu

Obj_gpu/%.o: %.CUF 
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_gpu

Obj_cpu/%.o: %.F90  $(AUX_DEP)
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_cpu

Obj_cpu/%.o: %.CUF 
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_cpu

Obj_hyb/%.o: %.F90  $(AUX_DEP)
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_hyb

Obj_hyb/%.o: %.CUF 
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_hyb

Obj_hyb/cpu/%.o: %.F90  $(AUX_DEP)
	$(FC) -c  $(FFLAGS) $(INCLUDE)  $< -o $@ -module=Obj_hyb

Obj_hyb/gpu/%.o: %.F90  $(AUX_DEP)
	$(FC) -c  $(FFLAGS) -DUSE_GPU $(INCLUDE)  $< -o $@ -module=Obj_hyb

afid_gpu: ${OBJ_GPU}
	$(LINKER)  $(FFLAGS) $(OBJ_GPU) $(FCLIBS) $(LDFLAGS) -o afid_gpu 

afid_cpu: $(OBJ_CPU)
	$(LINKER) $(FFLAGS) $(OBJ_CPU) $(FCLIBS) $(LDFLAGS) -o afid_cpu 

afid_hyb: $(OBJ_HYB)
	$(LINKER) $(FFLAGS) $(OBJ_HYB) $(FCLIBS) $(LDFLAGS) -o afid_hyb 

clean:
	rm  afid_cpu afid_gpu afid_hyb Obj_cpu/*.o  Obj_gpu/*.o  Obj_cpu/*.mod  Obj_gpu/*.mod Obj_hyb/*.o Obj_hyb/*.mod Obj_hyb/cpu/*.o Obj_hyb/gpu/*.o
