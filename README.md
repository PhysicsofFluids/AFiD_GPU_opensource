### AFiD_GPU
GPU version of AFID code, see also https://github.com/PhysicsofFluids/AFiD

###
When using this code please cite 

X. Zhu, E. Phillips, V. Spandan, J. Donners, G. Ruetsch, J. Romero, R. Ostilla-Mónico, Y. Yang, D. Lohse, R. Verzicco, M. Fatica, R.J.A.M. Stevens,
AFiD-GPU: a versatile Navier-Stokes Solver for Wall-Bounded Turbulent Flows on GPU Clusters,
Submitted to Computer Physics Communications, see ArXiv paper https://arxiv.org/abs/1705.01423

for a description of the GPU version of the code and 

E.P. van der Poel, R. Ostilla Mónico, J. Donners, R. Verzicco, A pencil distributed finite difference code for strongly turbulent wall-bounded flows, Computers and Fluids, 116, 10-16 (2015), see http://www.sciencedirect.com/science/article/pii/S0045793015001164

for a description about the CPU version of the development of the used numerical method.


### Building
* The code requires the PGI compiler version 15.7 or higher.

* To compile the code, use the provided `Makefile` (Cartesius) as reference and  `Makefile_PizDaint` (Piz Daint). Edit the variables at the top of the file for proper linking to external dependencies. Upon successful compilation, three executables will be generated, `afid_gpu`, `afid_cpu`, and `afid_hyb`. These correspond to the GPU version of the code, CPU version of the code, and the hybrid version respectively.  

### Running
* A couple of simple test input files have been provided within the `Run_test` directory. The procedure to run is identical to the standard CPU AFID code (https://github.com/PhysicsofFluids/AFiD). Simply rename or create a link to one of the test input files with the name 'bou.in' and run the desired executable from the directory containing this input file/link. For the GPU and CPU executables, the input files and variables within are identical to those used in the standard CPU AFID code. 

* To minimize repeating code for the various versions, many of the main solver subroutine files have been split into two files, a file that contains overloaded interfaces (one for CPU data, one for GPU data) to the subroutine and another included "common" file which contains the actual source. For example, `CalcLocalDivergence.F90` contains only the subroutine interfaces, while `CalcLocalDivergence_common.F90` contains the main routine source code, with necessary preprocessor directives to control the compilation. The generic interfaces themselves, for all the overloaded subroutines following this pattern can be found within `SolverInterfaces.F90` 

### Hybrid version (under development)
* For the hybrid code, an additional input variable, `GC_SPLIT_RATIO` is required. This variable must be set so that 0.0 < GC_SPLIT_RATIO < 1.0. This variable controls the distribution of the domain between the CPU and GPU. For example, setting GC_SPLIT_RATIO to 0.8 will associate 80% of the domain to the GPU and the remaining 20% to the CPU. In all cases, the solution of the Poisson system remains fully on the GPU. See the provided input files to see where to define this new input variable.

* When starting a simulation from scratch using the built-in initial condition, the hybrid implementation can show some stability issues, and may require a smaller timestep than the GPU or CPU code, at least initially. This is potentially due to the hybrid decomposition splitting across the extruded axis where the initial flow field is constant. Rotating the initial flow field seems to help the issue in this case. 

### Cartesius (SURFsara) compilation and run instructions (last check April 17, 2017)
To compile the code on Cartesius at SARA, first login on node gcn1 (the code requires some CUDA libraries that seems to be available only there), then:

1) Load the required modules:

   `module load mpi/mvapich2-gdr/2.1-cuda70 hdf5/mpich/pgi/1.8.16 mkl cuda/7.0.28 fftw3/intel/3.3.3-avx pgi/15.7`

2) Compile with:

  `MPICH_FC=pgf90 make -f Makefile `

The Makefile will build three executables, `afid_gpu`, `afid__cpu`, and `afid_hyb` (under development)

To submit a job on Cartesius, you can use a script similar to this one. The code is expecting a GPU for each MPI task and since Cartesius has two GPUs per node, if you want to run with 48 GPUs (`-n 48`) you will request half the nodes (`-N 24`). The nvidia-smi line will boost the clocks of the GPU. The example is using the gpu_short queue, for production runs you may want to use the gpu queue.

`################################`

`#!/bin/bash`

`#SBATCH -N 16`

`#SBATCH --tasks-per-node 2`

`#SBATCH -t 05:00`

`#SBATCH -p gpu_short`

`#SBATCH -o job_32GPU-%j.outerr`

`#SBATCH -e job_32GPU-%j.outerr`

`module load pgi/16.7 cuda`

`nvidia-smi -ac 3004,875`

`srun --mpi=pmi2  -n 32 ./afid_gpu`

`################################`

### Piz Daint (CSCS) compilation and run instructions (last check April 17, 2017)

1) Load the required modules:

   ` module load pgi cudatoolkit fftw cray-hdf5-parallel/1.8.16 intel`
   ` module load pgi cudatoolkit fftw cray-hdf5-parallel/1.8.16 intel`
   
2) Compile with:

  `MPICH_FC=pgf90 make -f Makefile `

To submit a job on Piz Daint, you can use the sample script below. 

`################################`

`#!/bin/bash -l`

`#SBATCH --nodes=12`

`#SBATCH --constraint=gpu`

`#SBATCH --ntasks-per-node=1`

`#SBATCH --gres=gpu:1`

`#SBATCH --cpus-per-task=12`

`#SBATCH -t 1-00:00:00`

`#SBATCH -p normal`

`#SBATCH -J test`

`echo "On which nodes it executes"`

`echo $SLURM_JOB_NODELIST`

`export OMP_NUM_THREADS=1`

`export PMI_NO_FORK=1`

`export OMP_SCHEDULE=static`

`echo "Now run the MPI tasks..."`

`srun -n 12 ./afid_gpu`

`################################`