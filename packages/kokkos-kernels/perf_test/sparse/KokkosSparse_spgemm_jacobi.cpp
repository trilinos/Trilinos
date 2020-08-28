/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Seher Acer (sacer@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <iostream>
#include "KokkosKernels_config.h"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_run_spgemm_jacobi.hpp"

void print_options(){
  std::cerr << "Options\n" << std::endl;
  std::cerr << "\t[Required] INPUT MATRIX: '--amtx [left_hand_side.mtx]' -- for C=AxA" << std::endl;
  std::cerr << "\t[Optional] BACKEND: '--threads [numThreads]' | '--openmp [numThreads]' | '--cuda [cudaDeviceIndex]' --> if none are specified, Serial is used (if enabled)" << std::endl;
  std::cerr << "\t[Optional] --bmtx [righ_hand_side.mtx]' for C = AxB" << std::endl;
  std::cerr << "\t[Optional] OUTPUT MATRICES: '--cmtx [output_matrix.mtx]' --> to write output C=AxB"  << std::endl;
  std::cerr << "\t[Optional] --DENSEACCMAX: on CPUs default algorithm may choose to use dense accumulators. This parameter defaults to 250k, which is max k value to choose dense accumulators. This can be increased with more memory bandwidth." << std::endl;
  std::cerr << "\t[Optional] The memory space used for each matrix: '--memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, check generate_makefile.sh" << std::endl;
  std::cerr << "\t[Optional] Loop scheduling: '--dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
  std::cerr << "\t[Optional] Verbose Output: '--verbose'" << std::endl;
}

static char* getNextArg(int& i, int argc, char** argv)
{
  i++;
  if(i >= argc)
    {
      std::cerr << "Error: expected additional command-line argument!\n";
      exit(1);
    }
  return argv[i];
}

int parse_inputs (KokkosKernels::Experiment::Parameters &params, int argc, char **argv){
  
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "--threads" ) ) {
      params.use_threads = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--openmp" ) ) {
      params.use_openmp = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--cuda" ) ) {
      params.use_cuda = atoi(getNextArg(i, argc, argv)) + 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--repeat" ) ) {
      params.repeat = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--hashscale" ) ) {
      params.minhashscale = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--chunksize" ) ) {
      params.chunk_size = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--teamsize" ) ) {
      params.team_size = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--vectorsize" ) ) {
      params.vector_size = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--compression2step" ) ) {
      params.compression2step =  true ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--shmem" ) ) {
      params.shmemsize = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--memspaces" ) ) {
      int memspaces = atoi(getNextArg(i, argc, argv));
      int memspaceinfo = memspaces;
      std::cout << "memspaceinfo:" << memspaceinfo << std::endl;
      if (memspaceinfo & 1){
        params.a_mem_space = 1;
        std::cout << "Using HBM for A" << std::endl;
      }
      else {
        params.a_mem_space = 0;
        std::cout << "Using DDR4 for A" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
      if (memspaceinfo & 1){
        params.b_mem_space = 1;
        std::cout << "Using HBM for B" << std::endl;
      }
      else {
        params.b_mem_space = 0;
        std::cout << "Using DDR4 for B" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
      if (memspaceinfo & 1){
        params.c_mem_space = 1;
        std::cout << "Using HBM for C" << std::endl;
      }
      else {
        params.c_mem_space = 0;
        std::cout << "Using DDR4 for C" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
      if (memspaceinfo & 1){
        params.work_mem_space = 1;
        std::cout << "Using HBM for work memory space" << std::endl;
      }
      else {
        params.work_mem_space = 0;
        std::cout << "Using DDR4 for work memory space" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--CRWC" ) ) {
      params.calculate_read_write_cost = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--CIF" ) ) {
      params.coloring_input_file = getNextArg(i, argc, argv);
    }
    else if ( 0 == strcasecmp( argv[i] , "--COF" ) ) {
      params.coloring_output_file = getNextArg(i, argc, argv);
    }
    else if ( 0 == strcasecmp( argv[i] , "--CCO" ) ) {
      //if 0.85 set, if compression does not reduce flops by at least 15% symbolic will run on original matrix.
      //otherwise, it will compress the graph and run symbolic on compressed one.
      params.compression_cut_off = atof(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--FLHCO" ) ) {
      //if linear probing is used as hash, what is the max occupancy percantage we allow in the hash.
      params.first_level_hash_cut_off = atof(getNextArg(i, argc, argv));
    }

    else if ( 0 == strcasecmp( argv[i] , "--flop" ) ) {
      //print flop statistics. only for the first repeat.
      params.calculate_read_write_cost = 1;
    }

    else if ( 0 == strcasecmp( argv[i] , "--checkoutput" ) ) {
      //check correctness
      params.check_output = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--amtx" ) ) {
      //A at C=AxB
      params.a_mtx_bin_file = getNextArg(i, argc, argv);
    }

    else if ( 0 == strcasecmp( argv[i] , "--bmtx" ) ) {
      //B at C=AxB.
      //if not provided, C = AxA will be performed.
      params.b_mtx_bin_file = getNextArg(i, argc, argv);
    }
    else if ( 0 == strcasecmp( argv[i] , "--cmtx" ) ) {
      //if provided, C will be written to given file.
      //has to have ".bin", or ".crs" extension.
      params.c_mtx_bin_file = getNextArg(i, argc, argv);
    }
    else if ( 0 == strcasecmp( argv[i] , "--dynamic" ) ) {
      //dynamic scheduling will be used for loops.
      //currently it is default already.
      //so has to use the dynamic schedulin.
      params.use_dynamic_scheduling = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--DENSEACCMAX" ) ) {
      //on CPUs and KNLs if DEFAULT algorithm or KKSPGEMM is chosen,
      //it uses dense accumulators for smaller matrices based on the size of column (k) in B.
      //Max column size is 250,000 for k to use dense accumulators.
      //this parameter overwrites this.
      //with cache mode, or CPUs with smaller thread count, where memory bandwidth is not an issue,
      //this cut-off can be increased to be more than 250,000
      params.MaxColDenseAcc = atoi(getNextArg(i, argc, argv));
    }
    else if ( 0 == strcasecmp( argv[i] , "--verbose" ) ) {
      //print the timing and information about the inner steps.
      params.verbose = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--algorithm" ) ) {
      char* algoStr = getNextArg(i, argc, argv);

      if ( 0 == strcasecmp( algoStr, "DEFAULT" ) ) {
	params.algorithm = KokkosSparse::SPGEMM_KK;
      }
      else if ( 0 == strcasecmp( algoStr, "KKDEFAULT" ) ) {
	params.algorithm = KokkosSparse::SPGEMM_KK;
      }
      else if ( 0 == strcasecmp( algoStr, "KKSPGEMM" ) ) {
	params.algorithm = KokkosSparse::SPGEMM_KK;
      }
      else if ( 0 == strcasecmp( algoStr, "KKMEM" ) ) {
	params.algorithm = KokkosSparse::SPGEMM_KK_MEMORY;
      }
      else if ( 0 == strcasecmp( algoStr, "KKDENSE" ) ) {
        params.algorithm =  KokkosSparse::SPGEMM_KK_DENSE;
      }
      else if ( 0 == strcasecmp( algoStr, "KKLP" ) ) {
	params.algorithm = KokkosSparse::SPGEMM_KK_LP;
      }
      else if ( 0 == strcasecmp( algoStr, "KKDEBUG" ) ) {
	params.algorithm = KokkosSparse::SPGEMM_KK_LP;
      }
      else {
	std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
	print_options();
	return 1;
      }
    }
    else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      print_options();
      return 1;
    }
  }
  return 0;
}

int main (int argc, char ** argv){
  using size_type = default_size_type;
  using lno_t =  default_lno_t;
  using scalar_t = default_scalar;

  KokkosKernels::Experiment::Parameters params;

  if (parse_inputs (params, argc, argv) ){
    return 1;
  }
  if (params.a_mtx_bin_file == NULL){
    std::cerr << "Provide a and b matrix files" << std::endl ;
    print_options();
    return 0;
  }
  if (params.b_mtx_bin_file == NULL){
    std::cout << "B is not provided. Multiplying AxA." << std::endl;
  }

  const int num_threads = std::max(params.use_openmp, params.use_threads);
  const int device_id = params.use_cuda - 1;

  Kokkos::initialize( Kokkos::InitArguments( num_threads, -1, device_id ) );
  Kokkos::print_configuration(std::cout);

#if defined( KOKKOS_ENABLE_OPENMP )

  if (params.use_openmp) {
#ifdef KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE
    KokkosKernels::Experiment::run_spgemm_jacobi
      <size_type, lno_t, scalar_t, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace, Kokkos::HostSpace>(params);
#else 
    KokkosKernels::Experiment::run_spgemm_jacobi
      <size_type, lno_t, scalar_t, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::OpenMP::memory_space>(params);
#endif
  }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
  if (params.use_cuda) {
#ifdef KOKKOSKERNELS_INST_MEMSPACE_CUDAHOSTPINNEDSPACE
    KokkosKernels::Experiment::run_spgemm_jacobi
      <size_type, lno_t, scalar_t, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::CudaHostPinnedSpace>(params);
#else
    KokkosKernels::Experiment::run_spgemm_jacobi
      <size_type, lno_t, scalar_t, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::Cuda::memory_space>(params);
#endif
  }
#endif

#if defined( KOKKOS_ENABLE_THREADS )
  //If only serial is enabled (or no other device was specified), run with serial
  if (params.use_threads)
  {
    KokkosKernels::Experiment::run_spgemm_jacobi
      <size_type, lno_t, scalar_t, Kokkos::Threads, Kokkos::HostSpace, Kokkos::HostSpace>(params);
  }
#endif

#if defined( KOKKOS_ENABLE_SERIAL )
  //If only serial is enabled (or no other device was specified), run with serial
  if (!params.use_openmp && !params.use_cuda && !params.use_threads)
  {
    KokkosKernels::Experiment::run_spgemm_jacobi
      <size_type, lno_t, scalar_t, Kokkos::Serial, Kokkos::HostSpace, Kokkos::HostSpace>(params);
  }
#endif

  Kokkos::finalize(); 

  return 0;
}

