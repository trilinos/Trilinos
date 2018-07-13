/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <iostream>

#include "KokkosGraph_multimem_triangle.hpp"
#include "KokkosKernels_IOUtils.hpp"







void print_options(){
  std::cerr << "Options\n" << std::endl;
  std::cerr << "Choose BackEnd                     : --openmp [numthreads] | --cuda" << std::endl;
  std::cerr << "Input Matrix                       : --amtx [path_to_input_matrix]" << std::endl;
  std::cerr << "\tInput Matrix format can be multiple formats. If it ends with:" << std::endl;
  std::cerr << "\t\t.mtx: it will read matrix market format." << std::endl;
  std::cerr << "\t\t.bin: it will read binary crs matrix format." << std::endl;
  std::cerr << "\t\t.crs: it will read text crs matrix format." << std::endl;
  std::cerr << "--algorithm                          :" << std::endl;
  std::cerr << "\tTRIANGLEAI: for Adj x Incidence" << std::endl;
  std::cerr << "\tTRIANGLEIA: for Incidence x Adj -- implementing set intersection (2D) -- 3rd fastest"  << std::endl;
  std::cerr << "\tTRIANGLEIAUNION: for Incidence x Adj -- implementing set union " << std::endl;
  std::cerr << "\tTRIANGLELL: Lower x Lower -- usually fastest " << std::endl;
  std::cerr << "\tTRIANGLELU: Lower x Upper -- usually 2nd fastest " << std::endl;
  std::cerr << "--FLOP                               : Calculate and print the number of operations. This will be calculated on the first run." << std::endl;
  std::cerr << "--COMPRESSION [0|1]                   : Enable disable compression. Default:1." << std::endl;
  std::cerr << "--RS [0|1|2]                         : Whether to sort lower triangular matrix. 0 - no sort, 1 - sort, 2 - algorithm decides based on max row size (default)" << std::endl;
  std::cerr << "--accumulator [default|dense|sparse] : what type of accumulator to use." << std::endl;
  std::cerr << "--RLT                                : If given, lower triangle will be used for AdjxIncidence or Incidence x Adj algorithms." << std::endl;
  std::cerr << "--dynamic                            : If set, dynamic schedule will be used. Currently default is dynamic scheduling as well." << std::endl;
  std::cerr << "--verbose                            : If set, the inner timer stats will be printed." << std::endl;
  std::cerr << "--repeat [repeatnum]                 : how many repeats will be run." << std::endl;
  std::cerr << "--chunksize [chunksize]              : how many vertices are executed with in a loop index. Default is 16." << std::endl;
  std::cerr << "--sort_option [0|1|2]                : How lower triangle will be sorted. 0: for largest to bottom, 1 for largest to top, 2 for interleaved." << std::endl;
  std::cerr << "--cache_flush [0|1|2]                : Flush between repetitions. 0 - no flush, 1 - soft flush, 2 - hard flush with random numbers." << std::endl;

  std::cerr << "\nSuggested use of LL: executable --amtx path_to_file.bin --algorithm TRIANGLELL --repeat 6 --verbose --chunksize [4|16]" << std::endl;
  std::cerr << "Suggested use of LU: executable --amtx path_to_file.bin --algorithm TRIANGLELU --repeat 6 --verbose --chunksize [4|16]" << std::endl;
  std::cerr << "Suggested use of AI: executable --amtx path_to_file.bin --algorithm TRIANGLEIA --repeat 6 --verbose --chunksize [4|16] rlt" << std::endl;

}


int parse_inputs (KokkosKernels::Experiment::Parameters &params, int argc, char **argv){
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "--threads" ) ) {
      params.use_threads = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--openmp" ) ) {
      params.use_openmp = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--cuda" ) ) {
      params.use_cuda = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--repeat" ) ) {
      params.repeat = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--triangle_operation" ) ) {
      params.triangle_options = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--chunksize" ) ) {
      params.chunk_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--teamsize" ) ) {
      params.team_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--vectorsize" ) ) {
      params.vector_size  = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--compression" ) ) {
      params.apply_compression = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--sort_option" ) ) {
      params.sort_option = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--memspaces" ) ) {
      int memspaces = atoi( argv[++i] ) ;
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
    else if ( 0 == strcasecmp( argv[i] , "--flop" ) ) {
      params.calculate_read_write_cost = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--CIF" ) ) {
      params.coloring_input_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "--COF" ) ) {
      params.coloring_output_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "--mhscale" ) ) {
      params.minhashscale = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--mcscale" ) ) {
      params.multi_color_scale = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--shmem" ) ) {
      params.shmemsize =  atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--compression2step" ) ) {
      params.compression2step =  true ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--mklsort" ) ) {
      params.mkl_sort_option = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--mklkeepout" ) ) {
      params.mkl_keep_output = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--checkoutput" ) ) {
      params.check_output = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--amtx" ) ) {
      params.a_mtx_bin_file = argv[++i];
    }
    /*
    else if ( 0 == strcasecmp( argv[i] , "cmtx" ) ) {
      params.c_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "bmtx" ) ) {
      params.b_mtx_bin_file = argv[++i];
    }
    */
    else if ( 0 == strcasecmp( argv[i] , "--dynamic" ) ) {
      params.use_dynamic_scheduling = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--cache_flush" ) ) {
      params.cache_flush = atoi(argv[++i]);
    }


    else if ( 0 == strcasecmp( argv[i] , "--RLT" ) ) {
      params.right_lower_triangle = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--RS" ) ) {
      params.right_sort = atoi(argv[++i]);
    }

    else if ( 0 == strcasecmp( argv[i] , "--verbose" ) ) {
      params.verbose = 1;
    }

    else if ( 0 == strcasecmp( argv[i] , "--accumulator" ) ) {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "default" ) ) {
        params.accumulator = 0;
      }
      else if ( 0 == strcasecmp( argv[i] , "dense" ) ) {
        params.accumulator = 1;
      }
      else if ( 0 == strcasecmp( argv[i] , "sparse" ) ) {
        params.accumulator = 2;
      }
      else {
        std::cerr << "1-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        print_options();
        return 1;
      }
    }
    else if ( 0 == strcasecmp( argv[i] , "--algorithm" ) ) {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "TRIANGLEAI" ) ) {
        params.algorithm = 16;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEIA" ) ) {
        params.algorithm = 17;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEIAUNION" ) ) {
        params.algorithm = 18;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLELL" ) ) {
        params.algorithm = 19;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLELU" ) ) {
        params.algorithm = 20;
      }
      else {
        std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        print_options();
        return 1;
      }
    }
    else {
      std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      print_options();
      return 1;
    }
  }
  return 0;
}

int main (int argc, char ** argv){

  typedef unsigned size_type;
  typedef int idx;


  KokkosKernels::Experiment::Parameters params;



  if (parse_inputs (params, argc, argv) ){
    return 1;
  }
  if (params.a_mtx_bin_file == NULL){
    std::cerr << "Provide a matrix file" << std::endl ;
    print_options();
    return 0;
  }

  std::cout << "Sizeof(idx):" << sizeof(idx) << " sizeof(size_type):" << sizeof(size_type) << std::endl;

  const int num_threads = params.use_openmp; // Assumption is that use_openmp variable is provided as number of threads
  const int device_id = 0;
  Kokkos::initialize( Kokkos::InitArguments( num_threads, -1, device_id ) );

#if !defined (KOKKOS_ENABLE_CUDA)
#if defined( KOKKOS_ENABLE_OPENMP )

  if (params.use_openmp) {
	  Kokkos::OpenMP::print_configuration(std::cout);
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_triangle
    <size_type, idx, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::HostSpace>(
        params
        );
#else
    KokkosKernels::Experiment::run_multi_mem_triangle
    <size_type, idx, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::OpenMP::memory_space>(
        params
        );
#endif
  }

#endif
#endif


#if defined( KOKKOS_ENABLE_CUDA1 )
  if (params.use_cuda) {
    Kokkos::Cuda::print_configuration(std::cout);
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_triangle
    <size_type, idx, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::CudaHostPinnedSpace>(
        params
        );
#else 
    KokkosKernels::Experiment::run_multi_mem_triangle
    <size_type, idx, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::Cuda::memory_space>(
        params
        );
#endif 
  }

#endif

  Kokkos::finalize();

  return 0;
}







