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

#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_multimem_spgemm.hpp"

#define SIZE_TYPE int
#define INDEX_TYPE int
#define SCALAR_TYPE double

void print_options(){
  std::cerr << "Options\n" << std::endl;

  std::cerr << "\t[Required] BACKEND: 'threads[numThreads]' | 'openmp [numThreads]' | 'cuda'" << std::endl;

  std::cerr << "\t[Required] INPUT MATRICES: 'amtx [left_hand_side.mtx]'  'bmtx [righ_hand_side.mtx]'" << std::endl;

  std::cerr << "\t[Required] 'algorithm [KKMEM|OUTER|KKSPEED|KKCOLOR|KKMULTICOLOR|KKMULTICOLOR2|MKL|CUSPARSE|CUSP|]'" << std::endl;
  std::cerr << "\t[Optional] OUTPUT MATRICES: 'cmtx [output_matrix.mtx]'" << std::endl;

  std::cerr << "\tThe memory space used for each matrix: 'memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, then compile SPGEMM executable with -DKOKKOSKERNELS_MULTILEVELMEM." << std::endl;
  std::cerr << "\t'CRWC': it will perform hypergraph analysis for memory accesses" << std::endl;
  std::cerr << "\t'CIF path_to_coloring_file': If coloring variants are used, colors will be read from this file." << std::endl;
  std::cerr << "\t'COF path_to_coloring_file': If coloring variants are used, first graph coloring will be performed, then writtent to this file." << std::endl;
  std::cerr << "\tLoop scheduling: 'dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
  std::cerr << "\tVerbose Output: 'verbose'" << std::endl;
}


int parse_inputs (KokkosKernels::Experiment::Parameters &params, int argc, char **argv){
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "threads" ) ) {
      params.use_threads = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "openmp" ) ) {
      params.use_openmp = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda" ) ) {
      params.use_cuda = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "repeat" ) ) {
      params.repeat = atoi( argv[++i] );
    }

    else if ( 0 == strcasecmp( argv[i] , "chunksize" ) ) {
      params.chunk_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "teamsize" ) ) {
      params.team_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "vectorsize" ) ) {
      params.vector_size  = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "memspaces" ) ) {
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
    else if ( 0 == strcasecmp( argv[i] , "CRWC" ) ) {
      params.calculate_read_write_cost = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "CIF" ) ) {
      params.coloring_input_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "COF" ) ) {
      params.coloring_output_file = argv[++i];
    }

    else if ( 0 == strcasecmp( argv[i] , "mcscale" ) ) {
      params.multi_color_scale = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "shmem" ) ) {
      params.shmemsize =  atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "compression2step" ) ) {
      params.compression2step =  true ;
    }
    else if ( 0 == strcasecmp( argv[i] , "mklsort" ) ) {
      params.mkl_sort_option = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "mklkeepout" ) ) {
      params.mkl_keep_output = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "checkoutput" ) ) {
      params.check_output = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "amtx" ) ) {
      params.a_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "cmtx" ) ) {
      params.c_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "bmtx" ) ) {
      params.b_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "dynamic" ) ) {
      params.use_dynamic_scheduling = 1;
    }

    else if ( 0 == strcasecmp( argv[i] , "LLT" ) ) {
      params.left_lower_triangle = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "RLT" ) ) {
      params.right_lower_triangle = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "LS" ) ) {
      params.left_sort = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "RS" ) ) {
      params.right_sort = 1;
    }

    else if ( 0 == strcasecmp( argv[i] , "verbose" ) ) {
      params.verbose = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "algorithm" ) ) {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "MKL" ) ) {
        params.algorithm = 1;
      }
      else if ( 0 == strcasecmp( argv[i] , "CUSPARSE" ) ) {
        params.algorithm = 2;
      }
      else if ( 0 == strcasecmp( argv[i] , "CUSP" ) ) {
        params.algorithm = 3;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKDEBUG" ) ) {
        params.algorithm = 4;
      }
      else if ( 0 == strcasecmp( argv[i] , "MKL2" ) ) {
        params.algorithm = 5;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEM2" ) ) {
        params.algorithm = 6;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEM" ) ) {
        params.algorithm = 7;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKSPEED" ) ) {
        params.algorithm = 8;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKCOLOR" ) ) {
        params.algorithm = 9;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMULTICOLOR" ) ) {
        params.algorithm = 10;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMULTICOLOR2" ) ) {
        params.algorithm = 11;
      }
      else if ( 0 == strcasecmp( argv[i] , "VIENNA" ) ) {
        params.algorithm = 12;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEMSPEED" ) ) {
        params.algorithm = 13;
      }
      else if ( 0 == strcasecmp( argv[i] , "MULTIMEM" ) ) {
        params.algorithm = 14;
      }
      else if ( 0 == strcasecmp( argv[i] , "OUTER" ) ) {
        params.algorithm = 15;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLE" ) ) {
        params.algorithm = 16;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEMEM" ) ) {
        params.algorithm = 17;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEDENSE" ) ) {
        params.algorithm = 18;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEIA" ) ) {
        params.algorithm = 19;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEIAMEM" ) ) {
        params.algorithm = 20;
      }
      else if ( 0 == strcasecmp( argv[i] , "TRIANGLEIADENSE" ) ) {
        params.algorithm = 21;
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

#if defined( KOKKOS_HAVE_OPENMP )

  if (params.use_openmp) {

    Kokkos::OpenMP::initialize( params.use_openmp );
	  Kokkos::OpenMP::print_configuration(std::cout);
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_spgemm
    <SIZE_TYPE, INDEX_TYPE, SCALAR_TYPE, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::HostSpace>(
        params
        );
#else 
    KokkosKernels::Experiment::run_multi_mem_spgemm
    <SIZE_TYPE, INDEX_TYPE, SCALAR_TYPE, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::OpenMP::memory_space>(
        params
        );
#endif
    Kokkos::OpenMP::finalize();
  }

#endif

#if defined( KOKKOS_HAVE_CUDA )
  if (params.use_cuda) {
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( 0 ) );
    Kokkos::Cuda::print_configuration(std::cout);

#ifdef KOKKOSKERNELS_MULTI_MEM

    KokkosKernels::Experiment::run_multi_mem_spgemm
    <SIZE_TYPE, INDEX_TYPE, SCALAR_TYPE, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::CudaHostPinnedSpace>(
        params
        );
#else
    KokkosKernels::Experiment::run_multi_mem_spgemm
    <SIZE_TYPE, INDEX_TYPE, SCALAR_TYPE, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::Cuda::memory_space>(
        params
        );
#endif

    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }

#endif


  return 0;

}







