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
#include "KokkosKernels_config.h"
#if defined(KOKKOSKERNELS_INST_DOUBLE) &&  \
    defined(KOKKOSKERNELS_INST_OFFSET_INT) && \
    defined(KOKKOSKERNELS_INST_ORDINAL_INT)
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_TestParameters.hpp"
#include "KokkosSparse_run_spadd.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#define SIZE_TYPE int
#define INDEX_TYPE int
#define SCALAR_TYPE double
//double

void print_options(){
  std::cerr << "Options\n" << std::endl;

  std::cerr << "\t[Required] BACKEND: '--threads[numThreads]' | '--openmp [numThreads]' | '--cuda [cudaDeviceIndex]'" << std::endl;

  std::cerr << "\t[Required] --amtx <path> :: 1st input matrix" << std::endl;
  std::cerr << "\t[Required] --bmtx <path> :: 2nd input matrix" << std::endl;
  std::cerr << "\t[Optional] --cmtx <path> :: output matrix for C = A+B"  << std::endl;
  std::cerr << "\t[Optional] Verbose Output: '--verbose'" << std::endl;
}


int parse_inputs (KokkosKernels::Experiment::Parameters &params, int argc, char **argv){
  params.assume_sorted = false;
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "--threads" ) ) {
      params.use_threads = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--openmp" ) ) {
      params.use_openmp = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--cuda" ) ) {
      params.use_cuda = atoi( argv[++i] ) + 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--sorted" ) ) {
      params.assume_sorted = true;
    }
    else if ( 0 == strcasecmp( argv[i] , "--unsorted" ) ) {
      params.assume_sorted = false;
    }
    else if ( 0 == strcasecmp( argv[i] , "--amtx" ) ) {
      //A at C=AxB
      params.a_mtx_bin_file = argv[++i];
    }

    else if ( 0 == strcasecmp( argv[i] , "--bmtx" ) ) {
      //B at C=AxB.
      //if not provided, C = AxA will be performed.
      params.b_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "--cmtx" ) ) {
      //if provided, C will be written to given file.
      //has to have ".bin", or ".crs" extension.
      params.c_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "--verbose" ) ) {
      params.verbose = true;
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
  if (!params.a_mtx_bin_file || !params.b_mtx_bin_file) {
    std::cerr << "Provide a and b matrix files" << std::endl;
    print_options();
    return 0;
  }

  const int num_threads = params.use_openmp; // Assumption is that use_openmp variable is provided as number of threads
  const int device_id = params.use_cuda - 1;

  Kokkos::initialize( Kokkos::InitArguments( num_threads, -1, device_id ) );
  Kokkos::print_configuration(std::cout);

  bool useOMP = params.use_openmp != 0;
  bool useCUDA = params.use_cuda != 0;
  bool useSerial = !useOMP && !useCUDA;

  if(useOMP)
  {
#if defined( KOKKOS_ENABLE_OPENMP )
    using crsMat_t = KokkosSparse::CrsMatrix<double, int, Kokkos::OpenMP, void, size_t>;
    KokkosKernels::Experiment::run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }
  if(useCUDA)
  {
#if defined( KOKKOS_ENABLE_CUDA )
    using crsMat_t = KokkosSparse::CrsMatrix<double, int, Kokkos::Cuda, void, size_t>;
    KokkosKernels::Experiment::run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }
  if(useSerial)
  {
#if defined( KOKKOS_ENABLE_SERIAL )
    using crsMat_t = KokkosSparse::CrsMatrix<double, int, Kokkos::Serial, void, size_t>;
    KokkosKernels::Experiment::run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }
  Kokkos::finalize(); 
  return 0;
}


#else
int main() {
#if !defined(KOKKOSKERNELS_INST_DOUBLE)
std::cout  << " not defined KOKKOSKERNELS_INST_DOUBLE"  << std::endl;
#endif

#if !defined(KOKKOSKERNELS_INST_OFFSET_INT)
std::cout  << " not defined KOKKOSKERNELS_INST_OFFSET_INT"  << std::endl;

#endif

#if !defined(KOKKOSKERNELS_INST_ORDINAL_INT)
std::cout  << " not defined KOKKOSKERNELS_INST_ORDINAL_INT"  << std::endl;

#endif
}
#endif

