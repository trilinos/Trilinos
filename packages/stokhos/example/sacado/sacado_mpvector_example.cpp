// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "sacado_mpvector_example.hpp"

#include "KokkosCore_config.h"
#include "Kokkos_Threads.hpp"
#include "Kokkos_hwloc.hpp"

int main(int argc, char **argv)
{
  try {

    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example explores operator overloading on CUDA.\n");

    const int num_storage_method = 12;
    const Storage_Method storage_method_values[] = {
      STATIC,
      STATIC_FIXED,
      LOCAL,
      DYNAMIC,
      DYNAMIC_STRIDED,
      DYNAMIC_THREADED,
      VIEW_STATIC,
      VIEW_STATIC_FIXED,
      VIEW_LOCAL,
      VIEW_DYNAMIC,
      VIEW_DYNAMIC_STRIDED,
      VIEW_DYNAMIC_THREADED };
    const char *storage_method_names[] = {
      "static",
      "static-fixed",
      "local",
      "dynamic",
      "dynamic-strided",
      "dynamic-threaded",
      "view-static",
      "view-static-fixed",
      "view-local",
      "view-dynamic",
      "view-dynamic-strided",
      "view-dynamic-threaded",};
    Storage_Method storage_method = STATIC_FIXED;
    CLP.setOption("storage_method", &storage_method,
                  num_storage_method, storage_method_values,
                  storage_method_names, "Storage method");

    int num_elements = 100;
    CLP.setOption("num_elements", &num_elements, "Number of elements");

    int num_samples = 32;
    CLP.setOption("num_samples", &num_samples, "Number of samples");

    bool reset = true;
    CLP.setOption("reset", "no-reset", &reset, "Reset initial vectors");

    bool print = false;
    CLP.setOption("print", "quiet", &print, "Print values");

#ifdef KOKKOS_HAVE_PTHREAD
    bool test_threads = true;
    CLP.setOption("threads", "no-threads", &test_threads, "Test Threads");

    int threads_team_size = -1;
    CLP.setOption("threads_team_size", &threads_team_size,
                  "Thread team size (use -1 for core capacity)");

    int threads_league_size = -1;
    CLP.setOption("threads_league_size", &threads_league_size,
                  "Thread league size (use -1 for # of cores)");
#endif

#ifdef KOKKOS_HAVE_CUDA
    bool test_cuda = true;
    CLP.setOption("cuda", "no-cuda", &test_cuda, "Test CUDA");

    int cuda_team_size = -1;
    CLP.setOption("cuda_team_size", &cuda_team_size,
                  "Cuda team size (use -1 for all samples)");

    int cuda_league_size = -1;
#else
    bool test_cuda = false;
#endif

    CLP.parse( argc, argv );

    std::cout << "Summary of command line options:" << std::endl
              << "\tstorage_method        = "
              << storage_method_names[storage_method] << std::endl
              << "\tnum_elements          = " << num_elements << std::endl
              << "\tnum_samples           = " << num_samples << std::endl
              << "\treset                 = " << reset << std::endl
              << "\tprint                 = " << print << std::endl
              << "\tthreads               = " << test_threads << std::endl
              << "\tthreads_team_size     = " << threads_team_size << std::endl
              << "\tthreads_league_size   = " << threads_league_size << std::endl
              << "\tcuda                  = " << test_cuda << std::endl
              << "\tcuda_team_size        = " << cuda_team_size << std::endl
              << std::endl;

#ifdef KOKKOS_HAVE_PTHREAD
    // Always initialize threads for Cuda::host_mirror
    if (threads_team_size == -1) {
      // So that default 'num_samples / threads_team_size' is 8 as required by the test
      threads_team_size = 4 ;
    }
    if (threads_league_size == -1) {
      // Don't oversubscribe the cores.
      threads_league_size =
        ( Kokkos::hwloc::get_available_numa_count() *
          Kokkos::hwloc::get_available_cores_per_numa() *
          Kokkos::hwloc::get_available_threads_per_core() ) / threads_team_size ;
    }

    Kokkos::Threads::initialize( threads_league_size * threads_team_size );

    if (test_threads) {
       bool status = MPVectorExample<MaxSize,Scalar,Kokkos::Threads>::run(
         storage_method, num_elements, num_samples,
         threads_team_size, threads_league_size,
         reset, print);

       if (status)
         std::cout << "Threads Test Passed!" << std::endl;
       else
         std::cout << "Threads Test Failed!" << std::endl;
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    if (test_cuda) {
      bool status = MPVectorExample<MaxSize,Scalar,Kokkos::Cuda>::run(
        storage_method, num_elements, num_samples,
        cuda_team_size, cuda_league_size,
        reset, print);

      if (status)
        std::cout << "CUDA Test Passed!" << std::endl;
      else
        std::cout << "CUDA Test Failed!" << std::endl;
    }
#endif

#ifdef KOKKOS_HAVE_PTHREAD
    // Finalize threads
    Kokkos::Threads::finalize();
#endif

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (std::string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std::endl;
  }
}
