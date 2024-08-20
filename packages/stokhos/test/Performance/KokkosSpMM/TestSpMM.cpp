// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Devices
#include "Kokkos_Core.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef KOKKOS_ENABLE_CUDA
#include "cuda_runtime_api.h"
#endif

template <typename Scalar, typename Ordinal, typename Device>
void performance_test_driver( const Ordinal nGrid,
                              const Ordinal nIter,
                              const Ordinal ensemble_min,
                              const Ordinal ensemble_max,
                              const Ordinal ensemble_step );

int main(int argc, char *argv[])
{
  bool success = true;
  bool verbose = false;
  try {

#ifdef KOKKOS_ENABLE_THREADS
    const size_t num_sockets = Kokkos::hwloc::get_available_numa_count();
    const size_t num_cores_per_socket =
      Kokkos::hwloc::get_available_cores_per_numa();
    const size_t num_threads_per_core =
      Kokkos::hwloc::get_available_threads_per_core();
#endif

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This test performance of MP::Vector multiply routines.\n");
    int nGrid = 32;
    CLP.setOption("n", &nGrid, "Number of mesh points in the each direction");
    int nIter = 10;
    CLP.setOption("ni", &nIter, "Number of multiply iterations");
    int ensemble_min = 4;
    CLP.setOption("emin", &ensemble_min, "Staring ensemble size");
    int ensemble_max = 24;
    CLP.setOption("emax", &ensemble_max, "Stoping ensemble size");
    int ensemble_step = 4;
    CLP.setOption("estep", &ensemble_step, "Ensemble increment");
#ifdef KOKKOS_ENABLE_THREADS
    bool threads = true;
    CLP.setOption("threads", "no-threads", &threads, "Enable Threads device");
    int num_cores = num_cores_per_socket * num_sockets;
    CLP.setOption("cores", &num_cores,
                  "Number of CPU cores to use (defaults to all)");
    int num_hyper_threads = num_threads_per_core;
    CLP.setOption("hyperthreads", &num_hyper_threads,
                  "Number of hyper threads per core to use (defaults to all)");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = true;
    CLP.setOption("cuda", "no-cuda", &cuda, "Enable Cuda device");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
#endif
    CLP.parse( argc, argv );

    typedef int Ordinal;
    typedef double Scalar;

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      typedef Kokkos::Threads Device;

      Kokkos::InitializationSettings init_args;
      init_args.set_num_threads(num_cores*num_hyper_threads);
      Kokkos::initialize( init_args );

      std::cout << std::endl
                << "Threads performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, ensemble_min, ensemble_max, ensemble_step);

      Kokkos::finalize();
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      typedef Kokkos::Cuda Device;

      Kokkos::InitializationSettings init_args;
      init_args.set_device_id(device_id);
      Kokkos::initialize( init_args );

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device_id);
      std::cout << std::endl
                << "CUDA performance for device " << device_id << " ("
                << deviceProp.name << "):"
                << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, ensemble_min, ensemble_max, ensemble_step);

      Kokkos::finalize();
    }
#endif

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
