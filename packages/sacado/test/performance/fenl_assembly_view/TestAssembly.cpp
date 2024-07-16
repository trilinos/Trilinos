// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Tests
#include "TestAssembly.hpp"

// Devices
#include "Kokkos_Core.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef KOKKOS_ENABLE_CUDA
#include "cuda_runtime_api.h"
#endif

// For vtune
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
  bool success = true;
  bool verbose = false;
  try {

    const size_t num_sockets = Kokkos::hwloc::get_available_numa_count();
    const size_t num_cores_per_socket =
      Kokkos::hwloc::get_available_cores_per_numa();
    const size_t num_threads_per_core =
      Kokkos::hwloc::get_available_threads_per_core();

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This test performance of MP::Vector FEM assembly.\n");
    int nGrid = 0;
    CLP.setOption("n", &nGrid, "Number of mesh points in each direction.  Set to zero to use a range");
    int nGridBegin = 8;
    CLP.setOption("n-begin", &nGridBegin, "Beginning number of mesh points in each direction.");
    int nGridEnd = 48;
    CLP.setOption("n-end", &nGridEnd, "Ending number of mesh points in each direction.");
    int nGridStep = 8;
    CLP.setOption("n-step", &nGridStep, "Increment in number of mesh points in each direction.");
    int nIter = 10;
    CLP.setOption("ni", &nIter, "Number of assembly iterations");
    bool print = false;
    CLP.setOption("print", "no-print", &print, "Print debugging output");
    bool check = false;
    CLP.setOption("check", "no-check", &check, "Check correctness");
    bool quadratic = false;
    CLP.setOption("quadratic", "linear", &quadratic, "Use quadratic basis functions");
    int num_cores = num_cores_per_socket * num_sockets;
    CLP.setOption("cores", &num_cores,
                  "Number of CPU cores to use (defaults to all)");
    int num_hyper_threads = num_threads_per_core;
    CLP.setOption("hyperthreads", &num_hyper_threads,
                  "Number of hyper threads per core to use (defaults to all)");
#ifdef KOKKOS_ENABLE_THREADS
    bool threads = true;
    CLP.setOption("threads", "no-threads", &threads, "Enable Threads device");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    bool openmp = true;
    CLP.setOption("openmp", "no-openmp", &openmp, "Enable OpenMP device");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = true;
    CLP.setOption("cuda", "no-cuda", &cuda, "Enable Cuda device");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID.");
#endif
    bool vtune = false;
    CLP.setOption("vtune", "no-vtune", &vtune, "connect to vtune");
    CLP.parse( argc, argv );

    if (nGrid > 0) {
      nGridBegin = nGrid;
      nGridEnd = nGrid;
    }

    // Connect to VTune if requested
    if (vtune) {
      std::stringstream cmd;
      pid_t my_os_pid=getpid();
      const std::string vtune_loc =
        "amplxe-cl";
      const std::string output_dir = "./vtune/vtune.0";
      cmd << vtune_loc
          << " -collect hotspots -result-dir " << output_dir
          << " -target-pid " << my_os_pid << " &";
      std::cout << cmd.str() << std::endl;
      system(cmd.str().c_str());
      system("sleep 10");
    }

    Kokkos::initialize(argc,argv);
#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      typedef Kokkos::Threads Device;

      std::cout << std::endl
                << "Threads performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      performance_test_driver<Device>(
        print, nIter, nGridBegin, nGridEnd, nGridStep, quadratic, check);
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      typedef Kokkos::OpenMP Device;

      std::cout << std::endl
                << "OpenMP performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      performance_test_driver<Device>(
        print, nIter, nGridBegin, nGridEnd, nGridStep, quadratic, check);
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      typedef Kokkos::Cuda Device;

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device_id);
      std::cout << std::endl
                << "CUDA performance performance with device " << device_id
                << " ("
                << deviceProp.name << "):"
                << std::endl;

      performance_test_driver<Device>(
        print, nIter, nGridBegin, nGridEnd, nGridStep, quadratic, check);

    }
#endif
    Kokkos::finalize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
