// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CHECKLAUNCHBLOCKING_HPP
#define TPETRA_DETAILS_CHECKLAUNCHBLOCKING_HPP

#include "TpetraCore_config.h"
#include <cstdlib>
#include <stdexcept>
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {
#ifdef HAVE_TPETRACORE_CUDA
  //Verify that for pre-Pascal CUDA architectures, $CUDA_LAUNCH_BLOCKING == 1
  inline void checkOldCudaLaunchBlocking()
  {
    if(!Kokkos::is_initialized())
      throw std::logic_error("Kokkos must be initialized in order to check CUDA_LAUNCH_BLOCKING setting.");
    size_t arch = Kokkos::Cuda::device_arch();
    if(arch < 600)
    {
      //Compiling for Kepler/Maxwell: require launch blocking.
      const char* launchBlockingEnv = std::getenv("CUDA_LAUNCH_BLOCKING");
      if(!launchBlockingEnv || strcmp(launchBlockingEnv, "1"))
      {
        throw std::runtime_error(
         "Tpetra::initialize(): Kokkos was compiled for an older CUDA architecture than Pascal, but\n"
         "the environment variable CUDA_LAUNCH_BLOCKING is either unset or is not \"1\".\n"
         "It must be set to \"1\" at runtime.\n");
      }
    }
  }
#else
  inline void checkOldCudaLaunchBlocking() {}
#endif
}}

#endif
