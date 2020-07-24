// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
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
