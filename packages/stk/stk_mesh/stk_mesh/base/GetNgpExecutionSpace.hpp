// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef NGP_EXECUTION_SPACE_HPP
#define NGP_EXECUTION_SPACE_HPP

#include <stk_util/util/GetEnv.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>

namespace stk {
namespace mesh {

template <typename ExecSpaceType = stk::ngp::ExecSpace>
class ExecSpaceWrapper {
 public:
  operator const ExecSpaceType&() { return space; }
  const char* name() const { return space.name(); }
  void fence() const { space.fence(); }
  const ExecSpaceType& get_execution_space() const { return space; }

 private:
  ExecSpaceType space;
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
class ExecSpaceWrapper<Kokkos::Cuda> {
 public:
  using ExecSpaceType = Kokkos::Cuda;

  ExecSpaceWrapper() {
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    space = ExecSpaceType(stream, Kokkos::Impl::ManageStream::yes);
  }

  operator const ExecSpaceType&() { return space; }
  const char* name() const { return space.name(); }
  void fence() const { space.fence(); }
  const ExecSpaceType& get_execution_space() const { return space; }

 private:
  ExecSpaceType space;
};

#elif defined(KOKKOS_ENABLE_HIP)
template<>
class ExecSpaceWrapper<Kokkos::HIP> {
 public:
  using ExecSpaceType = Kokkos::HIP;

  ExecSpaceWrapper() {
    hipStream_t stream;
    [[maybe_unused]] auto success = hipStreamCreate(&stream);
    space = ExecSpaceType(stream, Kokkos::Impl::ManageStream::yes);
  }

  operator const ExecSpaceType&() { return space; }
  const char* name() const { return space.name(); }
  void fence() const { space.fence(); }
  const ExecSpaceType& get_execution_space() const { return space; }

 private:
  ExecSpaceType space;
};
#endif

inline ExecSpaceWrapper<> get_execution_space_with_stream()
{
  std::string launchBlockingEnvVar;
  bool launchBlockingIsOn = false;
  static bool printedOnce = false;

#ifdef KOKKOS_ENABLE_CUDA
  launchBlockingEnvVar = "CUDA_LAUNCH_BLOCKING";
  launchBlockingIsOn = get_env_var_as_bool(launchBlockingEnvVar, false);
#elif defined(KOKKOS_ENABLE_HIP)
  launchBlockingEnvVar = "HIP_LAUNCH_BLOCKING";
  launchBlockingIsOn = get_env_var_as_bool(launchBlockingEnvVar, false);
#endif

  if (launchBlockingIsOn && !printedOnce) {
    sierra::Env::outputP0()
      << "stk::mesh::get_execution_space_with_stream: " << launchBlockingEnvVar
      << " is ON. Asynchronous operations will block." << std::endl;
    printedOnce = true;
  }

  return ExecSpaceWrapper<>();
}

}
}

#endif
