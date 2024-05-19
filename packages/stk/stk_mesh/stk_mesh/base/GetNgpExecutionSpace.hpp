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
namespace impl {
struct CudaStreamDeleter {
  void operator()(cudaStream_t* stream) const {
    cudaStreamDestroy(*stream);
    delete stream;
  }
};

struct ExecSpaceAndCudaStreamDeleter {
  ExecSpaceAndCudaStreamDeleter() {}
  ExecSpaceAndCudaStreamDeleter(cudaStream_t* streamPtr_)
    : streamPtr(streamPtr_)
  {}
  ExecSpaceAndCudaStreamDeleter(const ExecSpaceAndCudaStreamDeleter& deleter)
  {
    streamPtr = deleter.streamPtr;
  }

  void operator()(stk::ngp::ExecSpace* e) const {
    delete e;
    cudaStreamDestroy(*streamPtr);
    delete streamPtr;
  }

  cudaStream_t* streamPtr;
};
}

template<>
class ExecSpaceWrapper<Kokkos::Cuda> {
 public:
  using ExecSpaceType = Kokkos::Cuda;

  ExecSpaceWrapper()
   : stream(new cudaStream_t) {
    cudaStreamCreate(stream);
    space = std::shared_ptr<ExecSpaceType>(new ExecSpaceType(*stream), impl::ExecSpaceAndCudaStreamDeleter(stream));
  }

  operator const ExecSpaceType&() { return *space; }
  const char* name() const { return space->name(); }
  void fence() const { space->fence(); }
  const ExecSpaceType& get_execution_space() const { return *space; }

 private:
  std::shared_ptr<ExecSpaceType> space;
  cudaStream_t* stream;
};
#endif

inline ExecSpaceWrapper<> get_execution_space_with_stream()
{
  bool launchBlockingIsOn = get_env_var_as_bool("CUDA_LAUNCH_BLOCKING", false);
  static bool printedOnce = false;

  if(launchBlockingIsOn && !printedOnce) {
    sierra::Env::outputP0() << "stk::mesh::get_execution_space_with_stream: CUDA_LAUNCH_BLOCKING is ON. Asynchronous operations will block." << std::endl;
    printedOnce = true;
  }

  auto execSpace = ExecSpaceWrapper<>();
  return execSpace;
}

}
}

#endif
