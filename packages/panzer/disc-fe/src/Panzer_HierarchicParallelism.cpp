// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_HierarchicParallelism.hpp"

namespace panzer {

  HP::HP() :
    use_auto_team_size_(true),
    team_size_(-1),
    vector_size_(1),
    fad_vector_size_(1),
    use_shared_memory_(true),
    fad_use_shared_memory_(false)
  {
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
#if defined(KOKKOS_ENABLE_CUDA)
    fad_vector_size_ = 32;
#endif
#endif
  }

  HP& HP::inst()
  {
    static HP hp;
    return hp;
  }

  namespace {
    int roundDownToPowerOfTwo(int in) {
      int out=1;
      while (in > 1) {
        out *= 2;
        in /= 2;
      }
      return out;
    }
  }
  void HP::overrideSizes(const int& in_team_size,
			 const int& in_vector_size,
			 const int& in_fad_vector_size,
                         const bool force_override)
  {
    use_auto_team_size_ = false;
    if ( force_override ) {
      team_size_=in_team_size;
      vector_size_=in_vector_size;
      fad_vector_size_=in_fad_vector_size;
      return;
    }

    Kokkos::TeamPolicy<PHX::Device> policy(1, Kokkos::AUTO);
    auto blank_functor = KOKKOS_LAMBDA ( const Kokkos::TeamPolicy<PHX::exec_space>::member_type) {};

    int team_size_max = std::min(in_team_size, policy.team_size_max(blank_functor, Kokkos::ParallelForTag()));
    team_size_=roundDownToPowerOfTwo(team_size_max);

    int vec_size_max = policy.vector_length_max();
    vector_size_ = roundDownToPowerOfTwo(std::min(vec_size_max, in_vector_size));
    fad_vector_size_ = roundDownToPowerOfTwo(std::min(vec_size_max, in_fad_vector_size));
  }

  void HP::setUseSharedMemory(const bool& in_use_shared_memory,
			      const bool& in_fad_use_shared_memory)
  {
    use_shared_memory_ = in_use_shared_memory;
    fad_use_shared_memory_ = in_fad_use_shared_memory;
  }

}
