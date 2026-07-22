// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#if defined(KOKKOS_ENABLE_HIP)
    fad_vector_size_ = 64;
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
