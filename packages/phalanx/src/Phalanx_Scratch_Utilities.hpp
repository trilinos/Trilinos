// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_SCRATCH_VIEW_HPP
#define PHALANX_SCRATCH_VIEW_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"

namespace PHX {

  // template<typename ScalarT>
  // using ScratchView = Kokkos::View<ScalarT* ,typename PHX::DevLayout<ScalarT>::type,typename PHX::exec_space::scratch_memory_space,Kokkos::MemoryUnmanaged>;

  // Returns the Fad derivtive size of an already allocated View. This
  // can be called with non-FAD scalar types.
  template<typename ScalarT,typename... Props>
  KOKKOS_INLINE_FUNCTION
  std::size_t getFadSize(const Kokkos::View<ScalarT*,Props...>& view)
  { return Kokkos::dimension_scalar(view); }

  /*
  // Binds the shared memory to a view. 
  template<typename ScalarT,typename... Indices>
  KOKKOS_INLINE_FUNCTION
  void bindSharedMemory(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team,
			const PHX::ScratchView<ScalarT*>& view,
			const Indices & ... indices,
			const std::size_t fad_size = 0)
  {
    if (Sacado::IsADType<ScalarT>::value)
      view = PHX::ScratchView(team.team_shmem(),indices...,fad_size);
    else
      view = PHX::ScratchView(team.team_shmem(),indices...);
  }

  template<typename ScalarT,typename... Indices>
  std::size_t bytesToAllocate(const Indices & ... indices,
			      const std::size_t fad_size = 0)
  {
    if (Sacado::IsADType<ScalarT>::value)
      view = PHX::ScratchView::shmem_size(indices...,fad_size);
    else
      view = PHX::ScratchView::shmem_size(indices...);
  }

  template<typename ScalarT, typename ... TeamPolicyProperties, typename ExecSpace>
  Kokkos::TeamPolicy<ExecSpace, TeamPolicyProperties...> createTeamPolicy(ExecSpace exec_space, const int& league_size, scratch_bytes_per_team)
  {
#ifdef KOKKOS_HAVE_CUDA
    constexpr int vector_size = 32;
#elif KOKKOS_HAVE_HIP
    constexpr int vector_size = 64;
#else
    constexpr int vector_size = 1;
#endif

    return Kokkos::TeamPolicy<TeamPolicyProperties...>(exec_space,league_size,Kokkos::AUTO(),vector_size).set_scratch_size(0,Kokkos::PerTeam(scratch_bytes));
  }
  */
}

#endif
