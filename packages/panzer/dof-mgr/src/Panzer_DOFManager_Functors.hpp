// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_DOFManager_Functors_hpp__
#define __Panzer_DOFManager_Functors_hpp__

#include "Kokkos_Core.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {
namespace dof_functors {

//! Sums all entries of a Rank 2 Kokkos View 
template <typename ReductionDataType, typename view_t>
struct SumRank2 {
    using policy_t = Kokkos::MDRangePolicy<typename view_t::execution_space, Kokkos::Rank<2>>;

    const view_t values;

    void apply(ReductionDataType& sum) const
    {
      const auto& values_ref = values;

      Kokkos::parallel_reduce(
        policy_t({0, 0}, {values.extent(0), values.extent(1)}),
        KOKKOS_LAMBDA(const typename policy_t::index_type indexi, const typename policy_t::index_type indexj, ReductionDataType& local_sum)
        {
          local_sum += values_ref(indexi, indexj);
        },
        Kokkos::Sum<ReductionDataType>(sum)
      );
    }
};

}
}

#endif
