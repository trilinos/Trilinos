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
