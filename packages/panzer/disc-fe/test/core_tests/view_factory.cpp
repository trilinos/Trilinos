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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Sacado.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Phalanx_KokkosUtilities.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer_test {

  TEUCHOS_UNIT_TEST(view_factory, dyn_rank_views)
  {
    using FadType = Sacado::Fad::DFad<double>;
    const int derivative_dim_plus_one = 7;
				  
    // Test a DynRankView from a DynRankView
    {
      Kokkos::DynRankView<FadType,PHX::Device> a("a",10,4,13,derivative_dim_plus_one);
      TEST_EQUALITY(Kokkos::dimension_scalar(a),derivative_dim_plus_one);
      auto b = Kokkos::createDynRankView(a,"b",5,3,8);
      TEST_EQUALITY(Kokkos::dimension_scalar(b),derivative_dim_plus_one);
    }

    // Test a DynRankView from a View
    {
      Kokkos::View<FadType*,PHX::Device> a("a",8,derivative_dim_plus_one);
      TEST_EQUALITY(Kokkos::dimension_scalar(a),derivative_dim_plus_one);
      auto b = Kokkos::createDynRankView(a,"b",5,3,8);
      TEST_EQUALITY(Kokkos::dimension_scalar(b),derivative_dim_plus_one);
    }

  }
}
