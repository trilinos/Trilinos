// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Sacado.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer_test {

  TEUCHOS_UNIT_TEST(view_factory, dyn_rank_views)
  {
    using FadType = Sacado::Fad::DFad<double>;
    const int derivative_dim_plus_one = 7;
				  
    // Test a DynRankView from a DynRankView
    {
      Kokkos::DynRankView<FadType,PHX::Device> a("a",10,4,13,derivative_dim_plus_one);
      TEST_EQUALITY(static_cast<int>(Kokkos::dimension_scalar(a)),derivative_dim_plus_one);
      auto b = Kokkos::createDynRankView(a,"b",5,3,8);
      TEST_EQUALITY(static_cast<int>(Kokkos::dimension_scalar(b)),derivative_dim_plus_one);
    }

    // Test a DynRankView from a View
    {
      PHX::View<FadType*> a("a",8,derivative_dim_plus_one);
      TEST_EQUALITY(static_cast<int>(Kokkos::dimension_scalar(a)),derivative_dim_plus_one);
      auto b = Kokkos::createDynRankView(a,"b",5,3,8);
      TEST_EQUALITY(static_cast<int>(Kokkos::dimension_scalar(b)),derivative_dim_plus_one);
    }

  }
}
