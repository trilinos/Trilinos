// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Sacado.hpp"

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include "Kokkos_ViewFactory.hpp"

TEUCHOS_UNIT_TEST(view_factory, dyn_rank_views)
{
  using FadType = Sacado::Fad::DFad<double>;
  using Kokkos::View;
  using Kokkos::DynRankView;
  using Kokkos::createDynRankView;
  using Kokkos::createDynRankViewWithType;
  using Kokkos::createViewWithType;
  using Kokkos::dimension_scalar;
  using Kokkos::Experimental::view_alloc;
  using Kokkos::Experimental::WithoutInitializing;
  const unsigned derivative_dim_plus_one = 7;

  // Test a DynRankView from a DynRankView
  {
    DynRankView<FadType> a("a",10,4,13,derivative_dim_plus_one);
    TEST_EQUALITY(dimension_scalar(a),derivative_dim_plus_one);
    TEST_EQUALITY(a.rank(),3);

    auto b = createDynRankView(a,"b",5,3,8);
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);
    TEST_EQUALITY(b.rank(),3);

    auto c = createDynRankView(a,view_alloc("c",WithoutInitializing),5,3,8);
    TEST_EQUALITY(dimension_scalar(c),derivative_dim_plus_one);
    TEST_EQUALITY(c.rank(),3);

    using d_type = Kokkos::DynRankView<FadType,Kokkos::LayoutRight>;
    d_type d = createDynRankViewWithType<d_type>(a,"d",5,3,8);
    TEST_EQUALITY(dimension_scalar(d),derivative_dim_plus_one);
    TEST_EQUALITY(d.rank(),3);
  }

  // Test a DynRankView from a View
  {
    View<FadType*> a("a",8,derivative_dim_plus_one);
    TEST_EQUALITY(dimension_scalar(a),derivative_dim_plus_one);

    auto b = createDynRankView(a,"b",5,3,8);
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);
    TEST_EQUALITY(b.rank(),3);

    auto c = createDynRankView(a,view_alloc("c",WithoutInitializing),5,3,8);
    TEST_EQUALITY(dimension_scalar(c),derivative_dim_plus_one);
    TEST_EQUALITY(c.rank(),3);

    using d_type = Kokkos::DynRankView<FadType,Kokkos::LayoutRight>;
    d_type d = createDynRankViewWithType<d_type>(a,"d",5,3,8);
    TEST_EQUALITY(dimension_scalar(d),derivative_dim_plus_one);
    TEST_EQUALITY(d.rank(),3);
  }

  // Test a View from a View
  {
    View<FadType*> a("a",8,derivative_dim_plus_one);
    TEST_EQUALITY(dimension_scalar(a),derivative_dim_plus_one);

    using b_type = Kokkos::View<FadType***>;
    b_type b = createViewWithType<b_type>(a,"b",5,3,8);
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);

    b_type c = createViewWithType<b_type>(a,view_alloc("c",WithoutInitializing),5,3,8);
    TEST_EQUALITY(dimension_scalar(c),derivative_dim_plus_one);

    using d_type = Kokkos::View<FadType***,Kokkos::LayoutRight>;
    d_type d = createViewWithType<d_type>(a,"d",5,3,8);
    TEST_EQUALITY(dimension_scalar(d),derivative_dim_plus_one);
  }

  // Test a View from a DynRankView
  {
    DynRankView<FadType> a("a",10,4,13,derivative_dim_plus_one);
    TEST_EQUALITY(dimension_scalar(a),derivative_dim_plus_one);
    TEST_EQUALITY(a.rank(),3);

    using b_type = Kokkos::View<FadType***>;
    b_type b = createViewWithType<b_type>(a,"b",5,3,8);
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);

    b_type c = createViewWithType<b_type>(a,view_alloc("c",WithoutInitializing),5,3,8);
    TEST_EQUALITY(dimension_scalar(c),derivative_dim_plus_one);

    using d_type = Kokkos::View<FadType***,Kokkos::LayoutRight>;
    d_type d = createViewWithType<d_type>(a,"d",5,3,8);
    TEST_EQUALITY(dimension_scalar(d),derivative_dim_plus_one);
  }

  // Test creation of a Fad DynRankView from a double DynRankView
  {
    DynRankView<double> a("a",10,4,13);
    TEST_EQUALITY(dimension_scalar(a),0);
    TEST_EQUALITY(a.rank(),3);

    using b_type = Kokkos::DynRankView<FadType,Kokkos::LayoutRight>;
    b_type b = createDynRankViewWithType<b_type>(a,"b",5,3,8);
    TEST_EQUALITY(dimension_scalar(b),1);
    TEST_EQUALITY(b.rank(),3);
  }

  // Test a double DynRankView from a double DynRankView
  {
    DynRankView<double> a("a",10,4,13);
    TEST_EQUALITY(dimension_scalar(a),0);
    TEST_EQUALITY(a.rank(),3);

    auto b = createDynRankView(a,"b",5,3,8);
    TEST_EQUALITY(dimension_scalar(b),0);
    TEST_EQUALITY(b.rank(),3);
  }

  // Test double rank 0
  {
    DynRankView<double> a("a",10,4,13);
    TEST_EQUALITY(dimension_scalar(a),0);
    TEST_EQUALITY(a.rank(),3);

    auto b = createDynRankView(a,"b");
    TEST_EQUALITY(dimension_scalar(b),0);
    TEST_EQUALITY(b.rank(),0);
  }

  // Test Fad rank 0
  {
    DynRankView<FadType> a("a",10,4,13,derivative_dim_plus_one);
    TEST_EQUALITY(dimension_scalar(a),derivative_dim_plus_one);
    TEST_EQUALITY(a.rank(),3);

    auto b = createDynRankView(a,"b");
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);
    TEST_EQUALITY(b.rank(),0);
  }

}

#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize();

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  Kokkos::finalize();

  return res;
}
