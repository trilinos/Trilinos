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

  // Test unmanaged view of double
  {
    Kokkos::View<double*> a("a",5*3);
    using b_type = Kokkos::View<double**,Kokkos::MemoryUnmanaged>;
    b_type b = createViewWithType<b_type>(a,a.data(),5,3);
    TEST_EQUALITY(b.dimension_0(),5);
    TEST_EQUALITY(b.dimension_1(),3);
    TEST_EQUALITY(dimension_scalar(b),0);
  }

  // Test unmanaged view of Fad
  {
    Kokkos::View<FadType*> a("a",5*3,derivative_dim_plus_one);
    using b_type = Kokkos::View<FadType**,Kokkos::MemoryUnmanaged>;
    b_type b = createViewWithType<b_type>(a,a.data(),5,3);
    TEST_EQUALITY(b.dimension_0(),5);
    TEST_EQUALITY(b.dimension_1(),3);
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);
  }

  // Test LayoutStride view of double
  {
    Kokkos::DynRankView<double> a("a",10,13);
    auto b = Kokkos::subview(a, std::make_pair(4,8), std::make_pair(5,11));
    auto c = createDynRankView(b,"c",5,3);
    using b_type = decltype(b);
    using c_type = decltype(c);
    using b_layout = typename b_type::array_layout;
    using c_layout = typename c_type::array_layout;
    using default_layout = typename b_type::device_type::execution_space::array_layout;
    const bool is_b_layout_stride =
      std::is_same<b_layout,Kokkos::LayoutStride>::value;
    const bool is_c_default_layout =
      std::is_same<c_layout,default_layout>::value;
    TEST_EQUALITY(is_b_layout_stride,true);
    TEST_EQUALITY(is_c_default_layout,true);
    TEST_EQUALITY(c.rank(),2);
    TEST_EQUALITY(c.dimension_0(),5);
    TEST_EQUALITY(c.dimension_1(),3);
    TEST_EQUALITY(dimension_scalar(b),0);
  }

  // Test LayoutStride view of Fad
  {
    Kokkos::DynRankView<FadType> a("a",10,13,derivative_dim_plus_one);
    auto b = Kokkos::subview(a, std::make_pair(4,8), std::make_pair(5,11));
    auto c = createDynRankView(b,"c",5,3);
    using b_type = decltype(b);
    using c_type = decltype(c);
    using b_layout = typename b_type::array_layout;
    using c_layout = typename c_type::array_layout;
    using default_layout = typename b_type::device_type::execution_space::array_layout;
    const bool is_b_layout_stride =
      std::is_same<b_layout,Kokkos::LayoutStride>::value;
    const bool is_c_default_layout =
      std::is_same<c_layout,default_layout>::value;
    TEST_EQUALITY(is_b_layout_stride,true);
    TEST_EQUALITY(is_c_default_layout,true);
    TEST_EQUALITY(c.rank(),2);
    TEST_EQUALITY(c.dimension_0(),5);
    TEST_EQUALITY(c.dimension_1(),3);
    TEST_EQUALITY(dimension_scalar(b),derivative_dim_plus_one);
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize();

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  Kokkos::finalize();

  return res;
}
