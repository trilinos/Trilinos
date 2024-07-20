// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  const unsigned derivative_dim_plus_one = 7;


	// Test constructing a View pod using Kokkos::view_alloc with deduce_value_type
	{
		// Typedef View
		typedef View<double**, Kokkos::DefaultExecutionSpace> view_type;

		// Create two rank 2 Views that will be used for deducing types and Fad dims
		view_type v1("v1", 10, 4);
		view_type v2("v2", 10, 4);

		// Get common type of the Views
		using CommonValueType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::value_type;
		using ScalarArrayType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::scalar_array_type;
		// Create an instance of this returned type to pass to ViewCtorProp via view_alloc function
		auto cvt_for_ctorprop = Kokkos::common_view_alloc_prop(v1, v2);

		// Create a view with the common type and the max fad_dim of the views passed to deduce_value_type
		typedef View< CommonValueType** > ViewCommonType;
		ViewCommonType vct1( Kokkos::view_alloc("vct1", cvt_for_ctorprop), 10, 4 ); // fad_dim deduced and comes from the cvt_for_ctorprop

		TEST_EQUALITY(vct1.extent(0), v1.extent(0));
		TEST_EQUALITY(vct1.extent(1), v1.extent(1));
		TEST_EQUALITY(vct1.extent(2), v1.extent(2));
		TEST_EQUALITY( Kokkos::dimension_scalar(vct1), 0);
    bool check_eq_kokkos_type = std::is_same < CommonValueType, ScalarArrayType >::value;
    bool check_eq_scalar_double = std::is_same < double, ScalarArrayType >::value;
    TEST_EQUALITY(check_eq_kokkos_type, true);
    TEST_EQUALITY(check_eq_scalar_double, true);
	}
	// Test constructing a View of Fad using Kokkos::view_alloc with deduce_value_type
	{
		// Typedef View
		typedef View<FadType**, Kokkos::DefaultExecutionSpace> view_type;

		// Create two rank 2 Views that will be used for deducing types and Fad dims
		view_type v1("v1", 10, 4, derivative_dim_plus_one );
		view_type v2("v2", 10, 4, derivative_dim_plus_one );

		// Get common type of the Views
		using CommonValueType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::value_type;
		using ScalarArrayType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::scalar_array_type;
		// Create an instance of this returned type to pass to ViewCtorProp via view_alloc function
		auto cvt_for_ctorprop = Kokkos::common_view_alloc_prop(v1, v2);

		// Create a view with the common type and the max fad_dim of the views passed to deduce_value_type
		typedef View< CommonValueType** > ViewCommonType;
		ViewCommonType vct1( Kokkos::view_alloc("vct1", cvt_for_ctorprop), 10, 4 ); // fad_dim deduced and comes from the cvt_for_ctorprop

		TEST_EQUALITY(dimension_scalar(vct1), derivative_dim_plus_one);
		TEST_EQUALITY(vct1.extent(0), v1.extent(0));
		TEST_EQUALITY(vct1.extent(1), v1.extent(1));
		TEST_EQUALITY(vct1.extent(2), v1.extent(2));
    bool check_neq_kokkos_type = std::is_same < CommonValueType, ScalarArrayType >::value;
    bool check_eq_fad_type = std::is_same < CommonValueType, FadType >::value;
    bool check_eq_scalar_double = std::is_same < double, ScalarArrayType >::value;
    TEST_EQUALITY(check_neq_kokkos_type, false);
    TEST_EQUALITY(check_eq_fad_type, true);
    TEST_EQUALITY(check_eq_scalar_double, true);
	}
	// Test constructing a View from mix of View and Viewof Fads using Kokkos::view_alloc with deduce_value_type
	{
		// Typedef View
		typedef View<FadType**, Kokkos::DefaultExecutionSpace> view_of_fad_type;
		typedef View<double**, Kokkos::DefaultExecutionSpace> view_of_pod_type;

		// Create two rank 2 Views that will be used for deducing types and Fad dims
		view_of_fad_type v1("v1", 10, 4, derivative_dim_plus_one );
		view_of_pod_type v2("v2", 10, 4);

		// Get common type of the Views
		using CommonValueType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::value_type;
		using ScalarArrayType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::scalar_array_type;
		// Create an instance of this returned type to pass to ViewCtorProp via view_alloc function
		auto cvt_for_ctorprop = Kokkos::common_view_alloc_prop(v1, v2);

		// Create a view with the common type and the max fad_dim of the views passed to deduce_value_type
		typedef View< CommonValueType** > ViewCommonType;
		ViewCommonType vct1( Kokkos::view_alloc("vct1", cvt_for_ctorprop), 10, 4 ); // fad_dim deduced and comes from the cvt_for_ctorprop

		TEST_EQUALITY(dimension_scalar(vct1), derivative_dim_plus_one);
		TEST_EQUALITY(vct1.extent(0), v1.extent(0));
		TEST_EQUALITY(vct1.extent(1), v1.extent(1));
		TEST_EQUALITY(vct1.extent(2), v1.extent(2));
    bool check_neq_kokkos_type = std::is_same < CommonValueType, ScalarArrayType >::value;
    bool check_eq_fad_type = std::is_same < CommonValueType, FadType >::value;
    bool check_eq_scalar_double = std::is_same < double, ScalarArrayType >::value;
    TEST_EQUALITY(check_neq_kokkos_type, false);
    TEST_EQUALITY(check_eq_fad_type, true);
    TEST_EQUALITY(check_eq_scalar_double, true);
	}
	// Test constructing a DynRankView using Kokkos::view_alloc with deduce_value_type
	{
		// Typedef View
		typedef DynRankView<FadType, Kokkos::DefaultExecutionSpace> view_type;

		// Create two rank 2 Views that will be used for deducing types and Fad dims
		view_type v1("v1", 10, 4, derivative_dim_plus_one );
		view_type v2("v2", 10, 4, derivative_dim_plus_one );

		// Get common type of the Views
		using CommonValueType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::value_type;
		using ScalarArrayType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::scalar_array_type;
		// Create an instance of this returned type to pass to ViewCtorProp via view_alloc function
		auto cvt_for_ctorprop = Kokkos::common_view_alloc_prop(v1, v2);

		// Create a view with the common type and the max fad_dim of the views passed to deduce_value_type
		typedef DynRankView< CommonValueType > ViewCommonType;
		ViewCommonType vct1( Kokkos::view_alloc("vct1", cvt_for_ctorprop), 10, 4 ); // fad_dim deduced and comes from the cvt_for_ctorprop

		TEST_EQUALITY(dimension_scalar(vct1), derivative_dim_plus_one);
		TEST_EQUALITY(vct1.extent(0), v1.extent(0));
		TEST_EQUALITY(vct1.extent(1), v1.extent(1));
		TEST_EQUALITY(vct1.extent(2), v1.extent(2));
		TEST_EQUALITY(Kokkos::rank(vct1), 2);
    bool check_neq_kokkos_type = std::is_same < CommonValueType, ScalarArrayType >::value;
    bool check_eq_fad_type = std::is_same < CommonValueType, FadType >::value;
    bool check_eq_scalar_double = std::is_same < double, ScalarArrayType >::value;
    TEST_EQUALITY(check_neq_kokkos_type, false);
    TEST_EQUALITY(check_eq_fad_type, true);
    TEST_EQUALITY(check_eq_scalar_double, true);
	}
	// Test constructing a DynRankView from mix of DynRankView and DynRankView of Fads using Kokkos::view_alloc with deduce_value_type
	{
		// Typedef View
		typedef DynRankView<FadType, Kokkos::DefaultExecutionSpace> view_of_fad_type;
		typedef DynRankView<double, Kokkos::DefaultExecutionSpace> view_of_pod_type;

		// Create two rank 2 Views that will be used for deducing types and Fad dims
		view_of_fad_type v1("v1", 10, 4, derivative_dim_plus_one );
		view_of_pod_type v2("v2", 10, 4);

		// Get common type of the Views
		using CommonValueType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::value_type;
		using ScalarArrayType = typename decltype( Kokkos::common_view_alloc_prop( v1, v2 ) )::scalar_array_type;
		// Create an instance of this returned type to pass to ViewCtorProp via view_alloc function
		auto cvt_for_ctorprop = Kokkos::common_view_alloc_prop(v1, v2);

		// Create a view with the common type and the max fad_dim of the views passed to deduce_value_type
		typedef DynRankView< CommonValueType > ViewCommonType;
		ViewCommonType vct1( Kokkos::view_alloc("vct1", cvt_for_ctorprop), 10, 4 ); // fad_dim deduced and comes from the cvt_for_ctorprop

		TEST_EQUALITY(dimension_scalar(vct1), derivative_dim_plus_one);
		TEST_EQUALITY(vct1.extent(0), v1.extent(0));
		TEST_EQUALITY(vct1.extent(1), v1.extent(1));
		TEST_EQUALITY(vct1.extent(2), v1.extent(2));
		TEST_EQUALITY(Kokkos::rank(vct1), 2);
    bool check_neq_kokkos_type = std::is_same < CommonValueType, ScalarArrayType >::value;
    bool check_eq_fad_type = std::is_same < CommonValueType, FadType >::value;
    bool check_eq_scalar_double = std::is_same < double, ScalarArrayType >::value;
    TEST_EQUALITY(check_neq_kokkos_type, false);
    TEST_EQUALITY(check_eq_fad_type, true);
    TEST_EQUALITY(check_eq_scalar_double, true);
	}




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
    TEST_EQUALITY(b.extent(0),5);
    TEST_EQUALITY(b.extent(1),3);
    TEST_EQUALITY(dimension_scalar(b),0);
  }

  // Test unmanaged view of Fad
  {
    Kokkos::View<FadType*> a("a",5*3,derivative_dim_plus_one);
    using b_type = Kokkos::View<FadType**,Kokkos::MemoryUnmanaged>;
    b_type b = createViewWithType<b_type>(a,a.data(),5,3);
    TEST_EQUALITY(b.extent(0),5);
    TEST_EQUALITY(b.extent(1),3);
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
    TEST_EQUALITY(c.extent(0),5);
    TEST_EQUALITY(c.extent(1),3);
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
    TEST_EQUALITY(c.extent(0),5);
    TEST_EQUALITY(c.extent(1),3);
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
