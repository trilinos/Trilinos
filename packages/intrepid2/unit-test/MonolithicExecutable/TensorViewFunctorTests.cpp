// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   TensorViewFunctorTests.cpp
    \brief  Tests to verify the TensorViewFunctor used in TensorBasis.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_TensorBasis.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;

  template<typename ScalarViewType>
  void runTensorViewFunctorTest(ScalarViewType tensor_expected, ScalarViewType view1, ScalarViewType view2, double weight, bool tensorPoints,
                                Teuchos::FancyOStream &out, bool &success)
  {
    double tol = 1e-15;
    using namespace Intrepid2;
    
    using ExecutionSpace = typename ScalarViewType::execution_space;
    using Scalar         = typename ScalarViewType::value_type;
    
    const bool hasADType = false;
    const int vectorSize = hasADType ? FAD_VECTOR_SIZE : VECTOR_SIZE;
    
    auto policy = Kokkos::TeamPolicy<ExecutionSpace>(view1.extent_int(0),Kokkos::AUTO(),vectorSize);
    
    using FunctorType = TensorViewFunctor<ExecutionSpace, Scalar, ScalarViewType>;
    
    ScalarViewType tensor_actual;
    
    // TODO: figure out a better way to initialize tensor_actual to be the same size/shape as tensor_expected
    if (tensor_expected.rank() == 1)
    {
      tensor_actual = ScalarViewType("tensor actual", tensor_expected.extent_int(0));
    }
    else if (tensor_expected.rank() == 2)
    {
      tensor_actual = ScalarViewType("tensor actual", tensor_expected.extent_int(0),tensor_expected.extent_int(1));
    }
    else if (tensor_expected.rank() == 3)
    {
      tensor_actual = ScalarViewType("tensor actual", tensor_expected.extent_int(0),tensor_expected.extent_int(1), tensor_expected.extent_int(2));
    }
    else if (tensor_expected.rank() == 4)
    {
      tensor_actual = ScalarViewType("tensor actual", tensor_expected.extent_int(0),tensor_expected.extent_int(1), tensor_expected.extent_int(2), tensor_expected.extent_int(2));
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Test does not yet support this output rank");
    }
    
    FunctorType functor(tensor_actual, view1, view2, tensorPoints, weight);
    Kokkos::parallel_for("TensorViewFunctor", policy, functor);
    
    switch (tensor_expected.rank())
    {
      case 1: testFloatingEquality1(tensor_actual,tensor_expected,tol,tol,out,success); break;
      case 2: testFloatingEquality2(tensor_actual,tensor_expected,tol,tol,out,success); break;
      case 3: testFloatingEquality3(tensor_actual,tensor_expected,tol,tol,out,success); break;
      case 4: testFloatingEquality4(tensor_actual,tensor_expected,tol,tol,out,success); break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Test does not yet support this output rank");
    }
  }
  
  template<typename Scalar>
  void runTensorViewFunctorTests(Teuchos::FancyOStream &out, bool &success)
  {
    using DeviceType = DefaultTestDeviceType;
    using ScalarViewType = ViewType<Scalar,DeviceType>;
    
    // TEST 1: simple contraction
    // we'll use trivial fields so as to factor out problems in the tensor product logic
    out << "TEST 1: simple contraction.\n";
    int num_fields1 = 1;
    int num_fields2 = 1;
    int num_fields = num_fields1 * num_fields2;
    int num_points = 1;
    int space_dim = 3;
    ScalarViewType tensor_expected("expected_tensor",num_fields,num_points);
    ScalarViewType view1("view1",num_fields,num_points,space_dim);
    ScalarViewType view2("view2",num_fields,num_points,space_dim);
    
    auto tensor_expected_host = getHostCopy(tensor_expected);
    auto view1_host           = getHostCopy(view1);
    auto view2_host           = getHostCopy(view2);
    
    view1_host(0,0,0) = 3.0;
    view1_host(0,0,1) = 2.0;
    view1_host(0,0,2) = 0.0;
    view2_host(0,0,0) = 1.0;
    view2_host(0,0,1) = 0.5;
    view2_host(0,0,2) = 3.0;
    
    tensor_expected_host(0,0,0) = 0;
    for (int d=0; d<space_dim; d++)
    {
      tensor_expected_host(0,0,0) += view1_host(0,0,d) * view2_host(0,0,d);
    }
    Kokkos::deep_copy(view1,           view1_host);
    Kokkos::deep_copy(view2,           view2_host);
    Kokkos::deep_copy(tensor_expected, tensor_expected_host);
    
    double weight = 1.0;
    bool tensor_points = false; // does not matter for the single-point case
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
    
    // TEST 2: tensor product ordering
    out << "TEST 2: tensor product ordering.\n";
    num_fields1 = 2;
    num_fields2 = 2;
    num_fields = num_fields1 * num_fields2;
    num_points = 1;
    Kokkos::resize(tensor_expected,num_fields,num_points);
    Kokkos::resize(view1,          num_fields1,num_points);
    Kokkos::resize(view2,          num_fields2,num_points);
    
    Kokkos::resize(tensor_expected_host,num_fields,num_points);
    Kokkos::resize(view1_host,          num_fields1,num_points);
    Kokkos::resize(view2_host,          num_fields2,num_points);
    
    view1_host(0,0) = 3.0;
    view1_host(1,0) = 2.0;
    view2_host(0,0) = 1.0;
    view2_host(1,0) = 0.5;
    // view1 is the fastest-moving: tensor entry 1 corresponds to view1(1,0) and view2(0,0)
    tensor_expected_host(0,0) = view1_host(0,0) * view2_host(0,0);
    tensor_expected_host(1,0) = view1_host(1,0) * view2_host(0,0);
    tensor_expected_host(2,0) = view1_host(0,0) * view2_host(1,0);
    tensor_expected_host(3,0) = view1_host(1,0) * view2_host(1,0);
    
    Kokkos::deep_copy(view1,           view1_host);
    Kokkos::deep_copy(view2,           view2_host);
    Kokkos::deep_copy(tensor_expected, tensor_expected_host);
    
    weight = 1.0;
    tensor_points = false; // does not matter for the single-point case
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
    
    // TEST 3: like TEST 2, but include non-trivial weight
    out << "TEST 3: like TEST 2, but include non-trivial weight.\n";
    weight = 2.0;
    for (int i=0; i<num_fields; i++)
    {
      tensor_expected_host(i,0) *= weight;
    }
    Kokkos::deep_copy(tensor_expected, tensor_expected_host);
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
    
    // TEST 4: scalar times vector
    out << "TEST 4: scalar times vector.\n";
    num_fields1 = 2;
    num_fields2 = 2;
    num_fields = num_fields1 * num_fields2;
    num_points = 1;
    Kokkos::resize(tensor_expected,num_fields,  num_points, space_dim);
    Kokkos::resize(view1,          num_fields1, num_points);
    Kokkos::resize(view2,          num_fields2, num_points, space_dim);
    
    Kokkos::resize(tensor_expected_host, num_fields,  num_points,space_dim);
    Kokkos::resize(view1_host,           num_fields1, num_points);
    Kokkos::resize(view2_host,           num_fields2, num_points,space_dim);
    
    view1_host(0,0) = 3.0;
    view1_host(1,0) = 2.0;
    view2_host(0,0,0) = 1.0;
    view2_host(1,0,0) = 0.5;
    view2_host(0,0,1) = 2.0;
    view2_host(1,0,1) = 1.0;
    view2_host(0,0,2) = 4.0;
    view2_host(1,0,2) = 2.0;
    for (int d=0; d<space_dim; d++)
    {
      // view1 is the fastest-moving: tensor entry 1 corresponds to view1(1,…) and view2(0,…)
      tensor_expected_host(0,0,d) = view1_host(0,0) * view2_host(0,0,d);
      tensor_expected_host(1,0,d) = view1_host(1,0) * view2_host(0,0,d);
      tensor_expected_host(2,0,d) = view1_host(0,0) * view2_host(1,0,d);
      tensor_expected_host(3,0,d) = view1_host(1,0) * view2_host(1,0,d);
    }
    
    Kokkos::deep_copy(view1,           view1_host);
    Kokkos::deep_copy(view2,           view2_host);
    Kokkos::deep_copy(tensor_expected, tensor_expected_host);
    
    weight = 1.0;
    tensor_points = false; // does not matter for the single-point case
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
    
    // TEST 5: scalar times vector, nontrivial points, but still matching in point dimension
    out << "TEST 5: scalar times vector, nontrivial points, but still matching in point dimension.\n";
    num_fields1 = 2;
    num_fields2 = 2;
    num_fields = num_fields1 * num_fields2;
    num_points = 2;
    
    Kokkos::resize(tensor_expected,      num_fields,  num_points, space_dim);
    Kokkos::resize(view1,                num_fields1, num_points);
    Kokkos::resize(view2,                num_fields2, num_points, space_dim);
    
    Kokkos::resize(tensor_expected_host, num_fields,  num_points, space_dim);
    Kokkos::resize(view1_host,           num_fields1, num_points);
    Kokkos::resize(view2_host,           num_fields2, num_points, space_dim);
    
    view1_host(0,0)   = 3.0;
    view1_host(1,0)   = 2.0;
    view2_host(0,0,0) = 1.0;
    view2_host(1,0,0) = 0.5;
    view2_host(0,0,1) = 2.0;
    view2_host(1,0,1) = 1.0;
    view2_host(0,0,2) = 4.0;
    view2_host(1,0,2) = 2.0;
    view1_host(0,1)   = 1.0;
    view1_host(1,1)   = 1.0;
    view2_host(0,1,0) = 1.0;
    view2_host(1,1,0) = 2.0;
    view2_host(0,1,1) = 3.0;
    view2_host(1,1,1) = 4.0;
    view2_host(0,1,2) = 5.0;
    view2_host(1,1,2) = 6.0;
    for (int point_ordinal=0; point_ordinal<num_points; point_ordinal++)
    {
      for (int d=0; d<space_dim; d++)
      {
        // view1 is the fastest-moving: tensor entry 1 corresponds to view1(1,…) and view2(0,…)
        tensor_expected_host(0,point_ordinal,d) = view1_host(0,point_ordinal) * view2_host(0,point_ordinal,d);
        tensor_expected_host(1,point_ordinal,d) = view1_host(1,point_ordinal) * view2_host(0,point_ordinal,d);
        tensor_expected_host(2,point_ordinal,d) = view1_host(0,point_ordinal) * view2_host(1,point_ordinal,d);
        tensor_expected_host(3,point_ordinal,d) = view1_host(1,point_ordinal) * view2_host(1,point_ordinal,d);
      }
    }
    
    Kokkos::deep_copy(view1,           view1_host);
    Kokkos::deep_copy(view2,           view2_host);
    Kokkos::deep_copy(tensor_expected, tensor_expected_host);
    
    weight = 1.0;
    tensor_points = false; // does not matter for the single-point case
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
    
    // TEST 6: like TEST 2 above, but with different field counts
    out << "TEST 6: like TEST 2 above, but with different field counts.\n";
    num_fields1 = 2;
    num_fields2 = 3;
    num_fields = num_fields1 * num_fields2;
    num_points = 1;
    Kokkos::resize(tensor_expected,num_fields,num_points);
    Kokkos::resize(view1,          num_fields1,num_points);
    Kokkos::resize(view2,          num_fields2,num_points);
    
    Kokkos::resize(tensor_expected_host, num_fields,  num_points);
    Kokkos::resize(view1_host,           num_fields1, num_points);
    Kokkos::resize(view2_host,           num_fields2, num_points);
    
    view1_host(0,0) = 3.0;
    view1_host(1,0) = 2.0;
    view2_host(0,0) = 1.0;
    view2_host(1,0) = 0.5;
    view2_host(2,0) = 1.5;
    // view1 is the fastest-moving: tensor entry 1 corresponds to view1(1,0) and view2(0,0)
    tensor_expected_host(0,0) = view1_host(0,0) * view2_host(0,0);
    tensor_expected_host(1,0) = view1_host(1,0) * view2_host(0,0);
    tensor_expected_host(2,0) = view1_host(0,0) * view2_host(1,0);
    tensor_expected_host(3,0) = view1_host(1,0) * view2_host(1,0);
    tensor_expected_host(4,0) = view1_host(0,0) * view2_host(2,0);
    tensor_expected_host(5,0) = view1_host(1,0) * view2_host(2,0);
    
    Kokkos::deep_copy(view1,           view1_host);
    Kokkos::deep_copy(view2,           view2_host);
    Kokkos::deep_copy(tensor_expected, tensor_expected_host);
    
    weight = 1.0;
    tensor_points = false; // does not matter for the single-point case
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
  }
  
  TEUCHOS_UNIT_TEST( TensorViewFunctor, MultipleTests )
  {
    using Scalar = double;
    runTensorViewFunctorTests<Scalar>(out, success);
  }
} // namespace
