// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
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

  template<typename Scalar>
  void runTensorViewFunctorTest(ViewType<Scalar> tensor_expected, ViewType<Scalar> view1, ViewType<Scalar> view2, double weight, bool tensorPoints,
                                Teuchos::FancyOStream &out, bool &success)
  {
    double tol = 1e-15;
    using namespace Intrepid2;
    
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using ScalarViewType = ViewType<Scalar>;
    
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
    Kokkos::parallel_for( policy , functor, "TensorViewFunctor");
    
    auto tensor_actual_host   = getHostCopy(tensor_actual);
    auto tensor_expected_host = getHostCopy(tensor_expected);
    
    using ViewIterator = ViewIterator<ScalarViewType, Scalar>;
    ViewIterator it_actual(tensor_actual_host);
    ViewIterator it_expected(tensor_expected_host);
    
    do
    {
      auto actual_value   = it_actual.get();
      auto expected_value = it_expected.get();
      
      if (!approximatelyEqual(actual_value, expected_value, tol))
      {
        success = false;
        std::cout << "FAILURE: In entry " << it_actual.getEnumerationIndex() << ", ";
        std::cout << "actual (" << actual_value << ") differs from expected (" << expected_value << ")";
        std::cout << " by " << std::abs(actual_value-expected_value) << std::endl;
      }
      
    } while ((it_actual.increment() >= 0) && (it_expected.increment() >= 0));
  }
  
  template<typename Scalar>
  void runTensorViewFunctorTests(Teuchos::FancyOStream &out, bool &success)
  {
    using ScalarViewType = ViewType<Scalar>;
    
    // TEST 1: simple contraction
    // we'll use trivial fields so as to factor out problems in the tensor product logic
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
    weight = 2.0;
    for (int i=0; i<num_fields; i++)
    {
      tensor_expected(i,0) *= weight;
    }
    runTensorViewFunctorTest(tensor_expected, view1, view2, weight, tensor_points, out, success);
    
    // TEST 4: scalar times vector
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
    num_fields1 = 2;
    num_fields2 = 2;
    num_fields = num_fields1 * num_fields2;
    num_points = 2;
    Kokkos::resize(tensor_expected_host,num_fields, num_points,space_dim);
    Kokkos::resize(view1_host,          num_fields1,num_points);
    Kokkos::resize(view2_host,          num_fields2,num_points,space_dim);
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
