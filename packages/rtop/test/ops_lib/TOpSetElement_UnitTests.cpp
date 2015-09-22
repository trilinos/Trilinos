/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "RTOpPack_TOpSetElement.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"
#include "Teuchos_implicit_cast.hpp"

// Must come last!
#include "supportUnitTestsHelpers.hpp"


namespace RTOpPack {


const double g_tol = 100.0 * ScalarTraits<double>::eps();


template<class Scalar>
void setElementTestCase( TOpSetElement<Scalar> &setEleOp, 
  const SubVectorView<Scalar> sv, const Ordinal i,
  FancyOStream &out, bool &success)
{

  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::tuple;
  typedef ScalarTraits<Scalar> ST;

  out << "\nTest i="<<i<<":\n";

  const Scalar val_i = as<Scalar>(i+1);
  setEleOp.initialize(i, val_i);
  TEST_EQUALITY_CONST(setEleOp.range().lbound(), i);
  TEST_EQUALITY_CONST(setEleOp.range().ubound(), i);

  std::fill_n(sv.values().begin(), sv.values().size(), ST::zero());

  setEleOp.apply_op(null, tuple<SubVectorView<Scalar> >(sv)(), null);
  TEUCHOS_TEST_ARRAY_ELE_EQUALITY(sv, i, val_i, true, out, success);

  ROpSum<Scalar> sumOp;
  const RCP<ReductTarget> reduct_obj = sumOp.reduct_obj_create();
  sumOp.apply_op(tuple<ConstSubVectorView<Scalar> >(sv)(), null, reduct_obj.ptr());
  TEST_EQUALITY_CONST(sumOp(*reduct_obj), val_i);

}


template<class Scalar>
void basicTestTOpSetElement(const int size, const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::tuple;
  typedef ScalarTraits<Scalar> ST;
  
  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(size, stride, ST::zero());

  TOpSetElement<Scalar> setEleOp;

  setElementTestCase<Scalar>(setEleOp, sv, 0, out, success);

  setElementTestCase<Scalar>(setEleOp, sv, 1, out, success);

  setElementTestCase<Scalar>(setEleOp, sv, size-2, out, success);

  setElementTestCase<Scalar>(setEleOp, sv, size-1, out, success);

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetElement, basicTest_4_1, Scalar )
{
  basicTestTOpSetElement<Scalar>(4, 1, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetElement, basicTest_4_1)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetElement, basicTest_4_3, Scalar )
{
  basicTestTOpSetElement<Scalar>(4, 3, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetElement, basicTest_4_3)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetElement, basicTest_6_1, Scalar )
{
  basicTestTOpSetElement<Scalar>(6, 1, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetElement, basicTest_6_1)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetElement, basicTest_6_3, Scalar )
{
  basicTestTOpSetElement<Scalar>(6, 3, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetElement, basicTest_6_3)


} // namespace RTOpPack
