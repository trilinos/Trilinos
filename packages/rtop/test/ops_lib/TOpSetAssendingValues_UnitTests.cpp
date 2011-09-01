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


#include "RTOpPack_TOpSetAssendingValues.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"
#include "Teuchos_implicit_cast.hpp"

// Must come last!
#include "supportUnitTestsHelpers.hpp"


namespace RTOpPack {


const double g_tol = 100.0 * ScalarTraits<double>::eps();

template<class Scalar>
void basicTest(const int size, const int stride, const Scalar offset,
  FancyOStream &out, bool &success)
{
  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::tuple;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag ;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(size, stride, ST::zero());

  TOpSetAssendingValues<Scalar> op(offset);

  op.apply_op(null, tuple<SubVectorView<Scalar> >(sv)(), null);

  Scalar sum = ST::zero();
  for (int i = 0; i < size; ++i)
    sum += sv[i];

  //const Scalar ssize = size;
    
  TEST_FLOATING_EQUALITY( sum, offset*as<Scalar>(size) + as<Scalar>(0.5 * (size + 1) * size),
    as<ScalarMag>(g_tol * size) );
  
}

//
// size=1, offset=0
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetAssendingValues, basicTest_1_1_0, Scalar )
{
  basicTest<Scalar>(1, 1, 0.0, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetAssendingValues, basicTest_1_1_0)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetAssendingValues, basicTest_1_3_0, Scalar )
{
  basicTest<Scalar>(1, 3, 0.0, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetAssendingValues, basicTest_1_3_0)

//
// size=2, offset=5.0
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetAssendingValues, basicTest_2_1_5, Scalar )
{
  basicTest<Scalar>(2, 1, 5.0, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetAssendingValues, basicTest_2_1_5)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetAssendingValues, basicTest_2_3_5, Scalar )
{
  basicTest<Scalar>(2, 3, 5.0, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetAssendingValues, basicTest_2_3_5)


//
// size=4, offset=5.0
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetAssendingValues, basicTest_4_1_5, Scalar )
{
  basicTest<Scalar>(4, 1, 5.0, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetAssendingValues, basicTest_4_1_5)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetAssendingValues, basicTest_4_3_5, Scalar )
{
  basicTest<Scalar>(4, 3, 5.0, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetAssendingValues, basicTest_4_3_5)


} // namespace RTOpPack
