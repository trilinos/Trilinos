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


#include "RTOpPack_TOpPairWiseMax.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const int local_n = 4;

  SubVectorView<Scalar> v0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  SubVectorView<Scalar> v1 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  const Scalar alpha = 1.5*ST::one();

  SubVectorView<Scalar> expected_z0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  v0[0] = 1.0;  v1[0] = 2.0;  expected_z0[0] = 3.0;
  v0[1] = -1.0; v1[1] = -2.0; expected_z0[1] = -1.5;
  v0[2] = 1.0;  v1[2] = -2.0; expected_z0[2] = 1.5;
  v0[3] = 2.0;  v1[3] = 1.0;  expected_z0[3] = 3.0;
  
  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());

  RTOpPack::TOpPairWiseMax<Scalar> op(alpha);
  op.apply_op( tuple<ConstSubVectorView<Scalar>>(v0,v1), tuple(z0)(), null );

  if (verbose) {
    out << "alpha = " << alpha << "\n";
    dumpSubVectorView(v0, "v0", out);
    dumpSubVectorView(v1, "v1", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0),
    as<ScalarMag>(ST::eps() * errorTolSlack)
    );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpPairWiseMax, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TOpPairWiseMax, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpPairWiseMax, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TOpPairWiseMax, nonunitStride )


} // namespace
