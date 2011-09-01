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


#include "RTOpPack_TOpLinearCombination.hpp"
#include "Teuchos_Utils.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTestGuts(
  const Teuchos::ArrayView<const Scalar> &alpha,
  const Scalar &beta,
  const Teuchos::ArrayView<const ConstSubVectorView<Scalar> > &v,
  const ConstSubVectorView<Scalar> &orig_z0,
  FancyOStream &out, bool &success)
{

  typedef Teuchos::Utils TU;
  using Teuchos::as;
  using Teuchos::Array;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  Teuchos::OSTab tab(out);

  out << "\nTesting TOpLinearCombination for:\n"
      << "  alpha = " << alpha << "\n"
      << "  beta = " << beta << "\n";

  Teuchos::OSTab tab2(out);
  
  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(orig_z0.subDim(), orig_z0.stride(), ST::nan());
  RTOpPack::assign_entries<Scalar>( Teuchos::outArg(z0), orig_z0 );

  const RTOpPack::TOpLinearCombination<Scalar> op(alpha, beta);
  op.apply_op( v, tuple(z0)(), null );
  
  SubVectorView<Scalar> expected_z0 =
    newStridedSubVectorView<Scalar>(orig_z0.subDim(), orig_z0.stride(), ST::nan());
  RTOpPack::assign_entries<Scalar>( Teuchos::outArg(expected_z0), orig_z0 );

  const index_type subDim = z0.subDim();

  // expected_z0 *= beta
  if( beta == ST::zero() ) {
    for( int j = 0; j < subDim; ++j )
      expected_z0[j] = ST::zero();
  }
  else if( beta != ST::one() ) {
    for( int j = 0; j < subDim; ++j )
      expected_z0[j] *= beta;
  }
  // expected_z0 += sum( alpha[k]*v[k], k=0...num_vecs-1)
  for( int j = 0; j < subDim; ++j ) {
    for( int k = 0; k < as<int>(v.size()); ++k ) {
      expected_z0[j] += alpha[k] * v[k][j];
    }
  }

  if (verbose) {
    for (int k = 0; k < as<int>(v.size()); ++k) {
      dumpSubVectorView(v[k], "v"+TU::toString(k), out);
    }
    dumpSubVectorView(orig_z0, "orig_z0", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0),
    as<ScalarMag>(ST::eps() * errorTolSlack) * as<ScalarMag>(v.size())
    );

}


template<class Scalar>
void basicTest(const int stride, const int num_sub_vecs, FancyOStream &out,
  bool &success)
{

  using Teuchos::Array;
  typedef ScalarTraits<Scalar> ST;

  out << "\nTesting TOpLinearCombination for num_sub_vecs = " << num_sub_vecs << "\n";

  Array<ConstSubVectorView<Scalar> > v;
  for (int k = 0; k < num_sub_vecs; ++k) {
    v.push_back(newStridedRandomSubVectorView<Scalar>(n, stride));
  }

  ConstSubVectorView<Scalar> orig_z0 =
    newStridedRandomSubVectorView<Scalar>(n, stride);

  const int num_alpha_vals = 3;
  Array<Scalar> alpha_vals = tuple(ST::zero(), ST::one(), ST::random());
  Array<Scalar> alpha(num_sub_vecs);
  Scalar beta = ST::nan();

  for (int j = 0; j < num_alpha_vals; ++j) {
    beta = alpha_vals[j];
    for (int k0 = 0; k0 < num_alpha_vals; ++k0)
    {
      alpha[0] = alpha_vals[k0];
      if (num_sub_vecs > 1) {
        for (int k1 = 0; k1 < num_alpha_vals; ++k1)
        {
          alpha[1] = alpha_vals[k1];
          if (num_sub_vecs > 2) {
            for (int k2 = 0; k2 < num_alpha_vals; ++k2)
            {
              alpha[2] = alpha_vals[k2];
              if (num_sub_vecs > 3) {
                TEST_FOR_EXCEPT(true); // This is not scalable!
              }
              else {
                basicTestGuts<Scalar>(alpha(), beta, v(), orig_z0, out, success);
              }
            }
          }
          else {
            basicTestGuts<Scalar>(alpha(), beta, v(), orig_z0, out, success);
          }
        }
      }
      else {
        basicTestGuts<Scalar>(alpha(), beta, v(), orig_z0, out, success);
      }
    }
  }


}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpLinearCombination, unitStride_1, Scalar )
{
  basicTest<Scalar>(1, 1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpLinearCombination, unitStride_1 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpLinearCombination, nonunitStride_1, Scalar )
{
  basicTest<Scalar>(4, 1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpLinearCombination, nonunitStride_1 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpLinearCombination, unitStride_2, Scalar )
{
  basicTest<Scalar>(1, 2, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpLinearCombination, unitStride_2 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpLinearCombination, nonunitStride_2, Scalar )
{
  basicTest<Scalar>(4, 2, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpLinearCombination, nonunitStride_2 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpLinearCombination, unitStride_3, Scalar )
{
  basicTest<Scalar>(1, 3, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpLinearCombination, unitStride_3 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpLinearCombination, nonunitStride_3, Scalar )
{
  basicTest<Scalar>(4, 3, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpLinearCombination, nonunitStride_3 )


} // namespace
