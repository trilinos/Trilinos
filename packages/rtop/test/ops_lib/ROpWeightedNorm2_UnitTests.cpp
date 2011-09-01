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


#include "RTOpPack_ROpWeightedNorm2.hpp"
#include "Teuchos_implicit_cast.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const ScalarMag w = as<ScalarMag>(2.0);
  const Scalar v = ST::random();
  out << "w="<<w<<"\n";
  out << "v="<<v<<"\n";
  ConstSubVectorView<Scalar> svw = newStridedSubVectorView<Scalar>(n, stride, w);
  ConstSubVectorView<Scalar> svv = newStridedSubVectorView<Scalar>(n, stride, v);
  RTOpPack::ROpWeightedNorm2<Scalar> weightedNorm2Op;
  RCP<RTOpPack::ReductTarget> weightedNorm2 = weightedNorm2Op.reduct_obj_create();
  Teuchos::implicit_ref_cast<RTOpPack::RTOpT<Scalar> >(weightedNorm2Op).apply_op(
    tuple(svw, svv)(),
    Teuchos::null,
    weightedNorm2.ptr()
    );
  TEST_FLOATING_EQUALITY( weightedNorm2Op(*weightedNorm2),
    SMT::squareroot(n)*SMT::squareroot(w)*ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack / n) );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpWeightedNorm2, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpWeightedNorm2, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpWeightedNorm2, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpWeightedNorm2, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpWeightedNorm2, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  const ScalarMag three = as<ScalarMag>(3.0);
  const ScalarMag four = as<ScalarMag>(4.0);
  const ScalarMag two = as<ScalarMag>(2.0);

  RTOpPack::ROpWeightedNorm2<Scalar> weightedNorm2Op;

  RCP<RTOpPack::ReductTarget> reduct1 = weightedNorm2Op.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = weightedNorm2Op.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  weightedNorm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  weightedNorm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( weightedNorm2Op(*reduct2), SMT::squareroot(three+four+two),
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpWeightedNorm2, reduct )


} // namespace
