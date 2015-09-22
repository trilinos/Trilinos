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


#include "RTOpPack_ROpNorm2.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  ConstSubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, stride, v);
  RTOpPack::ROpNorm2<Scalar> norm2Op;
  RCP<RTOpPack::ReductTarget> norm2 = norm2Op.reduct_obj_create();
  norm2Op.apply_op(
    tuple(sv)(),
    Teuchos::null,
    norm2.ptr()
    );
  TEST_FLOATING_EQUALITY( norm2Op(*norm2),
    SMT::squareroot(n)*ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, reduct, Scalar )
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

  RTOpPack::ROpNorm2<Scalar> norm2Op;

  RCP<RTOpPack::ReductTarget> reduct1 = norm2Op.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = norm2Op.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  norm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  norm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( norm2Op(*reduct2), SMT::squareroot(three+four+two),
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, reduct )


} // namespace
