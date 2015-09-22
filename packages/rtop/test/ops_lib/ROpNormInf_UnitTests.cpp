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


#include "RTOpPack_ROpNormInf.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, v);
  RTOpPack::ROpNormInf<Scalar> normInfOp;
  RCP<RTOpPack::ReductTarget> normInf = normInfOp.reduct_obj_create();
  normInfOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    normInf.ptr()
    );
  TEST_FLOATING_EQUALITY( normInfOp(*normInf), ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, 3, v);
  RTOpPack::ROpNormInf<Scalar> normInfOp;
  RCP<RTOpPack::ReductTarget> normInf = normInfOp.reduct_obj_create();
  normInfOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    normInf.ptr()
    );
  TEST_FLOATING_EQUALITY( normInfOp(*normInf), ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, reduct, Scalar )
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

  RTOpPack::ROpNormInf<Scalar> normInfOp;

  RCP<RTOpPack::ReductTarget> reduct1 = normInfOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = normInfOp.reduct_obj_create();

  DefaultReductTarget<ScalarMag> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarMag> >(*reduct1); 
  DefaultReductTarget<ScalarMag> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<ScalarMag> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  normInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  normInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( normInfOp(*reduct2), four );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, reduct )


} // namespace
