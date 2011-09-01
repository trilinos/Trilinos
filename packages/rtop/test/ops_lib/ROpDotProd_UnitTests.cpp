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


#include "RTOpPack_ROpDotProd.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";
  SubVectorView<Scalar> sv1 = newSubVectorView<Scalar>(n, v1);
  SubVectorView<Scalar> sv2 = newSubVectorView<Scalar>(n, v2);
  RTOpPack::ROpDotProd<Scalar> dotProdOp;
  RCP<RTOpPack::ReductTarget> dotProd = dotProdOp.reduct_obj_create();
  dotProdOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv1, sv2)(),
    Teuchos::null,
    dotProd.ptr()
    );
  TEST_FLOATING_EQUALITY( dotProdOp(*dotProd), ST::conjugate(v1)*v2*as<Scalar>(n),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";
  SubVectorView<Scalar> sv1 = newStridedSubVectorView<Scalar>(n, 2, v1);
  SubVectorView<Scalar> sv2 = newStridedSubVectorView<Scalar>(n, 3, v2);
  RTOpPack::ROpDotProd<Scalar> dotProdOp;
  RCP<RTOpPack::ReductTarget> dotProd = dotProdOp.reduct_obj_create();
  dotProdOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv1, sv2)(),
    Teuchos::null,
    dotProd.ptr()
    );
  TEST_FLOATING_EQUALITY( dotProdOp(*dotProd), ST::conjugate(v1)*v2*as<Scalar>(n),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";

  RTOpPack::ROpDotProd<Scalar> dotProdOp;

  RCP<RTOpPack::ReductTarget> reduct1 = dotProdOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = dotProdOp.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(v1);
  scalarReduct2.set(v2);

  dotProdOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( scalarReduct2.get(), v1+v2, as<ScalarMag>(ST::eps()*errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, reduct )


} // namespace
