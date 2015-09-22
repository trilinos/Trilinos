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


#include "RTOpPack_ROpCountNanInf.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


using TestingSupportHelpers::print;


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(4, stride, as<Scalar>(0.0));

  sv(0) = ST::nan();
  sv(1) = ST::random();
  sv(2) = ST::one()/ST::zero();
  sv(3) = ST::zero();

  TEST_ASSERT(ST::isnaninf(sv(0)));
  TEST_ASSERT(!ST::isnaninf(sv(1)));
  TEST_ASSERT(ST::isnaninf(sv(2)));
  TEST_ASSERT(!ST::isnaninf(sv(3)));

  print(sv, "sv", out);

  RTOpPack::ROpCountNanInf<Scalar> countNanInfOp;
  RCP<RTOpPack::ReductTarget> countNanInf = countNanInfOp.reduct_obj_create();
  Teuchos::implicit_ref_cast<RTOpPack::RTOpT<Scalar> >(countNanInfOp).apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    countNanInf.ptr()
    );

  const index_type countNanInf_val = countNanInfOp(*countNanInf);
  TEST_EQUALITY( countNanInf_val, 2 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpCountNanInf, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpCountNanInf, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpCountNanInf, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpCountNanInf, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpCountNanInf, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::ROpCountNanInf<Scalar> countNanInfOp;

  RCP<ReductTarget> reduct1 = countNanInfOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = countNanInfOp.reduct_obj_create();

  DefaultReductTarget<index_type> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<index_type> >(*reduct1); 

  scalarReduct1.set(1);
  countNanInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(0);
  countNanInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(2);
  countNanInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( countNanInfOp(*reduct2), as<index_type>(3) )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpCountNanInf, reduct )


} // namespace
