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

#include "RTOpPack_RTOpSubRangeDecorator.hpp"
#include "RTOpPack_TOpSetAssendingValues.hpp"
#include "RTOpPack_TOpAssignVectors.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "Teuchos_as.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace RTOpPack {


// ////////////////////
// Utilities


template<class Scalar>
SubVectorView<Scalar>
createSubVectorView(const int n)
{
  const ArrayRCP<Scalar> values = Teuchos::arcp<Scalar>(n);
  std::fill_n(values.begin(), n, ScalarTraits<Scalar>::zero());
  return SubVectorView<Scalar>(0, n, values, 1);
}


template<class Scalar>
SubVectorView<Scalar>
createFilledSubVectorView(const int n, const Scalar lowerVal = static_cast<Scalar>(1.0))
{
  using Teuchos::as;
  SubVectorView<Scalar> sv = createSubVectorView<Scalar>(n);
  for (int i = 0; i < n; ++i)
    sv[i] = lowerVal + as<Scalar>(i);
  return sv;
}


template<class Scalar>
Scalar computeSum(const ConstSubVectorView<Scalar> &sv)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar sum = ST::zero();
  for (int i = 0; i < sv.subDim(); ++i) {
    sum += sv[i];
  }
  return sum;
}


template<class Scalar>
Scalar sumIntegers(const int n)
{
  return (0.5 * (n+1) * n);
}


// ////////////////////////////////////
// Test construction


//
// Test default constructor
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  constructDefault, Scalar )
{
  RTOpSubRangeDecorator<Scalar> srop;
  TEST_EQUALITY_CONST(srop.getNonconstOp(), Teuchos::null);
  TEST_EQUALITY_CONST(srop.getOp(), Teuchos::null);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  constructDefault )


//
// Test consturction from a TOp
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  constructTOp, Scalar )
{
  using Teuchos::outArg;
  TOpAssignVectors<Scalar> vecAssignOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(vecAssignOp)));
  TEST_ASSERT(nonnull(op.getNonconstOp()));
  TEST_ASSERT(nonnull(op.getOp()));
  int numValues = -1, numIndexes = -1, numChars = -1;
  op.get_reduct_type_num_entries(outArg(numValues), outArg(numIndexes), outArg(numChars));
  TEST_EQUALITY_CONST(numValues, 0);
  TEST_EQUALITY_CONST(numIndexes, 0);
  TEST_EQUALITY_CONST(numChars, 0);
  TEST_ASSERT(is_null(op.reduct_obj_create()));
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  constructTOp )


//
// Test consturction from a ROp
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  constructROp, Scalar )
{
  using Teuchos::null;
  using Teuchos::as;
  using Teuchos::outArg;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Tuple;
  using Teuchos::tuple;

  ROpSum<Scalar> sumOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(sumOp)));

  TEST_ASSERT(nonnull(op.getNonconstOp()));
  TEST_ASSERT(nonnull(op.getOp()));

  TEST_EQUALITY_CONST(op.op_name(),
    std::string("RTOpSubRangeDecorator{")+sumOp.op_name()+"}");

  TEST_ASSERT(op.coord_invariant());

  int numValues = -1, numIndexes = -1, numChars = -1;
  op.get_reduct_type_num_entries(outArg(numValues), outArg(numIndexes), outArg(numChars));
  TEST_EQUALITY_CONST(numValues, 1);
  TEST_EQUALITY_CONST(numIndexes, 0);
  TEST_EQUALITY_CONST(numChars, 0);

  const RCP<ReductTarget> reduct_obj = op.reduct_obj_create();
  TEST_ASSERT(nonnull(reduct_obj));

  const RCP<DefaultReductTarget<Scalar> > scalar_reduct_obj =
    rcp_dynamic_cast<DefaultReductTarget<Scalar> >(reduct_obj, true);
  TEST_ASSERT(nonnull(scalar_reduct_obj));

  const RCP<ReductTarget> reduct_obj_2 = op.reduct_obj_create();
  const RCP<DefaultReductTarget<Scalar> > scalar_reduct_obj_2 =
    rcp_dynamic_cast<DefaultReductTarget<Scalar> >(reduct_obj_2, true);
  const Scalar one = as<Scalar>(1.0), two = as<Scalar>(2.0), three = as<Scalar>(3.0);
  scalar_reduct_obj->set(one);
  scalar_reduct_obj_2->set(two);
  op.reduce_reduct_objs(*reduct_obj, reduct_obj_2.ptr());
  TEST_EQUALITY_CONST( scalar_reduct_obj_2->get(), three );

  Tuple<Scalar,1> valueData;
  op.extract_reduct_obj_state(*reduct_obj_2, valueData(), null, null);
  TEST_EQUALITY_CONST( valueData[0], three );

  op.reduct_obj_reinit(reduct_obj_2.ptr());
  const Scalar zero = as<Scalar>(0.0);
  TEST_EQUALITY_CONST( scalar_reduct_obj_2->get(), zero );

  op.load_reduct_obj_state( tuple<Scalar>(one)(), null, null, reduct_obj_2.ptr());
  TEST_EQUALITY_CONST( scalar_reduct_obj_2->get(), one );

}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( RTOpSubRangeDecorator,
  constructROp )


//
// Test non-const construction
//

//
// Test const construction
//


// ///////////////////////////////////////
// Test usage
//
// Here we need to test several different use cases and different types of
// ranges to make sure that this class works in all cases.  We need to test
// the following scenarios:
//
//  0) Equal matching ranges (fall through):
//     subvecs  [......]
//     rtop     [......]
//
//  1) Equal matching ranges:
//     subvecs  [......]
//     rtop     [......]
//
//  2) Match first part of range (inner subvec):
//     subvecs  [..]
//     rtop     [......]
//
//  3) Match middle part of range (inner subvec):
//     subvecs    [..]
//     rtop     [......]
//
//  4) Match last part of range (inner subvec):
//     subvecs      [..]
//     rtop     [......]
//
//  5) Match first part of range (inner rtop):
//     subvecs  [......]
//     rtop     [..]
//
//  6) Match middle part of range (inner rtop):
//     subvecs  [......]
//     rtop       [..]
//
//  7) Match last part of range (inner rtop):
//     subvecs  [......]
//     rtop         [..]
//
//  8) Non-overlapping subvec first
//     subvecs  [......]
//     rtop         [......]
//
//  9) Non-overlapping rtop first
//     subvecs      [......]
//     rtop     [......]
//
//  10) Zero-len subvecs, first
//     subvecs  []
//     rtop     [......]
//
//  11) Zero-len subvecs, middle
//     subvecs     []
//     rtop     [......]
//
//  12) Zero-len subvecs, end
//     subvecs        []
//     rtop     [......]
//


//
//  0) Equal matching ranges (fall through):
//     subvecs  [......]
//     rtop     [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  equalRanges_fallthrough, Scalar )
{
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;

  TOpAssignVectors<Scalar> vecAssignOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(vecAssignOp)));

  ConstSubVectorView<Scalar> x = createFilledSubVectorView<Scalar>(n);
  TEST_EQUALITY(computeSum<Scalar>(x), sumIntegers<Scalar>(n));
  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  op.apply_op(tuple(x), tuple(y), null);
  TEST_COMPARE_ARRAYS( y.values(), x.values() );
  TEST_EQUALITY(computeSum<Scalar>(x), sumIntegers<Scalar>(n));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(1, 2, 3, 4, 5, 6) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  equalRanges_fallthrough )


//
//  1) Equal matching ranges:
//     subvecs  [......]
//     rtop     [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  equalRanges, Scalar )
{
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;

  TOpAssignVectors<Scalar> vecAssignOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(vecAssignOp)), 0, n);

  ConstSubVectorView<Scalar> x = createFilledSubVectorView<Scalar>(n);
  TEST_EQUALITY(computeSum<Scalar>(x), sumIntegers<Scalar>(n));
  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  op.apply_op(tuple(x), tuple(y), null);
  TEST_COMPARE_ARRAYS( y.values(), x.values() );
  TEST_EQUALITY(computeSum<Scalar>(x), sumIntegers<Scalar>(n));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(1, 2, 3, 4, 5, 6) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  equalRanges )


//
//  2) Match first part of range (inner subvec):
//     subvecs  [..]
//     rtop     [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  matchFirst_inner_subvec, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int m = 2;

  TOpAssignVectors<Scalar> vecAssignOp;
  RTOpSubRangeDecorator<Scalar>
    op(RCP<RTOpT<Scalar> >(rcpFromRef(vecAssignOp)), 0, m);

  ConstSubVectorView<Scalar> x = createFilledSubVectorView<Scalar>(m);
  TEST_EQUALITY( computeSum<Scalar>(x), sumIntegers<Scalar>(m));
  SubVectorView<Scalar> y = createSubVectorView<Scalar>(m);
  op.apply_op(tuple(x), tuple(y), null);
  TEST_COMPARE_ARRAYS( y.values(), x.values() );
  TEST_EQUALITY( computeSum<Scalar>(x), sumIntegers<Scalar>(m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(1, 2) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  matchFirst_inner_subvec )


//
//  3) Match middle part of range (inner subvec):
//     subvecs    [..]
//     rtop     [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  matchMiddle_inner_subvec, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)));

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(m);
  y.setGlobalOffset(m);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(2*m) - sumIntegers<Scalar>(m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(3, 4) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  matchMiddle_inner_subvec )


//
//  4) Match last part of range (inner subvec):
//     subvecs      [..]
//     rtop     [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  matchLast_inner_subvec, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)));

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(m);
  y.setGlobalOffset(2*m);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(3*m) - sumIntegers<Scalar>(2*m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(5, 6) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  matchLast_inner_subvec )


//
//  5) Match first part of range (inner rtop):
//     subvecs  [......]
//     rtop     [..]

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  matchFirst_inner_rtop, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;
  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)), 0, m);

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(1, 2, 0, 0, 0, 0) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  matchFirst_inner_rtop )


//
//  6) Match middle part of range (inner rtop):
//     subvecs  [......]
//     rtop       [..]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  matchMiddle_inner_rtop, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;
  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)), m, m);

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(2*m) - sumIntegers<Scalar>(m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(0, 0, 3, 4, 0, 0) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  matchMiddle_inner_rtop )


//
//  7) Match last part of range (inner rtop):
//     subvecs  [......]
//     rtop         [..]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  matchLast_inner_rtop, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;
  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)), 2*m, m);

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(3*m) - sumIntegers<Scalar>(2*m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(0, 0, 0, 0, 5, 6) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  matchLast_inner_rtop )


//
//  8) Non-overlapping subvec first
//     subvecs  [......]
//     rtop         [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  nonoverlap_subvec_first, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;
  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)), 2*m, n);

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(3*m) - sumIntegers<Scalar>(2*m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(0, 0, 0, 0, 5, 6) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  nonoverlap_subvec_first )


//
//  9) Non-overlapping rtop first
//     subvecs      [......]
//     rtop     [......]
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RTOpSubRangeDecorator,
  nonoverlap_rtop_first, Scalar )
{
  using Teuchos::Range1D;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::as;

  const int n = 6;
  const int m = 2;

  TOpSetAssendingValues<Scalar> setAssigingOp;
  RTOpSubRangeDecorator<Scalar> op(RCP<RTOpT<Scalar> >(rcpFromRef(setAssigingOp)), 0, n);

  SubVectorView<Scalar> y = createSubVectorView<Scalar>(n);
  y.setGlobalOffset(2*m);
  op.apply_op(null, tuple(y), null);
  TEST_EQUALITY(computeSum<Scalar>(y), sumIntegers<Scalar>(3*m) - sumIntegers<Scalar>(2*m));
  TEST_COMPARE_ARRAYS( y.values(), tuple<Scalar>(5, 6, 0, 0, 0, 0) );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( RTOpSubRangeDecorator,
  nonoverlap_rtop_first )


} // namespace RTOpPack
