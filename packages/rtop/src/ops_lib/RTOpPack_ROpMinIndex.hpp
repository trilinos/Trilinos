// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_MIN_INDEX_HPP
#define RTOPPACK_ROP_MIN_INDEX_HPP


#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpMinIndexEleWiseReductionOp {
public:
  /** \brief . */
  void operator()(const index_type i, const Scalar &v0,
    ScalarIndex<Scalar> &reduct) const
    {
      if(
        v0 < reduct.scalar
        ||
        ( v0 == reduct.scalar && i < reduct.index )
        )
      {
        reduct = ScalarIndex<Scalar>(v0, i);
      }
    }
};


/** \brief. */
template<class Scalar>
class ROpMinIndexReductObjReductionOp {
public:
  /** \brief . */
  void operator()(
    const ScalarIndex<Scalar>& in_reduct, ScalarIndex<Scalar>& inout_reduct
    ) const
    {
      if(
        in_reduct.scalar < inout_reduct.scalar
        ||
        (
          in_reduct.scalar == inout_reduct.scalar
          &&
          in_reduct.index < inout_reduct.index
          )
        )
      {
        inout_reduct = in_reduct;
      }
    }
};


/** \brief Returns the minimum element and its index:
 * <tt>result.scalar = x(k)</tt> and <tt>result.index = k</tt> such
 * that <tt>x(k) <= x(i)</tt> for <tt>i=0...n-1</tt> and <tt>k</tt> is
 * the minimum index to break ties.
 */
template<class Scalar>
class ROpMinIndex
  : public ROp_1_CoordVariantScalarReduction<
      Scalar,
      ScalarIndex<Scalar>,
      ROpMinIndexEleWiseReductionOp<Scalar>,
      ROpMinIndexReductObjReductionOp<Scalar> >
{
public:
  /** \brief . */
  ROpMinIndex()
    {
      this->setOpNameBase("ROpMinIndex");
      this->initReductObjValue(
        ScalarIndex<Scalar>(+ScalarTraits<Scalar>::rmax(), -1));
    }
  /** \brief . */
  ScalarIndex<Scalar> operator()(const ReductTarget& reduct_obj) const
    { return this->getRawVal(reduct_obj); }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_MIN_INDEX_HPP
