// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_MIN_INDEX_GREATER_THAN_BOUND_HPP
#define RTOPPACK_ROP_MIN_INDEX_GREATER_THAN_BOUND_HPP

#include "RTOpPack_ROpMinIndex.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpMinIndexGreaterThanBoundEleWiseReductionOp {
public:
  /** \brief . */
  ROpMinIndexGreaterThanBoundEleWiseReductionOp(
    const Scalar &bound = ScalarTraits<Scalar>::zero()
    )
    :bound_(bound)
    {}
  /** \brief . */
  void operator()(const index_type i, const Scalar &v0,
    ScalarIndex<Scalar> &reduct) const
    {
      if(
        v0 >  bound_
        &&
        (
          v0 < reduct.scalar
          ||
          ( v0 == reduct.scalar && i < reduct.index )
          )
        )
      {
        reduct = ScalarIndex<Scalar>(v0, i);
      }
    }
private:
  Scalar bound_;
};


/** \brief Returns the minimum element greater than some bound along
 * with its index: <tt>result.scalar = x(k)</tt> and <tt>result.index
 * = k</tt> such that <tt>x(k) <= x(i)</tt> for all <tt>i</tt> where
 * <tt>x(i) > bound</tt> and <tt>k</tt> is the minimum index to break
 * ties.
 *
 * If no element is greater than <tt>bound</tt> then <tt>results.index
 * < 0</tt>.
 */
template<class Scalar>
class ROpMinIndexGreaterThanBound
  : public ROp_1_CoordVariantScalarReduction<
      Scalar,
      ScalarIndex<Scalar>,
      ROpMinIndexGreaterThanBoundEleWiseReductionOp<Scalar>,
      ROpMinIndexReductObjReductionOp<Scalar> >
{
public:

  /** \brief . */
  ROpMinIndexGreaterThanBound(
    const Scalar &bound_in = Teuchos::ScalarTraits<Scalar>::zero()
    )
    {
      this->setOpNameBase("ROpMinIndexGreaterThanBound");
      bound(bound_in);
      this->initReductObjValue(
        ScalarIndex<Scalar>(+ScalarTraits<Scalar>::rmax(), -1));
    }
  /** \brief . */
  void bound(const Scalar& bound_in)
    { 
      this->setEleWiseReduction(
        ROpMinIndexGreaterThanBoundEleWiseReductionOp<Scalar>(bound_in)
        );
    }
  /** \brief . */
  ScalarIndex<Scalar> operator()(const ReductTarget& reduct_obj ) const
    { return this->getRawVal(reduct_obj); }
  
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_MIN_INDEX_GREATER_THAN_BOUND_HPP
