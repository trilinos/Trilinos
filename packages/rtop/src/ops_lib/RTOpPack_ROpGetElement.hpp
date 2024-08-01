// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_GET_ELEMENT_HPP
#define RTOPPACK_ROP_GET_ELEMENT_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpGetElementEleWiseReductionOp {
public:
  /** \brief . */
  ROpGetElementEleWiseReductionOp(const Ordinal &global_i_in = -1)
    :global_i_(global_i_in)
    {}
  /** \brief . */
  Ordinal global_i() const
    {
      return global_i_;
    }
  /** \brief . */
  void operator()(const index_type i, const Scalar &v0, Scalar &reduct) const
    {
      if (i == global_i_) {
        reduct = v0;
      }
    }
private:
  Ordinal global_i_;
};


/** \brief Returns the value of the element at index <tt>global_i</tt>.
 *
 * Warning! If the element is not found, then 0 is returned!
 */
template<class Scalar>
class ROpGetElement
  : public ROp_1_CoordVariantScalarReduction<
        Scalar,
        Scalar,
        ROpGetElementEleWiseReductionOp<Scalar>
      >
{
public:
  /** \brief . */
  ROpGetElement(const Ordinal &global_i)
    {
      this->setOpNameBase("ROpGetElement");
      this->initialize(global_i);
      this->initReductObjValue(ScalarTraits<Scalar>::zero());
    }
  /** \brief . */
  void initialize(const Ordinal &global_i)
    { 
      this->setEleWiseReduction(ROpGetElementEleWiseReductionOp<Scalar>(global_i));
    }
  /** \brief . */
  Scalar operator()(const ReductTarget& reduct_obj ) const
    {
      return this->getRawVal(reduct_obj);
    }
protected:
  /** \brief . */
  virtual Range1D range_impl() const
    {
      const Ordinal i = this->getEleWiseReduction().global_i();
      return Range1D(i, i);
    }

};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_GET_ELEMENT_HPP
