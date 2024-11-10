// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_SET_ELEMENT_HPP
#define RTOPPACK_TOP_SET_ELEMENT_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Teuchos_as.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation for TOpSetElement. */
template<class Scalar>
class TOpSetElementEleWiseTransformation
{
public:
  /** \brief . */
  TOpSetElementEleWiseTransformation( const Ordinal &global_i_in = -1,
    const Scalar &val_i_in = static_cast<Scalar>(0.0) )
    :global_i_(global_i_in), val_i_(val_i_in)
    {}
  /** \brief . */
  Ordinal global_i() const
    {
      return global_i_;
    }
  /** \brief . */
  void operator()( const Ordinal global_i_in, Scalar &z0 ) const
    {
      if (global_i_in == global_i_) {
        z0 = val_i_;
      }
    }
private:
  Ordinal global_i_;
  Scalar val_i_;
};


/** \brief Set the elements of a vector to: <tt>z0[i] = i+global_i+1, i=0...n-1</tt>.
 */
template<class Scalar>
class TOpSetElement
  : public TOp_0_1_CoordVariantBase<Scalar, TOpSetElementEleWiseTransformation<Scalar> >
{
public:
  /** \brief . */
  TOpSetElement(const Ordinal &global_i_in = -1,
    const Scalar &val_i_in = static_cast<Scalar>(0.0))
    {
      this->setOpNameBase("TOpSetElement");
      this->setEleWiseTransformation(
        TOpSetElementEleWiseTransformation<Scalar>(global_i_in, val_i_in));
    }
  /** \brief . */
  void initialize(const Ordinal &global_i_in, const Scalar &val_i_in)
    { 
      this->setEleWiseTransformation(
        TOpSetElementEleWiseTransformation<Scalar>(global_i_in, val_i_in));
    }
protected:
  /** \brief . */
  virtual Range1D range_impl() const
    {
      const Ordinal i = this->getEleWiseTransformation().global_i();
      return Range1D(i, i);
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_SET_ELEMENT_HPP
