// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_SET_ASSENDING_VALUES_HPP
#define RTOPPACK_TOP_SET_ASSENDING_VALUES_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Teuchos_as.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation for TOpSetAssendingValues. */
template<class Scalar>
class TOpSetAssendingValuesEleWiseTransformation
{
public:
  TOpSetAssendingValuesEleWiseTransformation( const Scalar &offset )
    {
      offset_ = offset;
    }
  void operator()( const Ordinal global_i, Scalar &z0 ) const
    {
      z0 = offset_ + Teuchos::as<Scalar>(global_i + 1);
    }
private:
  Scalar offset_;
  TOpSetAssendingValuesEleWiseTransformation(); // Not defined
};


/** \brief Set the elements of a vector to: <tt>z0[i] = i+offset+1, i=0...n-1</tt>.
 */
template<class Scalar>
class TOpSetAssendingValues
  : public TOp_0_1_CoordVariantBase<Scalar, TOpSetAssendingValuesEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_0_1_CoordVariantBase<Scalar, TOpSetAssendingValuesEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpSetAssendingValues(const Scalar &offset = static_cast<Scalar>(0.0) )
    : base_t(TOpSetAssendingValuesEleWiseTransformation<Scalar>(offset))
    {
      this->setOpNameBase("TOpSetAssendingValues");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_SET_ASSENDING_VALUES_HPP
