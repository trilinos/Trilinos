// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ASSIGN_SCALAR_HPP
#define RTOPPACK_TOP_ASSIGN_SCALAR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpAssignScalar. */
template<class Scalar>
class TOpAssignScalarEleWiseTransformation
{
public:
  TOpAssignScalarEleWiseTransformation( const Scalar &alpha )
    {
      alpha_ = alpha;
    }
  void operator()( Scalar &z0 ) const
    {
      z0 = alpha_;
    }
private:
  Scalar alpha_;
  TOpAssignScalarEleWiseTransformation(); // Not defined
};


/** \brief Assign a scalar to a vector transformation operator: <tt>z0[i] =
 * alpha, i=0...n-1</tt>.
 */
template<class Scalar>
class TOpAssignScalar
  : public TOp_0_1_Base<Scalar, TOpAssignScalarEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_0_1_Base<Scalar, TOpAssignScalarEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpAssignScalar( const Scalar &alpha )
    : base_t(TOpAssignScalarEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpAssignScalar");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ASSIGN_SCALAR_HPP
