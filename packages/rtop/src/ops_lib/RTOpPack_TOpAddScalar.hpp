// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ADD_SCALAR_HPP
#define RTOPPACK_TOP_ADD_SCALAR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpAddScalar. */
template<class Scalar>
class TOpAddScalarEleWiseTransformation
{
public:
  TOpAddScalarEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( Scalar &z0 ) const
    {
      z0 += alpha_;
    }
private:
  Scalar alpha_;
};


/** \brief Add a scalar to a vector transformation operator: <tt>z0[i] +=
 * alpha, i=0...n-1</tt>.
 */
template<class Scalar>
class TOpAddScalar
  : public TOp_0_1_Base<Scalar, TOpAddScalarEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_0_1_Base<Scalar, TOpAddScalarEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpAddScalar( const Scalar &alpha )
    : base_t(TOpAddScalarEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpAddScalar");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ADD_SCALAR_HPP
