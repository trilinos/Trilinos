// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_AXPY_HPP
#define RTOPPACK_TOP_AXPY_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpAXPY. */
template<class Scalar>
class TOpAXPYEleWiseTransformation
{
public:
  TOpAXPYEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( const Scalar &v0, Scalar &z0 ) const
    {
      z0 += alpha_ * v0;
    }
private:
  Scalar alpha_;
};


/** \brief AXPY transformation operator: <tt>z0[i] += alpha*v0[i],
 * i=0...n-1</tt>.
 */
template<class Scalar>
class TOpAXPY
  : public TOp_1_1_Base<Scalar, TOpAXPYEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_1_1_Base<Scalar, TOpAXPYEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpAXPY( const Scalar &alpha )
    : base_t(TOpAXPYEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpAXPY");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_AXPY_HPP
