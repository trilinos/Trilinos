// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_SCALE_VECTOR_HPP
#define RTOPPACK_TOP_SCALE_VECTOR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpScaleVector. */
template<class Scalar>
class TOpScaleVectorEleWiseTransformation
{
public:
  TOpScaleVectorEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( Scalar &z0 ) const
    {
      z0 *= alpha_;
    }
private:
  Scalar alpha_;
};


/** \brief Simple transformation operator that scales every vector element by
 * a scalar <tt>alpha</tt>.
 */
template<class Scalar>
class TOpScaleVector
  : public TOp_0_1_Base<Scalar, TOpScaleVectorEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_0_1_Base<Scalar, TOpScaleVectorEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpScaleVector( const Scalar &alpha )
    : base_t(TOpScaleVectorEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpScaleVector");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_SCALE_VECTOR_HPP
