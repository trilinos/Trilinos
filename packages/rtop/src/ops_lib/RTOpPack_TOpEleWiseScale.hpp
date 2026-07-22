// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ELE_WISE_SCALE_HPP
#define RTOPPACK_TOP_ELE_WISE_SCALE_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise vector scaling op for TOpEleWiseScaling. */
template<class Scalar>
class TOpEleWiseScaleEleWiseTransformation
{
public:
  void operator()( const Scalar &v0, Scalar &z0 ) const
    {
      z0 *= v0;
    }
private:
  Scalar alpha_;
};


/** \brief Element-wise vector scaling: <tt>z0[i] *= v0[i], i=0...n-1</tt>.
 */
template<class Scalar>
class TOpEleWiseScale
  : public TOp_1_1_Base<Scalar, TOpEleWiseScaleEleWiseTransformation<Scalar> >
{
public:
  /** \brief . */
  TOpEleWiseScale()
    {
      this->setOpNameBase("TOpEleWiseScale");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ELE_WISE_SCALE_HPP
