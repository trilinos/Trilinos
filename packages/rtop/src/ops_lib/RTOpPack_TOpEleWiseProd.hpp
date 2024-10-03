// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ELE_WISE_PROD_HPP
#define RTOPPACK_TOP_ELE_WISE_PROD_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpEleWiseProd. */
template<class Scalar>
class TOpEleWiseProdEleWiseTransformation
{
public:
  TOpEleWiseProdEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( const Scalar &v0, const Scalar &v1, Scalar &z0 ) const
    {
      z0 += alpha_ * v0 * v1;
    }
private:
  Scalar alpha_;
};


/** \brief Element-wise product transformation operator: <tt>z0[i] +=
 * alpha*v0[i]*v1[i], i=0...n-1</tt>.
 */
template<class Scalar>
class TOpEleWiseProd
  : public TOp_2_1_Base<Scalar, TOpEleWiseProdEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_2_1_Base<Scalar, TOpEleWiseProdEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpEleWiseProd( const Scalar &alpha )
    : base_t(TOpEleWiseProdEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpEleWiseProd");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ELE_WISE_PROD_HPP
