// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ELE_WISE_PROD_UPDATE_HPP
#define RTOPPACK_TOP_ELE_WISE_PROD_UPDATE_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpEleWiseProdUpdate. */
template<class Scalar>
class TOpEleWiseProdUpdateEleWiseTransformation
{
public:
  TOpEleWiseProdUpdateEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( const Scalar &v0, Scalar &z0 ) const
    {
      z0 *= alpha_ * v0;
    }
private:
  Scalar alpha_;
};


/** \brief Element-wise product update transformation operator: 
 * <tt>z0[i] *= alpha*v0[i], i=0...n-1</tt>.
 */
template<class Scalar>
class TOpEleWiseProdUpdate
  : public TOp_1_1_Base<Scalar, TOpEleWiseProdUpdateEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_1_1_Base<Scalar, TOpEleWiseProdUpdateEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpEleWiseProdUpdate( const Scalar &alpha )
    : base_t(TOpEleWiseProdUpdateEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpEleWiseProdUpdate");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ELE_WISE_PROD_UPDATE_HPP
