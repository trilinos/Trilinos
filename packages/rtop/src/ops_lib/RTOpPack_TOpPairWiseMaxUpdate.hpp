// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_PAIR_WISE_MAX_UPDATE_HPP
#define RTOPPACK_TOP_PAIR_WISE_MAX_UPDATE_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Pair-wise transformation operator for TOpPairWiseMaxUpdate. */
template<class Scalar>
class TOpPairWiseMaxUpdatePairWiseTransformation
{
public:
  TOpPairWiseMaxUpdatePairWiseTransformation(const Scalar &alpha )
  : alpha_(alpha)
  {}
  void operator()( const Scalar &v0, Scalar &z0 ) const
    {
      z0 = alpha_ * std::max(z0, v0);
    }
private:
  Scalar alpha_;
};


/** \brief Pair-wise Maximum update transformation operator: 
 * <tt>z0[i] = alpha*max(z0[i],v0[i]), i=0...n-1</tt>.
 */
template<class Scalar>
class TOpPairWiseMaxUpdate
  : public TOp_1_1_Base<Scalar, TOpPairWiseMaxUpdatePairWiseTransformation<Scalar> >
{
public:
  typedef TOp_1_1_Base<Scalar, TOpPairWiseMaxUpdatePairWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpPairWiseMaxUpdate( const Scalar &alpha )
    : base_t(TOpPairWiseMaxUpdatePairWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpPairWiseMaxUpdate");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_PAIR_WISE_MAX_UPDATE_HPP
