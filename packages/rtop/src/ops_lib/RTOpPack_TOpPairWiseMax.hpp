// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_PAIR_WISE_MAX_HPP
#define RTOPPACK_TOP_PAIR_WISE_MAX_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Pair-wise transformation operator for TOpPairWiseMax. */
template<class Scalar>
class TOpPairWiseMaxPairWiseTransformation
{
public:
  TOpPairWiseMaxPairWiseTransformation(const Scalar &alpha)
  : alpha_(alpha)
  {}
  void operator()( const Scalar &v0, const Scalar &v1, Scalar &z0 ) const
    {
      z0 = alpha_ * std::max(v0, v1);
    }
private:
  Scalar alpha_;
};


/** \brief Pair-wise Maximum transformation operator:
 * <tt>z0[i] = alpha*max(v0[i],v1[i]), i=0...n-1</tt>.
 */
template<class Scalar>
class TOpPairWiseMax
  : public TOp_2_1_Base<Scalar, TOpPairWiseMaxPairWiseTransformation<Scalar> >
{
public:
  typedef TOp_2_1_Base<Scalar, TOpPairWiseMaxPairWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpPairWiseMax(const Scalar &alpha )
    : base_t(TOpPairWiseMaxPairWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpPAIRWiseMax");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_PAIR_WISE_MAX_HPP
