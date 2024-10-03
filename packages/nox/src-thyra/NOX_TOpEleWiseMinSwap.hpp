// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ELE_WISE_MIN_SWAP_HPP
#define RTOPPACK_TOP_ELE_WISE_MIN_SWAP_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Thyra_VectorBase.hpp"
#include <cmath>

namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpEleWiseMinSwap. */
template<class Scalar>
class TOpEleWiseMinSwapEleWiseTransformation
{
public:
  void operator()( const Scalar &v0, Scalar &z0 ) const
    {
      if (std::abs(z0) > v0)
    z0 = v0 * (z0/std::abs(z0));
    }
};


/** \brief Element-wise product update transformation operator:
 * <tt>z0[i] *= min(v0[i],abs(z0[i]) * z0[i]/abs(z0[i]), i=0...n-1</tt>.
 */
template<class Scalar>
class TOpEleWiseMinSwap
  : public TOp_1_1_Base<Scalar, TOpEleWiseMinSwapEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_1_1_Base<Scalar, TOpEleWiseMinSwapEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpEleWiseMinSwap()
    : base_t(TOpEleWiseMinSwapEleWiseTransformation<Scalar>())
    {
      this->setOpNameBase("TOpEleWiseMinSwap");
    }
};


} // namespace RTOpPack


namespace Thyra {

/** \brief Element-wise min swap:
 * <tt>y(i) *= min(x(i),abs(y(i)) * y(i)/abs(y(i), i = 0...y->space()->dim()-1</tt>.
 *
 * \relates VectorBase
 */
template<class Scalar>
void ele_wise_min_swap( const ::Thyra::VectorBase<Scalar>& x,
            const Teuchos::Ptr< ::Thyra::VectorBase<Scalar> > &y )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpEleWiseMinSwap<Scalar> ele_wise_min_swap_op;
  ::Thyra::applyOp<Scalar>( ele_wise_min_swap_op, tuple(ptrInArg(x)),
    tuple(y), null );
}

}

#endif // RTOPPACK_TOP_ELE_WISE_MIN_SWAP_HPP
