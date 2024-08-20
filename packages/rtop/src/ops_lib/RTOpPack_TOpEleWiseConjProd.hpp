// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ELE_WISE_CONJ_PROD_HPP
#define RTOPPACK_TOP_ELE_WISE_CONJ_PROD_HPP


#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpEleWiseConjProd. */
template<class Scalar>
class TOpEleWiseConjProdEleWiseTransformation
{
public:
  TOpEleWiseConjProdEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( const Scalar &v0, const Scalar &v1, Scalar &z0 ) const
    {
      z0 += alpha_ * ScalarTraits<Scalar>::conjugate(v0) * v1;
    }
private:
  Scalar alpha_;
};


/** \brief Element-wise product transformation operator: <tt>z0[i] +=
 * alpha*conj(v0[i])*v1[i], i=0...n-1</tt>.
 */
template<class Scalar>
class TOpEleWiseConjProd
  : public TOp_2_1_Base<Scalar, TOpEleWiseConjProdEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_2_1_Base<Scalar, TOpEleWiseConjProdEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpEleWiseConjProd( const Scalar &alpha )
    : base_t(TOpEleWiseConjProdEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpEleWiseConjProd");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ELE_WISE_CONJ_PROD_HPP
