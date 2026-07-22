// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_WEIGHTED_NORM2_HPP
#define RTOPPACK_ROP_WEIGHTED_NORM2_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpWeightedNorm2EleWiseReduction
{
public:
  void operator()( const Scalar &v0, const Scalar &v1, Scalar &reduct ) const
    {
      reduct += v0 * ScalarTraits<Scalar>::conjugate(v1)*v1;
    }
};


/** \brief Weighted Two (Euclidean) norm reduction operator: <tt>result =
 * sqrt( sum( v0[i]*conj(v1[i])*v1[i], i=0...n-1 ) )</tt>.
 */
template<class Scalar>
class ROpWeightedNorm2
  : public ROp_2_ScalarReduction<Scalar, Scalar,
      ROpWeightedNorm2EleWiseReduction<Scalar> >
{
public:
  /** \brief . */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief . */
  ROpWeightedNorm2()
    {
      this->setOpNameBase("ROpWeightedNorm2");
    }
  /** \brief . */
  typename ST::magnitudeType operator()(const ReductTarget& reduct_obj ) const
    { return ST::magnitude(ST::squareroot(this->getRawVal(reduct_obj))); }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_WEIGHTED_NORM2_HPP
