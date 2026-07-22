// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_NORM2_HPP
#define RTOPPACK_ROP_NORM2_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpNorm2EleWiseReduction
{
public:
  void operator()( const Scalar &v0, Scalar &reduct ) const
    {
      reduct += ScalarTraits<Scalar>::conjugate(v0)*v0;
    }
};


/** \brief Two (Euclidean) norm reduction operator: <tt>result = sqrt( sum(
 * conj(v0[i])*v0[i], i=0...n-1 ) )</tt>.
 */
template<class Scalar>
class ROpNorm2
  : public ROp_1_ScalarReduction<Scalar, Scalar, ROpNorm2EleWiseReduction<Scalar> >
{
public:
  /** \brief . */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief . */
  ROpNorm2()
    {
      this->setOpNameBase("ROpNorm2");
    }
  /** \brief . */
  typename ST::magnitudeType operator()(const ReductTarget& reduct_obj) const
    {
      const Scalar sqrt_reduct = ST::squareroot(this->getRawVal(reduct_obj));
      return ST::magnitude(sqrt_reduct);
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_NORM2_HPP
