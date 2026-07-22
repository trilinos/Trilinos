// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RANDOMVECTOR_H
#define ROL_RANDOMVECTOR_H

#include "ROL_Vector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Elementwise_Function.hpp"


namespace ROL {

/** @ingroup la_group
    \function RandomizeVector
    \brief Fill a ROL::Vector with uniformly-distributed random numbers
           in the interval [lower,upper]

    \f[ x_i \in \mathcal{U}(l,u)\;\forall\l i \f]

*/

template<class Real> 
void RandomizeVector( Vector<Real> &x, const Real &lower=0.0, const Real &upper=1.0 ) {

  Elementwise::UniformlyRandom<Real> ur(lower,upper);
  x.applyUnary(ur);
}

/** @ingroup la_group
    \function RandomizeFeasibleVector
    \brief Fill a ROL::Vector with uniformly-distributed random numbers
           which satisfy the supplied bound constraint

    \f[ x_i \in \mathcal{U}(l_i,u_i)\;\forall\l i \f]
*/

template<class Real>
void RandomizeFeasibleVector( Vector<Real> &x, BoundConstraint<Real> &bnd ) {
  const ROL::Ptr<const Vector<Real> > u = bnd.getUpperBound();
  const ROL::Ptr<const Vector<Real> > l = bnd.getLowerBound();

  Elementwise::UniformlyRandomMultiply<Real> urm;
  
  x.set(*u);
  x.axpy(-1.0,*l);
  x.applyUnary(urm);
  x.plus(*l);
}


} // namespace ROL

#endif // ROL_RANDOMVECTOR_H
