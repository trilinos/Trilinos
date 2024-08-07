// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_IDENTITYOPERATOR_H
#define ROL_IDENTITYOPERATOR_H

#include "ROL_LinearOperator.hpp"


/** @ingroup func_group
    \class ROL::IdentityOperator
    \brief Multiplication by unity
*/

namespace ROL {

template<class Real>
class IdentityOperator : public LinearOperator<Real> {

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v); 
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v); 
  }

};

} // namespace ROL




#endif // ROL_NULLOPERATOR_H

