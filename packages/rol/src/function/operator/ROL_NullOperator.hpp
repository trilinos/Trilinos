// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NULLOPERATOR_H
#define ROL_NULLOPERATOR_H

#include "ROL_LinearOperator.hpp"


/** @ingroup func_group
    \class ROL::NullOperator
    \brief Multiplication by zero
*/

namespace ROL {

template<class Real>
class NullOperator : public LinearOperator<Real> {

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.zero(); 
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
                                ">>> ERROR (ROL_NullOperator, applyInverse): "
                                "Null Operator has no inverse.");
  }

};

} // namespace ROL




#endif // ROL_NULLOPERATOR_H

