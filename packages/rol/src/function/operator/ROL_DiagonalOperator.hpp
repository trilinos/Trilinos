// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DIAGONALOPERATOR_H
#define ROL_DIAGONALOPERATOR_H

#include "ROL_Vector.hpp"
#include "ROL_Elementwise_Function.hpp"
#include "ROL_LinearOperator.hpp"


/** @ingroup func_group
    \class ROL::DiagonalOperator
    \brief Provides the interface to apply a diagonal operator which acts like
           elementwise multiplication when apply() is used and elementwise division
           when applyInverse() is used.
*/

namespace ROL {

template<class Real>
class DiagonalOperator : public LinearOperator<Real> {

private:
 
  ROL::Ptr<Vector<Real> >             diag_;

  const Elementwise::Multiply<Real>       mult_;
  const Elementwise::Divide<Real>         div_;   

public:

  DiagonalOperator( const Vector<Real> &diag ) : diag_(diag.clone()) {
    diag_->set(diag);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    diag_->set(x);
  }

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v); 
    Hv.applyBinary( mult_, *diag_ );
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v); 
    Hv.applyBinary( div_, *diag_ );
  }

};

} // namespace ROL




#endif // ROL_DIAGONALOPERATOR_H

