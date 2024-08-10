// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SECANTFACTORY_H
#include "ROL_SecantFactory.hpp"
#else

#ifndef ROL_BARZILAIBORWEIN_H
#define ROL_BARZILAIBORWEIN_H

/** \class ROL::BarzilaiBorwein
    \brief Provides definitions for Barzilai-Borwein operators.
*/

#include "ROL_Secant.hpp"

namespace ROL {

template<class Real>
class BarzilaiBorwein : public Secant<Real> {
private:

  int type_;

public:
  BarzilaiBorwein(int type = 1) : Secant<Real>(1), type_(type) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const {
    // Get Generic Secant State
    const ROL::Ptr<SecantState<Real> >& state = Secant<Real>::get_state();

    Hv.set(v.dual());
    if ( state->iter != 0 && state->current != -1 ) {
      if ( type_ == 1 ) {
        Real yy = state->gradDiff[state->current]->dot(*(state->gradDiff[state->current]));
        Hv.scale(state->product[state->current]/yy);
      }
      else if ( type_ == 2 ) {
        Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
        Hv.scale(ss/state->product[state->current]);
      }
    }
  }

  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const {
    // Get Generic Secant State
    const ROL::Ptr<SecantState<Real> >& state = Secant<Real>::get_state();

    Bv.set(v.dual());
    if ( state->iter != 0 && state->current != -1 ) {
      if ( type_ == 1 ) {
        Real yy = state->gradDiff[state->current]->dot(*(state->gradDiff[state->current]));
        Bv.scale(yy/state->product[state->current]);
      }
      else if ( type_ == 2 ) {
        Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
        Bv.scale(state->product[state->current]/ss);
      }
    }
  }
};

}

#endif
#endif
