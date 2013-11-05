//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_BARZILAIBORWEIN_H
#define ROL_BARZILAIBORWEIN_H

/** \class ROL::BarzilaiBorwein
    \brief Provides defintions for Barzilai-Borwein operators.
*/

namespace ROL {

template<class Real>
class BarzilaiBorwein : public Secant<Real> {
public:
  BarzilaiBorwein(void) : Secant<Real>(1) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, const int iter ) {
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Hv.set(v);
    if ( iter != 0 && state->current != -1 ) {
      Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
      Hv.scale(state->product[state->current]/ss);
    }
  }

  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x, const int iter ) { 
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Bv.set(v);
    if ( iter != 0 && state->current != -1 ) {
      Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
      Bv.scale(ss/state->product[state->current]);
    }
  }

};

}

#endif
