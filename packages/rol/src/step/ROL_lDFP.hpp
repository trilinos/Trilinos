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

#ifndef ROL_LDFP_H
#define ROL_LDFP_H

/** \class ROL::lDFP
    \brief Provides defintions for limited-memory DFP operators.
*/

namespace ROL {

template<class Real>
class lDFP : public Secant<Real> {
public:
  lDFP(int M) : Secant<Real>(M) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, const int iter ) { 
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    // Apply initial Hessian approximation to v   
    applyH0(Hv,v,x,iter);

    std::vector<Teuchos::RCP<Vector<Real> > > a(state->current+1);
    std::vector<Teuchos::RCP<Vector<Real> > > b(state->current+1);
    Real bv = 0.0, av = 0.0, bs = 0.0, as = 0.0;
    for (int i = 0; i <= state->current; i++) {
      b[i] = v.clone();
      b[i]->set(*(state->iterDiff[i]));
      b[i]->scale(1.0/sqrt(state->product[i]));
      bv = b[i]->dot(v);
      Hv.axpy(bv,*b[i]);

      a[i] = v.clone();
      applyH0(*a[i],*(state->gradDiff[i]),x,iter);

      for (int j = 0; j < i; j++) {
        bs = b[j]->dot(*(state->gradDiff[i]));
        a[i]->axpy(bs,*b[j]);
        as = a[j]->dot(*(state->gradDiff[i]));
        a[i]->axpy(-as,*a[j]);
      }
      as = a[i]->dot(*(state->gradDiff[i]));
      a[i]->scale(1.0/sqrt(as));
      av = a[i]->dot(v);
      Hv.axpy(-av,*a[i]);
    }
  }

  // Apply Initial Secant Approximate Hessian
  void applyH0( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, const int iter ) {
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Hv.set(v);
    if (iter != 0 && state->current != -1) {
      Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
      Hv.scale(state->product[state->current]/ss);
    }
  }


  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x, const int iter ) {
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Bv.set(v);
    std::vector<Real> alpha(state->current+1,0.0);
    for (int i = state->current; i>=0; i--) {
      alpha[i]  = state->gradDiff[i]->dot(Bv);
      alpha[i] /= state->product[i];
      Bv.axpy(-alpha[i],*(state->iterDiff[i]));
    }

    // Apply initial inverse Hessian approximation to v   
    Teuchos::RCP<Vector<Real> > tmp = Bv.clone();
    applyB0(*tmp,Bv,x,iter);
    Bv.set(*tmp);

    Real beta = 0.0;
    for (int i = 0; i <= state->current; i++) {
      beta  = state->iterDiff[i]->dot(Bv);
      beta /= state->product[i];
      Bv.axpy((alpha[i]-beta),*(state->gradDiff[i]));
    }
  }

  // Apply Initial Secant Approximate Hessian 
  void applyB0( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x, const int iter ) {
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Bv.set(v);
    if (iter != 0 && state->current != -1) {
      Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
      Bv.scale(ss/state->product[state->current]);
    }
  }

};

}

#endif
