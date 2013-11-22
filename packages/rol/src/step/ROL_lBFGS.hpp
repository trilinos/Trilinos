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

#ifndef ROL_LBFGS_H
#define ROL_LBFGS_H

/** \class ROL::lBFGS
    \brief Provides defintions for limited-memory BFGS operators.
*/

namespace ROL {

template<class Real>
class lBFGS : public Secant<Real> {
public:
  lBFGS(int M) : Secant<Real>(M) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x ) {
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Hv.set(v);
    std::vector<Real> alpha(state->current+1,0.0);
    for (int i = state->current; i>=0; i--) {
      alpha[i]  = state->iterDiff[i]->dot(Hv);
      alpha[i] /= state->product[i];
      Hv.axpy(-alpha[i],*(state->gradDiff[i]));
    }

    // Apply initial inverse Hessian approximation to v   
    Teuchos::RCP<Vector<Real> > tmp = Hv.clone();
    this->applyH0(*tmp,Hv,x);
    Hv.set(*tmp);

    Real beta = 0.0;
    for (int i = 0; i <= state->current; i++) {
      beta  = state->gradDiff[i]->dot(Hv);
      beta /= state->product[i];
      Hv.axpy((alpha[i]-beta),*(state->iterDiff[i]));
    }
  }

  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x ) { 
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    // Apply initial Hessian approximation to v   
    this->applyB0(Bv,v,x);

    std::vector<Teuchos::RCP<Vector<Real> > > a(state->current+1);
    std::vector<Teuchos::RCP<Vector<Real> > > b(state->current+1);
    Real bv = 0.0, av = 0.0, bs = 0.0, as = 0.0;
    for (int i = 0; i <= state->current; i++) {
      b[i] = v.clone();
      b[i]->set(*(state->gradDiff[i]));
      b[i]->scale(1.0/sqrt(state->product[i]));
      bv = b[i]->dot(v);
      Bv.axpy(bv,*b[i]);

      a[i] = v.clone();
      this->applyB0(*a[i],*(state->iterDiff[i]),x);

      for (int j = 0; j < i; j++) {
        bs = b[j]->dot(*(state->iterDiff[i]));
        a[i]->axpy(bs,*b[j]);
        as = a[j]->dot(*(state->iterDiff[i]));
        a[i]->axpy(-as,*a[j]);
      }
      as = a[i]->dot(*(state->iterDiff[i]));
      a[i]->scale(1.0/sqrt(as));
      av = a[i]->dot(v);
      Bv.axpy(-av,*a[i]);
    }
  }

};

}

#endif
