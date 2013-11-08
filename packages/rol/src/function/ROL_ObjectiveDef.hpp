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

#ifndef ROL_OBJECTIVE_DEF_H
#define ROL_OBJECTIVE_DEF_H

#include <Teuchos_ScalarTraits.hpp>

/** \class ROL::Objective
    \brief Provides the definition of the objective function interface.
*/

namespace ROL {

template <class Real>
Real Objective<Real>::dirDeriv( const Vector<Real> &x, const Vector<Real> &d, const Real eta) {
  Teuchos::RCP<Vector<Real> > xd = d.clone();
  xd->set(x);
  xd->axpy(eta, d);
  return (this->value(*xd) - this->value(x)) / eta;
}

template <class Real>
void Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x ) {
  g.zero();
  Real deriv = 0.0;
  Real eta   = 1e-2*sqrt(Teuchos::ScalarTraits<Real>::eps());
  Real h     = 0.0;
  for (int i = 0; i < g.dimension(); i++) {
    h     = x.dot(*g.basis(i))*eta;
    deriv = dirDeriv(x,*g.basis(i),h);
    g.axpy(deriv,*g.basis(i));
  }
}

template <class Real>
void Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x ) {
  // Get Step Length
  Real eta = 1e-4*sqrt(Teuchos::ScalarTraits<Real>::eps());
  Real h   = x.dot(v)*eta;

  // Compute New Step x + h*v
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  
  // Compute Gradient at x
  Teuchos::RCP<Vector<Real> > g = x.clone();
  gradient(*g,x);

  // Compute Gradient at x + h*v
  hv.zero();
  gradient(hv,*xnew);
  
  // Compute Newton Quotient
  hv.axpy(-1.0,*g);
  hv.scale(1.0/h);
} 

} // namespace ROL

#endif
