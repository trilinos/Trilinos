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

template <class Real>
std::vector<std::vector<Real> > Objective<Real>::checkGradient( const Vector<Real> &x,
                                                                const Vector<Real> &d,
                                                                const bool printToScreen ) {
  int numSteps = 13;
  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > gCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = 1.0;

  std::ios::fmtflags f( std::cout.flags() );

  // Compute gradient at x.
  Teuchos::RCP<Vector<Real> > g = x.clone();
  this->gradient(*g, x);

  // Temporary vectors.
  Teuchos::RCP<Vector<Real> > xnew = x.clone();

  // Evaluate objective value at x.
  Real fval_at_x = this->value(x);

  Real fval_at_xnew = 0;
  for (int i=0; i<numSteps; i++) {
    // Evaluate objective value at x+eta*d.
    xnew->set(x);
    xnew->axpy(eta, d);
    fval_at_xnew = this->value(*xnew);

    // Compute gradient, finite-difference gradient, and absolute error.
    gCheck[i][0] = eta;
    gCheck[i][1] = d.dot(*g);
    gCheck[i][2] = (fval_at_xnew - fval_at_x) / eta;
    gCheck[i][3] = gCheck[i][2] - gCheck[i][1];

    if (printToScreen) {
      if (i==0) {
      std::cout << std::right
                << std::setw(20) << "Step size"
                << std::setw(20) << "grad'*dir"
                << std::setw(20) << "FD approx"
                << std::setw(20) << "abs error"
                << "\n";
      }
      std::cout << std::scientific << std::setprecision(8) << std::right
                << std::setw(20) << gCheck[i][0]
                << std::setw(20) << gCheck[i][1]
                << std::setw(20) << gCheck[i][2]
                << std::setw(20) << gCheck[i][3]
                << "\n";
    }

    // Update eta.
    eta = eta*eta_factor;
  }

  std::cout.flags( f );

  return gCheck;
} // checkGradient


} // namespace ROL

#endif
