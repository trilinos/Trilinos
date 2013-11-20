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

/** \class ROL::Objective
    \brief Provides the definition of the objective function interface.
*/

namespace ROL {

template <class Real>
Real Objective<Real>::dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol) {
  Real ftol = std::sqrt(ROL_EPSILON);

  Teuchos::RCP<Vector<Real> > xd = d.clone();
  xd->set(x);
  xd->axpy(tol, d);
  return (this->value(*xd,ftol) - this->value(x,ftol)) / tol;
}

template <class Real>
void Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  g.zero();
  Real deriv = 0.0;
  Real h     = 0.0;
  for (int i = 0; i < g.dimension(); i++) {
    h     = x.dot(*g.basis(i))*tol;
    deriv = this->dirDeriv(x,*g.basis(i),h);
    g.axpy(deriv,*g.basis(i));
  }
}

template <class Real>
void Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  Real gtol = std::sqrt(ROL_EPSILON);

  // Get Step Length
  Real h = std::max(1.0,x.norm()/v.norm())*tol;

  // Compute New Step x + h*v
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  
  // Compute Gradient at x
  Teuchos::RCP<Vector<Real> > g = x.clone();
  this->gradient(*g,x,gtol);

  // Compute Gradient at x + h*v
  hv.zero();
  this->gradient(hv,*xnew,gtol);
  
  // Compute Newton Quotient
  hv.axpy(-1.0,*g);
  hv.scale(1.0/h);
} 

template <class Real>
std::vector<std::vector<Real> > Objective<Real>::checkGradient( const Vector<Real> &x,
                                                                const Vector<Real> &d,
                                                                const bool printToScreen,
                                                                const int numSteps ) {
  Real tol = std::sqrt(ROL_EPSILON);

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > gCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = 1.0;

  std::ios::fmtflags f( std::cout.flags() );

  // Compute gradient at x.
  Teuchos::RCP<Vector<Real> > g = x.clone();
  this->gradient(*g, x, tol);
  Real dtg = d.dot(*g);

  // Temporary vectors.
  Teuchos::RCP<Vector<Real> > xnew = x.clone();

  // Evaluate objective value at x.
  Real fval_at_x = this->value(x,tol);

  Real fval_at_xnew = 0;
  for (int i=0; i<numSteps; i++) {
    // Evaluate objective value at x+eta*d.
    xnew->set(x);
    xnew->axpy(eta, d);
    fval_at_xnew = this->value(*xnew,tol);

    // Compute gradient, finite-difference gradient, and absolute error.
    gCheck[i][0] = eta;
    gCheck[i][1] = dtg;
    gCheck[i][2] = (fval_at_xnew - fval_at_x) / eta;
    gCheck[i][3] = std::abs(gCheck[i][2] - gCheck[i][1]);

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

template <class Real>
std::vector<std::vector<Real> > Objective<Real>::checkHessVec( const Vector<Real> &x,
                                                               const Vector<Real> &v,
                                                               const bool printToScreen,
                                                               const int numSteps ) {
  Real tol = std::sqrt(ROL_EPSILON);

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > hvCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = 1.0;

  std::ios::fmtflags f( std::cout.flags() );

  // Compute gradient at x.
  Teuchos::RCP<Vector<Real> > g = x.clone();
  this->gradient(*g, x, tol);

  // Compute (Hessian at x) times (vector v).
  Teuchos::RCP<Vector<Real> > Hv = x.clone();
  this->hessVec(*Hv, v, x, tol);
  Real normHv = Hv->norm();

  // Temporary vectors.
  Teuchos::RCP<Vector<Real> > gnew = x.clone();
  Teuchos::RCP<Vector<Real> > xnew = x.clone();

  for (int i=0; i<numSteps; i++) {
    // Evaluate objective value at x+eta*d.
    xnew->set(x);
    xnew->axpy(eta, v);
    this->gradient(*gnew, *xnew, tol);
    gnew->axpy(-1.0, *g);
    gnew->scale(1.0/eta);

    // Compute norms of hessvec, finite-difference hessvec, and error.
    hvCheck[i][0] = eta;
    hvCheck[i][1] = normHv;
    hvCheck[i][2] = gnew->norm();
    gnew->axpy(-1.0, *Hv);
    hvCheck[i][3] = gnew->norm();

    if (printToScreen) {
      if (i==0) {
      std::cout << std::right
                << std::setw(20) << "Step size"
                << std::setw(20) << "norm(Hess*vec)"
                << std::setw(20) << "norm(FD approx)"
                << std::setw(20) << "norm(abs error)"
                << "\n";
      }
      std::cout << std::scientific << std::setprecision(8) << std::right
                << std::setw(20) << hvCheck[i][0]
                << std::setw(20) << hvCheck[i][1]
                << std::setw(20) << hvCheck[i][2]
                << std::setw(20) << hvCheck[i][3]
                << "\n";
    }

    // Update eta.
    eta = eta*eta_factor;
  }

  std::cout.flags( f );

  return hvCheck;
} // checkHessVec


} // namespace ROL

#endif
