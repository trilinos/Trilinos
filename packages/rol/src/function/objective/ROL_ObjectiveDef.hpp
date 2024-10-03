// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OBJECTIVE_DEF_H
#define ROL_OBJECTIVE_DEF_H

/** \class ROL::Objective
    \brief Provides the definition of the objective function interface.
*/

namespace ROL {

template<typename Real>
Real Objective<Real>::dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol) {
  if (dual_ == nullPtr) dual_ = x.dual().clone();
  gradient(*dual_,x,tol);
  //return d.dot(dual_->dual());
  return d.apply(*dual_);
  //Real dnorm = d.norm(), zero(0);
  //if ( dnorm == zero ) {
  //  return zero;
  //}
  //Real cbrteps = std::cbrt(ROL_EPSILON<Real>()), one(1), v0(0), v1(0);
  //Real xnorm = x.norm(), h = cbrteps * std::max(xnorm/dnorm,one);
  //v0 = value(x,tol);
  //prim_->set(x); prim_->axpy(h, d);
  //update(*prim_,UpdateType::Temp);
  //v1 = value(*prim_,tol);
  //update(x,UpdateType::Revert);
  //return (v1 - v0) / h;
}

template<typename Real>
void Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  if (prim_ == nullPtr) prim_ = x.clone();
  if (basis_ == nullPtr) basis_ = x.clone();

  const Real cbrteps = std::cbrt(ROL_EPSILON<Real>()), zero(0), one(1);
  Real f0 = value(x,tol), h(0), xi(0), gi(0);
  g.zero();
  for (int i = 0; i < x.dimension(); i++) {
    basis_->set(*x.basis(i));
    xi = x.dot(*basis_);
    h  = cbrteps * std::max(std::abs(xi),one) * (xi < zero ? -one : one);
    prim_->set(x); prim_->axpy(h,*basis_);
    h  = prim_->dot(*basis_) - xi;
    update(*prim_,UpdateType::Temp);
    gi = (value(*prim_,tol) - f0) / h;
    g.axpy(gi,*g.basis(i));
  }
  update(x,UpdateType::Revert);
}

template<typename Real>
void Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  const Real zero(0), vnorm = v.norm();
  // Get Step Length
  if ( vnorm == zero ) {
    hv.zero();
  }
  else {
    if (prim_ == nullPtr) prim_ = x.clone();
    if (dual_ == nullPtr) dual_ = hv.clone();

    //Real h = 2.0/(v.norm()*v.norm())*tol;
    const Real one(1), h(std::max(one,x.norm()/vnorm)*tol);

    gradient(*dual_,x,tol);           // Compute gradient at x
    prim_->set(x); prim_->axpy(h,v);  // Set prim = x + hv
    update(*prim_,UpdateType::Temp);       // Temporarily update objective at x + hv
    gradient(hv,*prim_,tol);          // Compute gradient at x + hv
    hv.axpy(-one,*dual_);             // Compute difference (f'(x+hv)-f'(x))
    hv.scale(one/h);                  // Compute Newton quotient (f'(x+hv)-f'(x))/h
    update(x,UpdateType::Revert);          // Reset objective to x
  }
}

template<typename Real>
std::vector<std::vector<Real>> Objective<Real>::checkGradient( const Vector<Real> &x,
                                                               const Vector<Real> &g,
                                                               const Vector<Real> &d,
                                                               const bool printToStream,
                                                               std::ostream & outStream,
                                                               const int numSteps,
                                                               const int order ) {

  const Real ten(10);
  std::vector<Real> steps(numSteps);
  for(int i=0;i<numSteps;++i) {
    steps[i] = pow(ten,static_cast<Real>(-i));
  }

  return checkGradient(x,g,d,steps,printToStream,outStream,order);

} // checkGradient

template<typename Real>
std::vector<std::vector<Real>> Objective<Real>::checkGradient( const Vector<Real> &x,
                                                               const Vector<Real> &g,
                                                               const Vector<Real> &d,
                                                               const std::vector<Real> &steps,
                                                               const bool printToStream,
                                                               std::ostream & outStream,
                                                               const int order ) {

  ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                              "Error: finite difference order must be 1,2,3, or 4" );

  using Finite_Difference_Arrays::shifts;
  using Finite_Difference_Arrays::weights;

  Real tol = std::sqrt(ROL_EPSILON<Real>());

  int numSteps = steps.size();
  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real>> gCheck(numSteps, tmp);

  // Save the format state of the original outStream.
  nullstream oldFormatState;
  oldFormatState.copyfmt(outStream);

  // Evaluate objective value at x.
  update(x,UpdateType::Temp);
  Real val = value(x,tol);

  // Compute gradient at x.
  Ptr<Vector<Real>> gtmp = g.clone();
  gradient(*gtmp, x, tol);
  //Real dtg = d.dot(gtmp->dual());
  Real dtg = d.apply(*gtmp);

  // Temporary vectors.
  Ptr<Vector<Real>> xnew = x.clone();

  for (int i=0; i<numSteps; i++) {

    Real eta = steps[i];

    xnew->set(x);

    // Compute gradient, finite-difference gradient, and absolute error.
    gCheck[i][0] = eta;
    gCheck[i][1] = dtg;

    gCheck[i][2] = weights[order-1][0] * val;

    for(int j=0; j<order; ++j) {
      // Evaluate at x <- x+eta*c_i*d.
      xnew->axpy(eta*shifts[order-1][j], d);

      // Only evaluate at shifts where the weight is nonzero  
      if( weights[order-1][j+1] != 0 ) {
        update(*xnew,UpdateType::Temp);
        gCheck[i][2] += weights[order-1][j+1] * this->value(*xnew,tol);
      }
    }

    gCheck[i][2] /= eta;

    gCheck[i][3] = std::abs(gCheck[i][2] - gCheck[i][1]);

    if (printToStream) {
      if (i==0) {
        outStream << std::right
                  << std::setw(20) << "Step size"
                  << std::setw(20) << "grad'*dir"
                  << std::setw(20) << "FD approx"
                  << std::setw(20) << "abs error"
                  << "\n"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "---------"
                  << "\n";
      }
      outStream << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << gCheck[i][0]
                << std::setw(20) << gCheck[i][1]
                << std::setw(20) << gCheck[i][2]
                << std::setw(20) << gCheck[i][3]
                << "\n";
    }

  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return gCheck;
} // checkGradient

template<typename Real>
std::vector<std::vector<Real>> Objective<Real>::checkHessVec( const Vector<Real> &x,
                                                              const Vector<Real> &hv,
                                                              const Vector<Real> &v,
                                                              const bool printToStream,
                                                              std::ostream & outStream,
                                                              const int numSteps,
                                                              const int order ) {
  const Real ten(10);
  std::vector<Real> steps(numSteps);
  for(int i=0;i<numSteps;++i) {
    steps[i] = pow(ten,static_cast<Real>(-i));
  }

  return checkHessVec(x,hv,v,steps,printToStream,outStream,order);
} // checkHessVec



template<typename Real>
std::vector<std::vector<Real>> Objective<Real>::checkHessVec( const Vector<Real> &x,
                                                              const Vector<Real> &hv,
                                                              const Vector<Real> &v,
                                                              const std::vector<Real> &steps,
                                                              const bool printToStream,
                                                              std::ostream & outStream,
                                                              const int order ) {

  ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                              "Error: finite difference order must be 1,2,3, or 4" );

  using Finite_Difference_Arrays::shifts;
  using Finite_Difference_Arrays::weights;

  const Real one(1);
  Real tol = std::sqrt(ROL_EPSILON<Real>());

  int numSteps = steps.size();
  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real>> hvCheck(numSteps, tmp);

  // Save the format state of the original outStream.
  nullstream oldFormatState;
  oldFormatState.copyfmt(outStream);

  // Compute gradient at x.
  Ptr<Vector<Real>> g = hv.clone();
  update(x,UpdateType::Temp);
  gradient(*g, x, tol);

  // Compute (Hessian at x) times (vector v).
  Ptr<Vector<Real>> Hv = hv.clone();
  hessVec(*Hv, v, x, tol);
  Real normHv = Hv->norm();

  // Temporary vectors.
  Ptr<Vector<Real>> gdif = hv.clone();
  Ptr<Vector<Real>> gnew = hv.clone();
  Ptr<Vector<Real>> xnew = x.clone();

  for (int i=0; i<numSteps; i++) {
    Real eta = steps[i]; 
    // Evaluate objective value at x+eta*d.
    xnew->set(x);
    gdif->set(*g);
    gdif->scale(weights[order-1][0]);
    for (int j=0; j<order; ++j) {
      // Evaluate at x <- x+eta*c_i*d.
      xnew->axpy(eta*shifts[order-1][j], v);
      // Only evaluate at shifts where the weight is nonzero  
      if ( weights[order-1][j+1] != 0 ) {
        update(*xnew,UpdateType::Temp);
        gradient(*gnew, *xnew, tol); 
        gdif->axpy(weights[order-1][j+1],*gnew);
      }
    }
    gdif->scale(one/eta);    

    // Compute norms of hessvec, finite-difference hessvec, and error.
    hvCheck[i][0] = eta;
    hvCheck[i][1] = normHv;
    hvCheck[i][2] = gdif->norm();
    gdif->axpy(-one, *Hv);
    hvCheck[i][3] = gdif->norm();

    if (printToStream) {
      if (i==0) {
      outStream << std::right
                << std::setw(20) << "Step size"
                << std::setw(20) << "norm(Hess*vec)"
                << std::setw(20) << "norm(FD approx)"
                << std::setw(20) << "norm(abs error)"
                << "\n"
                << std::setw(20) << "---------"
                << std::setw(20) << "--------------"
                << std::setw(20) << "---------------"
                << std::setw(20) << "---------------"
                << "\n";
      }
      outStream << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << hvCheck[i][0]
                << std::setw(20) << hvCheck[i][1]
                << std::setw(20) << hvCheck[i][2]
                << std::setw(20) << hvCheck[i][3]
                << "\n";
    }

  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return hvCheck;
} // checkHessVec

template<typename Real>
std::vector<Real> Objective<Real>::checkHessSym( const Vector<Real> &x,
                                                 const Vector<Real> &hv,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &w,
                                                 const bool printToStream,
                                                 std::ostream & outStream ) {

  Real tol = std::sqrt(ROL_EPSILON<Real>());
  
  // Compute (Hessian at x) times (vector v).
  Ptr<Vector<Real>> h = hv.clone();
  update(x,UpdateType::Temp);
  hessVec(*h, v, x, tol);
  //Real wHv = w.dot(h->dual());
  Real wHv = w.apply(*h);

  hessVec(*h, w, x, tol);
  //Real vHw = v.dot(h->dual());
  Real vHw = v.apply(*h);

  std::vector<Real> hsymCheck(3, 0);

  hsymCheck[0] = wHv;
  hsymCheck[1] = vHw;
  hsymCheck[2] = std::abs(vHw-wHv);

  // Save the format state of the original outStream.
  nullstream oldFormatState;
  oldFormatState.copyfmt(outStream);

  if (printToStream) {
    outStream << std::right
              << std::setw(20) << "<w, H(x)v>"
              << std::setw(20) << "<v, H(x)w>"
              << std::setw(20) << "abs error"
              << "\n";
    outStream << std::scientific << std::setprecision(11) << std::right
              << std::setw(20) << hsymCheck[0]
              << std::setw(20) << hsymCheck[1]
              << std::setw(20) << hsymCheck[2]
              << "\n";
  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return hsymCheck;

} // checkHessSym

} // namespace ROL

#endif
