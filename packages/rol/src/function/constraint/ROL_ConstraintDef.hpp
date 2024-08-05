// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_DEF_H
#define ROL_CONSTRAINT_DEF_H

#include "ROL_LinearAlgebra.hpp"
#include "ROL_LAPACK.hpp"

namespace ROL {

template <class Real>
void Constraint<Real>::applyJacobian(Vector<Real> &jv,
                                     const Vector<Real> &v,
                                     const Vector<Real> &x,
                                     Real &tol) {
  // By default we compute the finite-difference approximation.
  const Real one(1);
  Real ctol = std::sqrt(ROL_EPSILON<Real>());

  // Get step length.
  Real h = std::max(one,x.norm()/v.norm())*tol;
  //Real h = 2.0/(v.norm()*v.norm())*tol;

  // Compute constraint at x.
  ROL::Ptr<Vector<Real> > c = jv.clone();
  this->value(*c,x,ctol);

  // Compute perturbation x + h*v.
  ROL::Ptr<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  this->update(*xnew,UpdateType::Temp);

  // Compute constraint at x + h*v.
  jv.zero();
  this->value(jv,*xnew,ctol);

  // Compute Newton quotient.
  jv.axpy(-one,*c);
  jv.scale(one/h);
}


template <class Real>
void Constraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x,
                                            Real &tol) {
  applyAdjointJacobian(ajv,v,x,v.dual(),tol);
}





template <class Real>
void Constraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x,
                                            const Vector<Real> &dualv,
                                            Real &tol) { 
  // By default we compute the finite-difference approximation.
  // This requires the implementation of a vector-space basis for the optimization variables.
  // The default implementation requires that the constraint space is equal to its dual.
  const Real one(1);
  Real h(0), ctol = std::sqrt(ROL_EPSILON<Real>());

  ROL::Ptr<Vector<Real> > xnew = x.clone();
  ROL::Ptr<Vector<Real> > ex   = x.clone();
  ROL::Ptr<Vector<Real> > eajv = ajv.clone();
  ROL::Ptr<Vector<Real> > cnew = dualv.clone();  // in general, should be in the constraint space
  ROL::Ptr<Vector<Real> > c0   = dualv.clone();  // in general, should be in the constraint space
  this->value(*c0,x,ctol);
  
  ajv.zero();
  for ( int i = 0; i < ajv.dimension(); i++ ) {
    ex = x.basis(i);
    eajv = ajv.basis(i);
    h = std::max(one,x.norm()/ex->norm())*tol;
    xnew->set(x);
    xnew->axpy(h,*ex);
    this->update(*xnew,UpdateType::Temp);
    this->value(*cnew,*xnew,ctol);
    cnew->axpy(-one,*c0);
    cnew->scale(one/h);
    //ajv.axpy(cnew->dot(v.dual()),*eajv);
    ajv.axpy(cnew->apply(v),*eajv);
  }
}


/*template <class Real>
void Constraint<Real>::applyHessian(Vector<Real> &huv,
                                    const Vector<Real> &u,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    Real &tol ) {
  Real jtol = std::sqrt(ROL_EPSILON<Real>());

  // Get step length.
  Real h = std::max(1.0,x.norm()/v.norm())*tol;
  //Real h = 2.0/(v.norm()*v.norm())*tol;

  // Compute constraint Jacobian at x.
  ROL::Ptr<Vector<Real> > ju = huv.clone();
  this->applyJacobian(*ju,u,x,jtol);

  // Compute new step x + h*v.
  ROL::Ptr<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  this->update(*xnew);

  // Compute constraint Jacobian at x + h*v.
  huv.zero();
  this->applyJacobian(huv,u,*xnew,jtol);

  // Compute Newton quotient.
  huv.axpy(-1.0,*ju);
  huv.scale(1.0/h);
}*/


template <class Real>
void Constraint<Real>::applyAdjointHessian(Vector<Real> &huv,
                                           const Vector<Real> &u,
                                           const Vector<Real> &v,
                                           const Vector<Real> &x,
                                           Real &tol ) {
  // Get step length.
  Real h = std::max(static_cast<Real>(1),x.norm()/v.norm())*tol;

  // Compute constraint Jacobian at x.
  ROL::Ptr<Vector<Real> > aju = huv.clone();
  applyAdjointJacobian(*aju,u,x,tol);

  // Compute new step x + h*v.
  ROL::Ptr<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  update(*xnew,UpdateType::Temp);

  // Compute constraint Jacobian at x + h*v.
  huv.zero();
  applyAdjointJacobian(huv,u,*xnew,tol);

  // Compute Newton quotient.
  huv.axpy(static_cast<Real>(-1),*aju);
  huv.scale(static_cast<Real>(1)/h);
}


template <class Real>
std::vector<Real> Constraint<Real>::solveAugmentedSystem(Vector<Real> &v1,
                                                         Vector<Real> &v2,
                                                         const Vector<Real> &b1,
                                                         const Vector<Real> &b2,
                                                         const Vector<Real> &x,
                                                         Real &tol) {

  /*** Initialization. ***/
  const Real zero(0), one(1);
  int m = 200;           // Krylov space size.
  Real zerotol = zero;
  int i = 0;
  int k = 0;
  Real temp = zero;
  Real resnrm = zero;

  //tol = std::sqrt(b1.dot(b1)+b2.dot(b2))*1e-8;
  tol = std::sqrt(b1.dot(b1)+b2.dot(b2))*tol;

  // Set initial guess to zero.
  v1.zero(); v2.zero();

  // Allocate static memory.
  ROL::Ptr<Vector<Real> > r1 = b1.clone();
  ROL::Ptr<Vector<Real> > r2 = b2.clone();
  ROL::Ptr<Vector<Real> > z1 = v1.clone();
  ROL::Ptr<Vector<Real> > z2 = v2.clone();
  ROL::Ptr<Vector<Real> > w1 = b1.clone();
  ROL::Ptr<Vector<Real> > w2 = b2.clone();
  std::vector<ROL::Ptr<Vector<Real> > > V1;
  std::vector<ROL::Ptr<Vector<Real> > > V2;
  ROL::Ptr<Vector<Real> > V2temp = b2.clone();
  std::vector<ROL::Ptr<Vector<Real> > > Z1;
  std::vector<ROL::Ptr<Vector<Real> > > Z2;
  ROL::Ptr<Vector<Real> > w1temp = b1.clone();
  ROL::Ptr<Vector<Real> > Z2temp = v2.clone();

  std::vector<Real> res(m+1, zero); 
  LA::Matrix<Real> H(m+1,m);
  LA::Vector<Real> cs(m);
  LA::Vector<Real> sn(m);
  LA::Vector<Real> s(m+1);
  LA::Vector<Real> y(m+1);
  LA::Vector<Real> cnorm(m);
  ROL::LAPACK<int, Real> lapack;

  // Compute initial residual.
  applyAdjointJacobian(*r1, v2, x, zerotol);
  r1->scale(-one); r1->axpy(-one, v1.dual()); r1->plus(b1);
  applyJacobian(*r2, v1, x, zerotol);
  r2->scale(-one); r2->plus(b2);
  res[0] = std::sqrt(r1->dot(*r1) + r2->dot(*r2));

  // Check if residual is identically zero.
  if (res[0] == zero) {
    res.resize(0);
    return res;
  }

  V1.push_back(b1.clone()); (V1[0])->set(*r1); (V1[0])->scale(one/res[0]);
  V2.push_back(b2.clone()); (V2[0])->set(*r2); (V2[0])->scale(one/res[0]);

  s(0) = res[0];

  for (i=0; i<m; i++) {

    // Apply right preconditioner.
    V2temp->set(*(V2[i]));
    applyPreconditioner(*Z2temp, *V2temp, x, b1, zerotol);
    Z2.push_back(v2.clone()); (Z2[i])->set(*Z2temp);
    Z1.push_back(v1.clone()); (Z1[i])->set((V1[i])->dual());

    // Apply operator.
    applyJacobian(*w2, *(Z1[i]), x, zerotol);
    applyAdjointJacobian(*w1temp, *Z2temp, x, zerotol);
    w1->set(*(V1[i])); w1->plus(*w1temp);
    
    // Evaluate coefficients and orthogonalize using Gram-Schmidt.
    for (k=0; k<=i; k++) {
      H(k,i) = w1->dot(*(V1[k])) + w2->dot(*(V2[k]));
      w1->axpy(-H(k,i), *(V1[k]));
      w2->axpy(-H(k,i), *(V2[k]));
    }
    H(i+1,i) = std::sqrt(w1->dot(*w1) + w2->dot(*w2));
    
    V1.push_back(b1.clone()); (V1[i+1])->set(*w1); (V1[i+1])->scale(one/H(i+1,i));
    V2.push_back(b2.clone()); (V2[i+1])->set(*w2); (V2[i+1])->scale(one/H(i+1,i));

    // Apply Givens rotations.
    for (k=0; k<=i-1; k++) {
      temp     = cs(k)*H(k,i) + sn(k)*H(k+1,i);
      H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
      H(k,i)   = temp;
    }

    // Form i-th rotation matrix.
    if ( H(i+1,i) == zero ) {
      cs(i) = one;
      sn(i) = zero;
    }
    else if ( std::abs(H(i+1,i)) > std::abs(H(i,i)) ) {
      temp = H(i,i) / H(i+1,i);
      sn(i) = one / std::sqrt( one + temp*temp );
      cs(i) = temp * sn(i);
    }
    else {
      temp = H(i+1,i) / H(i,i);
      cs(i) = one / std::sqrt( one + temp*temp );
      sn(i) = temp * cs(i);
    }

    // Approximate residual norm.
    temp     = cs(i)*s(i);
    s(i+1)   = -sn(i)*s(i);
    s(i)     = temp;
    H(i,i)   = cs(i)*H(i,i) + sn(i)*H(i+1,i);
    H(i+1,i) = zero;
    resnrm   = std::abs(s(i+1));
    res[i+1] = resnrm;

    // Update solution approximation.
    const char uplo = 'U';
    const char trans = 'N';
    const char diag = 'N';
    const char normin = 'N';
    Real scaling = zero;
    int info = 0;
    y = s;
    lapack.LATRS(uplo, trans, diag, normin, i+1, H.values(), m+1, y.values(), &scaling, cnorm.values(), &info);
    z1->zero();
    z2->zero();
    for (k=0; k<=i; k++) {
      z1->axpy(y(k), *(Z1[k]));
      z2->axpy(y(k), *(Z2[k]));
    }

    // Evaluate special stopping condition.
    //tol = ???;

//    std::cout << "  " << i+1 << ": " << res[i+1]/res[0] << std::endl;
    if (res[i+1] <= tol) {
//      std::cout << "  solved in " << i+1 << " iterations to " << res[i+1] << " (" << res[i+1]/res[0] << ")" << std::endl;
      // Update solution vector.
      v1.plus(*z1);
      v2.plus(*z2);
      break;
    }

  } // for (int i=0; i++; i<m)

  res.resize(i+2);

  /*
  std::stringstream hist;
  hist << std::scientific << std::setprecision(8);
  hist << "\n    Augmented System Solver:\n";
  hist << "    Iter Residual\n";
  for (unsigned j=0; j<res.size(); j++) {
    hist << "    " << std::left << std::setw(14) << res[j] << "\n";
  }
  hist << "\n";
  std::cout << hist.str();
  */

  return res;
}


template <class Real>
std::vector<std::vector<Real> > Constraint<Real>::checkApplyJacobian(const Vector<Real> &x,
                                                                     const Vector<Real> &v,
                                                                     const Vector<Real> &jv,
                                                                     const bool printToStream,
                                                                     std::ostream & outStream,
                                                                     const int numSteps,
                                                                     const int order) {
  std::vector<Real> steps(numSteps);
  for(int i=0;i<numSteps;++i) {
    steps[i] = pow(10,-i);
  }
 
  return checkApplyJacobian(x,v,jv,steps,printToStream,outStream,order);
}




template <class Real>
std::vector<std::vector<Real> > Constraint<Real>::checkApplyJacobian(const Vector<Real> &x,
                                                                     const Vector<Real> &v,
                                                                     const Vector<Real> &jv,
                                                                     const std::vector<Real> &steps, 
                                                                     const bool printToStream,
                                                                     std::ostream & outStream,
                                                                     const int order) {
  ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                              "Error: finite difference order must be 1,2,3, or 4" );

  const Real one(1.0);

  using Finite_Difference_Arrays::shifts;
  using Finite_Difference_Arrays::weights;

  Real tol = std::sqrt(ROL_EPSILON<Real>());

  int numSteps = steps.size();
  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > jvCheck(numSteps, tmp);

  // Save the format state of the original outStream.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(outStream);

  // Compute constraint value at x.
  ROL::Ptr<Vector<Real> > c = jv.clone();
  this->update(x,UpdateType::Temp);
  this->value(*c, x, tol);

  // Compute (Jacobian at x) times (vector v).
  ROL::Ptr<Vector<Real> > Jv = jv.clone();
  this->applyJacobian(*Jv, v, x, tol);
  Real normJv = Jv->norm();

  // Temporary vectors.
  ROL::Ptr<Vector<Real> > cdif = jv.clone();
  ROL::Ptr<Vector<Real> > cnew = jv.clone();
  ROL::Ptr<Vector<Real> > xnew = x.clone();

  for (int i=0; i<numSteps; i++) {

    Real eta = steps[i];

    xnew->set(x);
 
    cdif->set(*c);
    cdif->scale(weights[order-1][0]);

    for(int j=0; j<order; ++j) {

       xnew->axpy(eta*shifts[order-1][j], v);

       if( weights[order-1][j+1] != 0 ) {
           this->update(*xnew,UpdateType::Temp);
           this->value(*cnew,*xnew,tol);
           cdif->axpy(weights[order-1][j+1],*cnew);    
       }

    }

    cdif->scale(one/eta);    

    // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
    jvCheck[i][0] = eta;
    jvCheck[i][1] = normJv;
    jvCheck[i][2] = cdif->norm();
    cdif->axpy(-one, *Jv);
    jvCheck[i][3] = cdif->norm();

    if (printToStream) {
      std::stringstream hist;
      if (i==0) {
      hist << std::right
           << std::setw(20) << "Step size"
           << std::setw(20) << "norm(Jac*vec)"
           << std::setw(20) << "norm(FD approx)"
           << std::setw(20) << "norm(abs error)"
           << "\n"
           << std::setw(20) << "---------"
           << std::setw(20) << "-------------"
           << std::setw(20) << "---------------"
           << std::setw(20) << "---------------"
           << "\n";
      }
      hist << std::scientific << std::setprecision(11) << std::right
           << std::setw(20) << jvCheck[i][0]
           << std::setw(20) << jvCheck[i][1]
           << std::setw(20) << jvCheck[i][2]
           << std::setw(20) << jvCheck[i][3]
           << "\n";
      outStream << hist.str();
    }

  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return jvCheck;
} // checkApplyJacobian


template <class Real>
std::vector<std::vector<Real> > Constraint<Real>::checkApplyAdjointJacobian(const Vector<Real> &x,
                                                                            const Vector<Real> &v,
                                                                            const Vector<Real> &c,
                                                                            const Vector<Real> &ajv,
                                                                            const bool printToStream,
                                                                            std::ostream & outStream,
                                                                            const int numSteps) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  const Real one(1);

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > ajvCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = one;

  // Temporary vectors.
  ROL::Ptr<Vector<Real> > c0   = c.clone();
  ROL::Ptr<Vector<Real> > cnew = c.clone();
  ROL::Ptr<Vector<Real> > xnew = x.clone();
  ROL::Ptr<Vector<Real> > ajv0 = ajv.clone();
  ROL::Ptr<Vector<Real> > ajv1 = ajv.clone();
  ROL::Ptr<Vector<Real> > ex   = x.clone();
  ROL::Ptr<Vector<Real> > eajv = ajv.clone();

  // Save the format state of the original outStream.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(outStream);

  // Compute constraint value at x.
  this->update(x,UpdateType::Temp);
  this->value(*c0, x, tol);

  // Compute (Jacobian at x) times (vector v).
  this->applyAdjointJacobian(*ajv0, v, x, tol);
  Real normAJv = ajv0->norm();

  for (int i=0; i<numSteps; i++) {

    ajv1->zero();

    for ( int j = 0; j < ajv.dimension(); j++ ) {
      ex = x.basis(j);
      eajv = ajv.basis(j);
      xnew->set(x);
      xnew->axpy(eta,*ex);
      this->update(*xnew,UpdateType::Temp);
      this->value(*cnew,*xnew,tol);
      cnew->axpy(-one,*c0);
      cnew->scale(one/eta);
      //ajv1->axpy(cnew->dot(v.dual()),*eajv);
      ajv1->axpy(cnew->apply(v),*eajv);
    }

    // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
    ajvCheck[i][0] = eta;
    ajvCheck[i][1] = normAJv;
    ajvCheck[i][2] = ajv1->norm();
    ajv1->axpy(-one, *ajv0);
    ajvCheck[i][3] = ajv1->norm();

    if (printToStream) {
      std::stringstream hist;
      if (i==0) {
      hist << std::right
           << std::setw(20) << "Step size"
           << std::setw(20) << "norm(adj(Jac)*vec)"
           << std::setw(20) << "norm(FD approx)"
           << std::setw(20) << "norm(abs error)"
           << "\n"
           << std::setw(20) << "---------"
           << std::setw(20) << "------------------"
           << std::setw(20) << "---------------"
           << std::setw(20) << "---------------"
           << "\n";
      }
      hist << std::scientific << std::setprecision(11) << std::right
           << std::setw(20) << ajvCheck[i][0]
           << std::setw(20) << ajvCheck[i][1]
           << std::setw(20) << ajvCheck[i][2]
           << std::setw(20) << ajvCheck[i][3]
           << "\n";
      outStream << hist.str();
    }

    // Update eta.
    eta = eta*eta_factor;
  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return ajvCheck;
} // checkApplyAdjointJacobian

template <class Real>
Real Constraint<Real>::checkAdjointConsistencyJacobian(const Vector<Real> &w,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &x,
                                                       const Vector<Real> &dualw,
                                                       const Vector<Real> &dualv,
                                                       const bool printToStream,
                                                       std::ostream & outStream) {
  Real tol = ROL_EPSILON<Real>();

  ROL::Ptr<Vector<Real> > Jv = dualw.clone();
  ROL::Ptr<Vector<Real> > Jw = dualv.clone();
  
  this->update(x,UpdateType::Temp);
  applyJacobian(*Jv,v,x,tol);
  applyAdjointJacobian(*Jw,w,x,tol);

  //Real vJw = v.dot(Jw->dual());
  Real vJw = v.apply(*Jw);
  //Real wJv = w.dot(Jv->dual());
  Real wJv = w.apply(*Jv);

  Real diff = std::abs(wJv-vJw);

  if ( printToStream ) {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(8);
    hist << "\nTest Consistency of Jacobian and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = " 
         << diff << "\n";
    hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
    hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW<Real>()) << "\n";
    outStream << hist.str();
  }
  return diff;
} // checkAdjointConsistencyJacobian

template <class Real>
std::vector<std::vector<Real> > Constraint<Real>::checkApplyAdjointHessian(const Vector<Real> &x,
                                                                           const Vector<Real> &u,
                                                                           const Vector<Real> &v,
                                                                           const Vector<Real> &hv,
                                                                           const bool printToStream,
                                                                           std::ostream & outStream,
                                                                           const int numSteps,  
                                                                           const int order) {
  std::vector<Real> steps(numSteps);
  for(int i=0;i<numSteps;++i) {
    steps[i] = pow(10,-i);
  }
 
  return checkApplyAdjointHessian(x,u,v,hv,steps,printToStream,outStream,order);
}


template <class Real>
std::vector<std::vector<Real> > Constraint<Real>::checkApplyAdjointHessian(const Vector<Real> &x,
                                                                           const Vector<Real> &u,
                                                                           const Vector<Real> &v,
                                                                           const Vector<Real> &hv,
                                                                           const std::vector<Real> &steps,  
                                                                           const bool printToStream,
                                                                           std::ostream & outStream,
                                                                           const int order) {
  using Finite_Difference_Arrays::shifts;
  using Finite_Difference_Arrays::weights;

  const Real one(1);
  Real tol = std::sqrt(ROL_EPSILON<Real>());

  int numSteps = steps.size();
  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > ahuvCheck(numSteps, tmp);

  // Temporary vectors.
  ROL::Ptr<Vector<Real> > AJdif = hv.clone();
  ROL::Ptr<Vector<Real> > AJu = hv.clone();
  ROL::Ptr<Vector<Real> > AHuv = hv.clone();
  ROL::Ptr<Vector<Real> > AJnew = hv.clone();
  ROL::Ptr<Vector<Real> > xnew = x.clone();

  // Save the format state of the original outStream.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(outStream);

  // Apply adjoint Jacobian to u.
  this->update(x,UpdateType::Temp);
  this->applyAdjointJacobian(*AJu, u, x, tol);

  // Apply adjoint Hessian at x, in direction v, to u.
  this->applyAdjointHessian(*AHuv, u, v, x, tol);
  Real normAHuv = AHuv->norm();

  for (int i=0; i<numSteps; i++) {

    Real eta = steps[i];

    // Apply adjoint Jacobian to u at x+eta*v.
    xnew->set(x);

    AJdif->set(*AJu);
    AJdif->scale(weights[order-1][0]);     

    for(int j=0; j<order; ++j) {

        xnew->axpy(eta*shifts[order-1][j],v); 

        if( weights[order-1][j+1] != 0 ) {    
            this->update(*xnew,UpdateType::Temp);
            this->applyAdjointJacobian(*AJnew, u, *xnew, tol);
            AJdif->axpy(weights[order-1][j+1],*AJnew);
        }
    }

    AJdif->scale(one/eta);

    // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
    ahuvCheck[i][0] = eta;
    ahuvCheck[i][1] = normAHuv;
    ahuvCheck[i][2] = AJdif->norm();
    AJdif->axpy(-one, *AHuv);
    ahuvCheck[i][3] = AJdif->norm();

    if (printToStream) {
      std::stringstream hist;
      if (i==0) {
      hist << std::right
           << std::setw(20) << "Step size"
           << std::setw(20) << "norm(adj(H)(u,v))"
           << std::setw(20) << "norm(FD approx)"
           << std::setw(20) << "norm(abs error)"
           << "\n"
           << std::setw(20) << "---------"
           << std::setw(20) << "-----------------"
           << std::setw(20) << "---------------"
           << std::setw(20) << "---------------"
           << "\n";
      }
      hist << std::scientific << std::setprecision(11) << std::right
           << std::setw(20) << ahuvCheck[i][0]
           << std::setw(20) << ahuvCheck[i][1]
           << std::setw(20) << ahuvCheck[i][2]
           << std::setw(20) << ahuvCheck[i][3]
           << "\n";
      outStream << hist.str();
    }

  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return ahuvCheck;
} // checkApplyAdjointHessian

} // namespace ROL

#endif
