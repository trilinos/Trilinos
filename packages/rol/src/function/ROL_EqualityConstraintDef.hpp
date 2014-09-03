// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_EQUALITYCONSTRAINT_DEF_H
#define ROL_EQUALITYCONSTRAINT_DEF_H

namespace ROL {

template <class Real>
void EqualityConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &x,
                                             Real &tol) {

  // By default we compute the finite-difference approximation.

  Real ctol = std::sqrt(ROL_EPSILON);

  // Get step length.
  Real h = std::max(1.0,x.norm()/v.norm())*tol;
  //Real h = 2.0/(v.norm()*v.norm())*tol;

  // Compute constraint at x.
  Teuchos::RCP<Vector<Real> > c = jv.clone();
  this->value(*c,x,ctol);

  // Compute perturbation x + h*v.
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  this->update(*xnew);

  // Compute constraint at x + h*v.
  jv.zero();
  this->value(jv,*xnew,ctol);

  // Compute Newton quotient.
  jv.axpy(-1.0,*c);
  jv.scale(1.0/h);
}


template <class Real>
void EqualityConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                                    const Vector<Real> &v,
                                                    const Vector<Real> &x,
                                                    Real &tol) {

  // By default we compute the finite-difference approximation.
  // This requires the implementation of a vector-space basis for the optimization variables.

  Real ctol = std::sqrt(ROL_EPSILON);

  Real h = 0.0;
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
  Teuchos::RCP<Vector<Real> > e    = x.clone();
  Teuchos::RCP<Vector<Real> > J    = v.clone();
  Teuchos::RCP<Vector<Real> > c    = v.clone();
  this->value(*c,x,ctol);
  ajv.zero();
  for ( unsigned i = 0; i < (unsigned)x.dimension(); i++ ) {
    e = x.basis(i);
    h = std::max(1.0,x.norm()/e->norm())*tol;
    xnew->set(x);
    xnew->axpy(h,*e);
    this->update(*xnew);
    this->value(*J,*xnew,ctol);
    J->axpy(-1.0,*c);
    J->scale(1.0/h);
    ajv.axpy(J->dot(v),*e);
  }
}


/*template <class Real>
void EqualityConstraint<Real>::applyHessian(Vector<Real> &huv,
                                            const Vector<Real> &u,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x,
                                            Real &tol ) {
  Real jtol = std::sqrt(ROL_EPSILON);

  // Get step length.
  Real h = std::max(1.0,x.norm()/v.norm())*tol;
  //Real h = 2.0/(v.norm()*v.norm())*tol;

  // Compute constraint Jacobian at x.
  Teuchos::RCP<Vector<Real> > ju = huv.clone();
  this->applyJacobian(*ju,u,x,jtol);

  // Compute new step x + h*v.
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
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
void EqualityConstraint<Real>::applyAdjointHessian(Vector<Real> &huv,
                                                   const Vector<Real> &u,
                                                   const Vector<Real> &v,
                                                   const Vector<Real> &x,
                                                   Real &tol ) {
}


template <class Real>
void EqualityConstraint<Real>::solveAugmentedSystem(Vector<Real> &v1,
                                                    Vector<Real> &v2,
                                                    const Vector<Real> &b1,
                                                    const Vector<Real> &b2,
                                                    const Vector<Real> &x,
                                                    Real &tol) {

  // Initialization.
  Real zero = 0.0;
  Real one  = 1.0;
  int m = 50;
  Real zerotol = zero;
  int iter = 0;
  int flag = 0;
  int i = 0;
  int k = 0;
  Real temp = zero;
  Real resnrm = zero;

  Teuchos::RCP<Vector<Real> > r1 = v1.clone();
  Teuchos::RCP<Vector<Real> > r2 = v2.clone();
  Teuchos::RCP<Vector<Real> > r2temp = v2.clone();
  Teuchos::RCP<Vector<Real> > v1trial = v1.clone();
  Teuchos::RCP<Vector<Real> > v2trial = v2.clone();
  Teuchos::RCP<Vector<Real> > z1 = v1.clone();
  Teuchos::RCP<Vector<Real> > z2 = v2.clone();
  Teuchos::RCP<Vector<Real> > w1 = v1.clone();
  Teuchos::RCP<Vector<Real> > w2 = v2.clone();
  Teuchos::RCP<Vector<Real> > w1temp = v1.clone();
  Teuchos::RCP<Vector<Real> > w2temp = v2.clone();
  std::vector<Real> res(m+1, zero);
  std::vector<Teuchos::RCP<Vector<Real> > > V1;
  std::vector<Teuchos::RCP<Vector<Real> > > V2;
  Teuchos::SerialDenseMatrix<int, Real> H(m+1,m);
  Teuchos::SerialDenseVector<int, Real> cs(m);
  Teuchos::SerialDenseVector<int, Real> sn(m);
  Teuchos::SerialDenseVector<int, Real> s(m+1);
  Teuchos::SerialDenseVector<int, Real> y(m);
  Teuchos::SerialDenseVector<int, Real> cnorm(m);
  Teuchos::LAPACK<int, Real> lapack;

  // Compute initial residual.
  applyAdjointJacobian(*r1, v2, x, zerotol);
  r1->scale(-one); r1->axpy(-one, v1); r1->plus(b1);
  applyJacobian(*r2, v1, x, zerotol);
  r2->scale(-one); r2->plus(b2);
  res[0] = std::sqrt(r1->dot(*r1) + r2->dot(*r2));

  // Check if residual is identically zero.
  if (res[0] == zero) {
    iter = 0;
    flag = 0;
    return;
  }

  // Apply left preconditioner to constraint residual.
  r2temp->set(*r2);
  applyPreconditioner(*r2, *r2temp, x, zerotol);

  // Compute norm of preconditioned residual.
  res[0] = std::sqrt(r1->dot(*r1) + r2->dot(*r2));

  // Evaluate special stopping condition and check convergence.
  tol = tol;
  if ( res[0] <= tol ) {
    return;
  }

  V1.push_back(r1->clone()); (V1[0])->set(*r1); (V1[0])->scale(one/res[0]);
  V2.push_back(r2->clone()); (V2[0])->set(*r2); (V2[0])->scale(one/res[0]);

  s(0) = res[0];

  for (i=0; i<m; i++) {

    w1->set(*(V1[i]));
    applyAdjointJacobian(*w1temp, *(V2[i]), x, zerotol);
    applyJacobian(*w2temp, *w1, x, zerotol);
    w1->plus(*w1temp);
    applyPreconditioner(*w2, *w2temp, x, zerotol);
    
    // Evaluate coefficients and orthogonalize using Gram-Schmidt.
    for (k=0; k<i+1; k++) {
      H(k,i) = w1->dot(*(V1[k])) + w2->dot(*(V2[k]));
      w1->axpy(-H(k,i), *(V1[k]));
      w2->axpy(-H(k,i), *(V2[k]));
    } // for (int k=0; k<i; k++)
    H(i+1,i) = std::sqrt(w1->dot(*w1) + w2->dot(*w2));
    V1.push_back(w1->clone()); (V1[i+1])->set(*w1); (V1[i+1])->scale(one/H(i+1,i));
    V2.push_back(w2->clone()); (V2[i+1])->set(*w2); (V2[i+1])->scale(one/H(i+1,i));

    // Apply Givens rotations.
    for (k=0; k<i; k++) {
      temp     = cs(k)*H(k,i) + sn(k)*H(k+1,i);
      H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
      H(k,i)   = temp;
    } // for (int k=0; k<i-1; k++)

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
    lapack.LATRS(uplo, trans, diag, normin, i, H.values(), m, y.values(), &scaling, cnorm.values(), &info);
    z1->zero();
    z2->zero();
    for (k=0; k<i; k++) {
      z1->axpy(y(k), *(V1[k]));
      z2->axpy(y(k), *(V2[k]));
    }

    // Compute true residual.
    v1trial->set(v1); v1trial->plus(*z1);
    v2trial->set(v2); v2trial->plus(*z2);
    applyAdjointJacobian(*r1, *v2trial, x, zerotol);
    r1->scale(-one); r1->axpy(-one, *v1trial); r1->plus(b1);
    applyJacobian(*r2, *v1trial, x, zerotol);
    r2->scale(-one); r2->plus(b2);
    resnrm = std::sqrt(r1->dot(*r1) + r2->dot(*r2));

    // Evaluate special stopping condition.
    tol = tol;

std::cout << "resnrm = " << resnrm << "\n";
std::cout << "tol = " << tol << "\n";

    if (resnrm <= tol) {
      // Update solution vector.
      v1.plus(*z1);
      v2.plus(*z2);
      break;
    }

  } // for (int i=0; i++; i<m)

  for (int j=0; j<m; j++) {
    std::cout << "res(" << j << ")" << res[j] << "\n";
  }

}


template <class Real>
std::vector<std::vector<Real> > EqualityConstraint<Real>::checkApplyJacobian(const Vector<Real> &x,
                                                                             const Vector<Real> &v,
                                                                             const Vector<Real> &jv,
                                                                             const bool printToScreen,
                                                                             const int numSteps) {
  Real tol = std::sqrt(ROL_EPSILON);

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > jvCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = 1.0;

  std::ios::fmtflags f( std::cout.flags() );

  // Compute constraint value at x.
  Teuchos::RCP<Vector<Real> > c = jv.clone();
  this->value(*c, x, tol);

  // Compute (Jacobian at x) times (vector v).
  Teuchos::RCP<Vector<Real> > Jv = jv.clone();
  this->applyJacobian(*Jv, v, x, tol);
  Real normJv = Jv->norm();

  // Temporary vectors.
  Teuchos::RCP<Vector<Real> > cnew = jv.clone();
  Teuchos::RCP<Vector<Real> > xnew = x.clone();

  for (int i=0; i<numSteps; i++) {
    // Evaluate constraint value at x+eta*v.
    xnew->set(x);
    xnew->axpy(eta, v);
    this->update(*xnew);
    this->value(*cnew, *xnew, tol);
    cnew->axpy(-1.0, *c);
    cnew->scale(1.0/eta);

    // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
    jvCheck[i][0] = eta;
    jvCheck[i][1] = normJv;
    jvCheck[i][2] = cnew->norm();
    cnew->axpy(-1.0, *Jv);
    jvCheck[i][3] = cnew->norm();

    if (printToScreen) {
      if (i==0) {
      std::cout << std::right
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
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << jvCheck[i][0]
                << std::setw(20) << jvCheck[i][1]
                << std::setw(20) << jvCheck[i][2]
                << std::setw(20) << jvCheck[i][3]
                << "\n";
    }

    // Update eta.
    eta = eta*eta_factor;
  }

  std::cout.flags( f );

  return jvCheck;
} // checkApplyJacobian


template <class Real>
std::vector<std::vector<Real> > EqualityConstraint<Real>::checkApplyAdjointJacobian(const Vector<Real> &x,
                                                                                    const Vector<Real> &v,
                                                                                    const bool printToScreen,
                                                                                    const int numSteps) {
  Real tol = std::sqrt(ROL_EPSILON);

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > ajvCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = 1.0;

  std::ios::fmtflags f( std::cout.flags() );

  // Compute constraint value at x.
  Teuchos::RCP<Vector<Real> > c = v.clone();
  this->value(*c, x, tol);

  // Compute (Jacobian at x) times (vector v).
  Teuchos::RCP<Vector<Real> > AJv = x.clone();
  this->applyAdjointJacobian(*AJv, v, x, tol);
  Real normAJv = AJv->norm();

  // Temporary vectors.
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
  Teuchos::RCP<Vector<Real> > e    = x.clone();
  Teuchos::RCP<Vector<Real> > cnew = v.clone();
  Teuchos::RCP<Vector<Real> > ajv  = x.clone();

  for (int i=0; i<numSteps; i++) {

    ajv->zero();
    for ( unsigned j = 0; j < (unsigned)x.dimension(); j++ ) {
      e = x.basis(j);
      xnew->set(x);
      xnew->axpy(eta,*e);
      this->update(*xnew);
      this->value(*cnew,*xnew,tol);
      cnew->axpy(-1.0,*c);
      cnew->scale(1.0/eta);
      ajv->axpy(cnew->dot(v),*e);
    }

    // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
    ajvCheck[i][0] = eta;
    ajvCheck[i][1] = normAJv;
    ajvCheck[i][2] = ajv->norm();
    ajv->axpy(-1.0, *AJv);
    ajvCheck[i][3] = ajv->norm();

    if (printToScreen) {
      if (i==0) {
      std::cout << std::right
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
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << ajvCheck[i][0]
                << std::setw(20) << ajvCheck[i][1]
                << std::setw(20) << ajvCheck[i][2]
                << std::setw(20) << ajvCheck[i][3]
                << "\n";
    }

    // Update eta.
    eta = eta*eta_factor;
  }

  std::cout.flags( f );

  return ajvCheck;
} // checkApplyAdjointJacobian


template <class Real>
std::vector<std::vector<Real> > EqualityConstraint<Real>::checkApplyAdjointHessian(const Vector<Real> &x,
                                                                                   const Vector<Real> &u,
                                                                                   const Vector<Real> &v,
                                                                                   const bool printToScreen,
                                                                                   const int numSteps) {
  Real tol = std::sqrt(ROL_EPSILON);

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real> > ahuvCheck(numSteps, tmp);
  Real eta_factor = 1e-1;
  Real eta = 1.0;

  std::ios::fmtflags f( std::cout.flags() );

  // Apply adjoint Jacobian to u.
  Teuchos::RCP<Vector<Real> > AJu = v.clone();
  this->applyAdjointJacobian(*AJu, u, x, tol);

  // Apply adjoint Hessian at x, in direction v, to u.
  Teuchos::RCP<Vector<Real> > AHuv = v.clone();
  this->applyAdjointHessian(*AHuv, u, v, x, tol);
  Real normAHuv = AHuv->norm();

  // Temporary vectors.
  Teuchos::RCP<Vector<Real> > AJnew = v.clone();
  Teuchos::RCP<Vector<Real> > xnew = x.clone();

  for (int i=0; i<numSteps; i++) {
    // Apply adjoint Jacobian to u at x+eta*v.
    xnew->set(x);
    xnew->axpy(eta, v);
    this->update(*xnew);
    this->applyAdjointJacobian(*AJnew, u, *xnew, tol);
    AJnew->axpy(-1.0, *AJu);
    AJnew->scale(1.0/eta);

    // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
    ahuvCheck[i][0] = eta;
    ahuvCheck[i][1] = normAHuv;
    ahuvCheck[i][2] = AJnew->norm();
    AJnew->axpy(-1.0, *AHuv);
    ahuvCheck[i][3] = AJnew->norm();

    if (printToScreen) {
      if (i==0) {
      std::cout << std::right
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
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << ahuvCheck[i][0]
                << std::setw(20) << ahuvCheck[i][1]
                << std::setw(20) << ahuvCheck[i][2]
                << std::setw(20) << ahuvCheck[i][3]
                << "\n";
    }

    // Update eta.
    eta = eta*eta_factor;
  }

  std::cout.flags( f );

  return ahuvCheck;
} // checkApplyAdjointHessian

} // namespace ROL

#endif
