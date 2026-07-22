// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLE_U_AS_DEF_H
#define ROL_BUNDLE_U_AS_DEF_H

namespace ROL {

template<typename Real>
Bundle_U_AS<Real>::Bundle_U_AS(const unsigned maxSize,
            const Real coeff,
            const Real omega,
            const unsigned remSize) 
  : Bundle_U<Real>(maxSize,coeff,omega,remSize), isInitialized_(false) {}

template<typename Real>
void Bundle_U_AS<Real>::initialize(const Vector<Real> &g) {
  Bundle_U<Real>::initialize(g);
  if ( !isInitialized_ ) {
    tG_ = g.clone();
    yG_ = g.clone();
    eG_ = g.clone();
    gx_ = g.clone();
    ge_ = g.clone();
    isInitialized_ = true;
  }
}

template<typename Real>
unsigned Bundle_U_AS<Real>::solveDual(const Real t, const unsigned maxit, const Real tol) {
  unsigned iter = 0;
  if (Bundle_U<Real>::size() == 1) {
    iter = Bundle_U<Real>::solveDual_dim1(t,maxit,tol);
  }
  else if (Bundle_U<Real>::size() == 2) {
    iter = Bundle_U<Real>::solveDual_dim2(t,maxit,tol);
  }
  else {
    iter = solveDual_arbitrary(t,maxit,tol);
  }
  return iter;
}

template<typename Real>
void Bundle_U_AS<Real>::initializeDualSolver(void) {
  const Real zero(0);
  Real sum(0), err(0), tmp(0), y(0);
  for (unsigned i = 0; i < Bundle_U<Real>::size(); ++i) {
    // Compute sum of dualVariables_ using Kahan's compensated sum
    //sum += Bundle_U<Real>::getDualVariable(i);
    y   = Bundle_U<Real>::getDualVariable(i) - err;
    tmp = sum + y;
    err = (tmp - sum) - y;
    sum = tmp;
  }
  for (unsigned i = 0; i < Bundle_U<Real>::size(); ++i) {
    tmp = Bundle_U<Real>::getDualVariable(i)/sum;
    Bundle_U<Real>::setDualVariable(i,tmp);
  }
  nworkingSet_.clear();
  workingSet_.clear();
  for (unsigned i = 0; i < Bundle_U<Real>::size(); ++i) {
    if ( Bundle_U<Real>::getDualVariable(i) == zero ) {
      workingSet_.insert(i);
    }
    else {
      nworkingSet_.insert(i);
    }
  }
}

template<typename Real>
void Bundle_U_AS<Real>::computeLagMult(std::vector<Real> &lam, const Real mu, const std::vector<Real> &g) const {
  const Real zero(0);
  unsigned n = workingSet_.size();
  if ( n > 0 ) {
    lam.resize(n,zero);
    typename std::set<unsigned>::iterator it = workingSet_.begin();
    for (unsigned i = 0; i < n; ++i) {
      lam[i] = g[*it] - mu; it++;
    }
  }
  else {
    lam.clear();
  }
}
 
template<typename Real>
bool Bundle_U_AS<Real>::isNonnegative(unsigned &ind, const std::vector<Real> &x) const {
  bool flag = true;
  unsigned n = workingSet_.size();
  ind = Bundle_U<Real>::size();
  if ( n > 0 ) {
    Real min = ROL_OVERFLOW<Real>();
    typename std::set<unsigned>::iterator it = workingSet_.begin();
    for (unsigned i = 0; i < n; ++i) {
      if ( x[i] < min ) {
        ind = *it;
        min = x[i];
      }
      it++;
    }
    flag = ((min < -ROL_EPSILON<Real>()) ? false : true);
  }
  return flag;
}

template<typename Real>
Real Bundle_U_AS<Real>::computeStepSize(unsigned &ind, const std::vector<Real> &x, const std::vector<Real> &p) const {
  const Real zero(0);
  Real alpha(1), tmp(0);
  ind = Bundle_U<Real>::size();
  typename std::set<unsigned>::iterator it;
  for (it = nworkingSet_.begin(); it != nworkingSet_.end(); it++) {
    if ( p[*it] < -ROL_EPSILON<Real>() ) {
      tmp = -x[*it]/p[*it];
      if ( alpha >= tmp ) {
        alpha = tmp;
        ind = *it;
      }
    }
  }
  return std::max(zero,alpha);
}

template<typename Real>
unsigned Bundle_U_AS<Real>::solveEQPsubproblem(std::vector<Real> &s, Real &mu,
                      const std::vector<Real> &g, const Real tol) const {
  // Build reduced QP information
  const Real zero(0);
  unsigned n = nworkingSet_.size(), cnt = 0;
  mu = zero;
  s.assign(Bundle_U<Real>::size(),zero);
  if ( n > 0 ) {
    std::vector<Real> gk(n,zero);
    typename std::set<unsigned>::iterator it = nworkingSet_.begin();
    for (unsigned i = 0; i < n; ++i) {
      gk[i] = g[*it]; it++;
    }
    std::vector<Real> sk(n,zero);
    cnt = projectedCG(sk,mu,gk,tol);
    it  = nworkingSet_.begin();
    for (unsigned i = 0; i < n; ++i) {
      s[*it] = sk[i]; it++;
    }
  }
  return cnt;
}

template<typename Real>
void Bundle_U_AS<Real>::applyPreconditioner(std::vector<Real> &Px, const std::vector<Real> &x) const {
  const Real zero(0);
  int type = 0;
  std::vector<Real> tmp(Px.size(),zero);
  switch (type) {
    case 0: applyPreconditioner_Identity(tmp,x); break;
    case 1: applyPreconditioner_Jacobi(tmp,x);   break;
    case 2: applyPreconditioner_SymGS(tmp,x);    break;
  }
  applyPreconditioner_Identity(Px,tmp);
}

template<typename Real>
void Bundle_U_AS<Real>::applyG(std::vector<Real> &Gx, const std::vector<Real> &x) const {
  int type = 0;
  switch (type) {
    case 0: applyG_Identity(Gx,x); break;
    case 1: applyG_Jacobi(Gx,x);   break;
    case 2: applyG_SymGS(Gx,x);    break;
  }
}

template<typename Real>
void Bundle_U_AS<Real>::applyPreconditioner_Identity(std::vector<Real> &Px, const std::vector<Real> &x) const {
  unsigned dim = nworkingSet_.size();
  Real sum(0), err(0), tmp(0), y(0);
  for (unsigned i = 0; i < dim; ++i) {
    // Compute sum of x using Kahan's compensated sum
    //sum += x[i];
    y   = x[i] - err;
    tmp = sum + y;
    err = (tmp - sum) - y;
    sum = tmp;
  }
  sum /= (Real)dim;
  for (unsigned i = 0; i < dim; ++i) {
    Px[i] = x[i] - sum;
  }
}

template<typename Real>
void Bundle_U_AS<Real>::applyG_Identity(std::vector<Real> &Gx, const std::vector<Real> &x) const {
  Gx.assign(x.begin(),x.end());
}

template<typename Real>
void Bundle_U_AS<Real>::applyPreconditioner_Jacobi(std::vector<Real> &Px, const std::vector<Real> &x) const {
  const Real zero(0), one(1);
  unsigned dim = nworkingSet_.size();
  Real eHe(0), sum(0);
  Real errX(0), tmpX(0), yX(0), errE(0), tmpE(0), yE(0);
  std::vector<Real> gg(dim,zero);
  typename std::set<unsigned>::iterator it = nworkingSet_.begin(); 
  for (unsigned i = 0; i < dim; ++i) {
    gg[i] = one/std::abs(Bundle_U<Real>::GiGj(*it,*it)); it++;
    // Compute sum of inv(D)x using Kahan's aggregated sum
    //sum += x[i]*gg[i];
    yX   = x[i]*gg[i] - errX;
    tmpX = sum + yX;
    errX = (tmpX - sum) - yX;
    sum  = tmpX;
    // Compute sum of inv(D)e using Kahan's aggregated sum
    //eHe += gg[i];
    yE   = gg[i] - errE;
    tmpE = eHe + yE;
    errE = (tmpE - eHe) - yE;
    eHe  = tmpE;
  }
  sum /= eHe;
  for (unsigned i = 0; i < dim; ++i) {
    Px[i] = (x[i]-sum)*gg[i];
  }
}

template<typename Real>
void Bundle_U_AS<Real>::applyG_Jacobi(std::vector<Real> &Gx, const std::vector<Real> &x) const {
  unsigned dim = nworkingSet_.size();
  typename std::set<unsigned>::iterator it = nworkingSet_.begin();
  for (unsigned i = 0; i < dim; ++i) {
    Gx[i] = std::abs(Bundle_U<Real>::GiGj(*it,*it))*x[i]; it++;
  }
}

template<typename Real>
void Bundle_U_AS<Real>::applyPreconditioner_SymGS(std::vector<Real> &Px, const std::vector<Real> &x) const {
  int dim = nworkingSet_.size();
  //unsigned cnt = 0;
  gx_->zero(); ge_->zero();
  Real eHx(0), eHe(0), one(1);
  // Forward substitution
  std::vector<Real> x1(dim,0), e1(dim,0),gg(dim,0);
  typename std::set<unsigned>::iterator it, jt;
  it = nworkingSet_.begin(); 
  for (int i = 0; i < dim; ++i) {
    gx_->zero(); ge_->zero(); jt = nworkingSet_.begin();
    for (int j = 0; j < i; ++j) {
      Bundle_U<Real>::addGi(*jt,x1[j],*gx_);
      Bundle_U<Real>::addGi(*jt,e1[j],*ge_);
      jt++;
    }
    gg[i] = Bundle_U<Real>::GiGj(*it,*it);
    x1[i] = (x[i] - Bundle_U<Real>::dotGi(*it,*gx_))/gg[i];
    e1[i] = (one  - Bundle_U<Real>::dotGi(*it,*ge_))/gg[i];
    it++;
  }
  // Apply diagonal
  for (int i = 0; i < dim; ++i) {
    x1[i] *= gg[i];
    e1[i] *= gg[i];
  }
  // Back substitution
  std::vector<Real> Hx(dim,0), He(dim,0); it = nworkingSet_.end();
  for (int i = dim-1; i >= 0; --i) {
    it--;
    gx_->zero(); ge_->zero(); jt = nworkingSet_.end();
    for (int j = dim-1; j >= i+1; --j) {
      jt--;
      Bundle_U<Real>::addGi(*jt,Hx[j],*gx_);
      Bundle_U<Real>::addGi(*jt,He[j],*ge_);
    }
    Hx[i] = (x1[i] - Bundle_U<Real>::dotGi(*it,*gx_))/gg[i];
    He[i] = (e1[i] - Bundle_U<Real>::dotGi(*it,*ge_))/gg[i];
    // Compute sums
    eHx += Hx[i];
    eHe += He[i];
  }
  // Accumulate preconditioned vector
  for (int i = 0; i < dim; ++i) {
    Px[i] = Hx[i] - (eHx/eHe)*He[i];
  }
}

template<typename Real>
void Bundle_U_AS<Real>::applyG_SymGS(std::vector<Real> &Gx, const std::vector<Real> &x) const {
  unsigned dim = nworkingSet_.size();
  typename std::set<unsigned>::iterator it = nworkingSet_.begin();
  for (unsigned i = 0; i < dim; ++i) {
    Gx[i] = std::abs(Bundle_U<Real>::GiGj(*it,*it))*x[i]; it++;
  }
}

template<typename Real>
void Bundle_U_AS<Real>::computeResidualUpdate(std::vector<Real> &r, std::vector<Real> &g) const {
  unsigned n = g.size();
  std::vector<Real> Gg(n,0);
  Real y(0), ytmp(0), yprt(0), yerr(0);
  applyPreconditioner(g,r);
  applyG(Gg,g);
  // Compute multiplier using Kahan's compensated sum
  for (unsigned i = 0; i < n; ++i) {
    //y += (r[i] - Gg[i]);
    yprt = (r[i] - Gg[i]) - yerr;
    ytmp = y + yprt;
    yerr = (ytmp - y) - yprt;
    y    = ytmp;
  }
  y /= (Real)n;
  for (unsigned i = 0; i < n; ++i) {
    r[i] -= y;
  }
  applyPreconditioner(g,r);
}

template<typename Real>
void Bundle_U_AS<Real>::applyFullMatrix(std::vector<Real> &Hx, const std::vector<Real> &x) const {
  const Real one(1);
  gx_->zero(); eG_->zero();
  for (unsigned i = 0; i < Bundle_U<Real>::size(); ++i) {
    // Compute Gx using Kahan's compensated sum
    //gx_->axpy(x[i],Bundle_U<Real>::subgradient(i));
    yG_->set(Bundle_U<Real>::subgradient(i)); yG_->scale(x[i]); yG_->axpy(-one,*eG_);
    tG_->set(*gx_); tG_->plus(*yG_);
    eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->axpy(-one,*yG_);
    gx_->set(*tG_);
  }
  for (unsigned i = 0; i < Bundle_U<Real>::size(); ++i) {
    // Compute < g_i, Gx >
    Hx[i] = Bundle_U<Real>::dotGi(i,*gx_);
  }
}

template<typename Real>
void Bundle_U_AS<Real>::applyMatrix(std::vector<Real> &Hx, const std::vector<Real> &x) const {
  const Real one(1);
  gx_->zero(); eG_->zero();
  unsigned n = nworkingSet_.size();
  typename std::set<unsigned>::iterator it = nworkingSet_.begin(); 
  for (unsigned i = 0; i < n; ++i) {
    // Compute Gx using Kahan's compensated sum
    //gx_->axpy(x[i],Bundle_U<Real>::subgradient(*it));
    yG_->set(Bundle_U<Real>::subgradient(*it)); yG_->scale(x[i]); yG_->axpy(-one,*eG_);
    tG_->set(*gx_); tG_->plus(*yG_);
    eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->axpy(-one,*yG_);
    gx_->set(*tG_);
    it++;
  }
  it = nworkingSet_.begin();
  for (unsigned i = 0; i < n; ++i) {
    // Compute < g_i, Gx >
    Hx[i] = Bundle_U<Real>::dotGi(*it,*gx_); it++;
  }
}

template<typename Real>
unsigned Bundle_U_AS<Real>::projectedCG(std::vector<Real> &x, Real &mu, const std::vector<Real> &b, const Real tol) const {
  const Real one(1), zero(0);
  unsigned n = nworkingSet_.size();
  std::vector<Real> r(n,0), r0(n,0), g(n,0), d(n,0), Ad(n,0);
  // Compute residual Hx + g = g with x = 0
  x.assign(n,0);
  scale(r,one,b);
  r0.assign(r.begin(),r.end());
  // Precondition residual
  computeResidualUpdate(r,g);
  Real rg = dot(r,g), rg0(0);
  // Get search direction
  scale(d,-one,g);
  Real alpha(0), kappa(0), beta(0), TOL(1.e-2);
  Real CGtol = std::min(tol,TOL*rg);
  unsigned cnt = 0;
  while (rg > CGtol && cnt < 2*n+1) {
    applyMatrix(Ad,d);
    kappa = dot(d,Ad);
    alpha = rg/kappa;
    axpy(alpha,d,x);
    axpy(alpha,Ad,r);
    axpy(alpha,Ad,r0);
    computeResidualUpdate(r,g);
    rg0 = rg;
    rg  = dot(r,g);
    beta = rg/rg0;
    scale(d,beta);
    axpy(-one,g,d);
    cnt++;
  }
  // Compute multiplier for equality constraint using Kahan's compensated sum
  mu = zero;
  Real err(0), tmp(0), y(0);
  for (unsigned i = 0; i < n; ++i) {
    //mu += r0[i];
    y   = r0[i] - err;
    tmp = mu + y;
    err = (tmp - mu) - y;
    mu  = tmp;
  }
  mu /= static_cast<Real>(n);
  // Return iteration count
  return cnt;
}

template<typename Real>
Real Bundle_U_AS<Real>::dot(const std::vector<Real> &x, const std::vector<Real> &y) const {
  // Compute dot product of two vectors using Kahan's compensated sum
  Real val(0), err(0), tmp(0), y0(0);
  unsigned n = std::min(x.size(),y.size());
  for (unsigned i = 0; i < n; ++i) {
    //val += x[i]*y[i];
    y0  = x[i]*y[i] - err;
    tmp = val + y0;
    err = (tmp - val) - y0;
    val = tmp;
  }
  return val;
}

template<typename Real>
Real Bundle_U_AS<Real>::norm(const std::vector<Real> &x) const {
  return std::sqrt(dot(x,x));
}

template<typename Real>
void Bundle_U_AS<Real>::axpy(const Real a, const std::vector<Real> &x, std::vector<Real> &y) const {
  unsigned n = std::min(y.size(),x.size());
  for (unsigned i = 0; i < n; ++i) {
    y[i] += a*x[i];
  }
}

template<typename Real>
void Bundle_U_AS<Real>::scale(std::vector<Real> &x, const Real a) const {
  for (unsigned i = 0; i < x.size(); ++i) {
    x[i] *= a;
  }
}

template<typename Real>
void Bundle_U_AS<Real>::scale(std::vector<Real> &x, const Real a, const std::vector<Real> &y) const {
  unsigned n = std::min(x.size(),y.size());
  for (unsigned i = 0; i < n; ++i) {
    x[i] = a*y[i];
  }
}

template<typename Real>
unsigned Bundle_U_AS<Real>::solveDual_arbitrary(const Real t, const unsigned maxit, const Real tol) {
  const Real zero(0), one(1);
  initializeDualSolver();
  bool nonneg = false;
  unsigned ind = 0, i = 0, CGiter = 0;
  Real snorm(0), alpha(0), mu(0);
  std::vector<Real> s(Bundle_U<Real>::size(),0), Hs(Bundle_U<Real>::size(),0);
  std::vector<Real> g(Bundle_U<Real>::size(),0), lam(Bundle_U<Real>::size()+1,0);
  std::vector<Real> dualVariables(Bundle_U<Real>::size(),0);
  for (unsigned j = 0; j < Bundle_U<Real>::size(); ++j) {
    dualVariables[j] = Bundle_U<Real>::getDualVariable(j);
  }
  //Real val = Bundle_U<Real>::evaluateObjective(g,dualVariables,t);
  Bundle_U<Real>::evaluateObjective(g,dualVariables,t);
  for (i = 0; i < maxit; ++i) {
    CGiter += solveEQPsubproblem(s,mu,g,tol);
    snorm = norm(s);
    if ( snorm < ROL_EPSILON<Real>() ) {
      computeLagMult(lam,mu,g);
      nonneg = isNonnegative(ind,lam);
      if ( nonneg ) {
        break;
      }
      else {
        alpha = one;
        if ( ind < Bundle_U<Real>::size() ) {
          workingSet_.erase(ind);
          nworkingSet_.insert(ind);
        }
      }
    }
    else {
      alpha = computeStepSize(ind,dualVariables,s);
      if ( alpha > zero ) {
        axpy(alpha,s,dualVariables);
        applyFullMatrix(Hs,s);
        axpy(alpha,Hs,g);
      }
      if (ind < Bundle_U<Real>::size()) {
        workingSet_.insert(ind);
        nworkingSet_.erase(ind);
      }
    }
    //std::cout << "iter = " << i << "  snorm = " << snorm << "  alpha = " << alpha << "\n";
  }
  //Real crit = computeCriticality(g,dualVariables);
  //std::cout << "Criticality Measure: " << crit << "\n";
  //std::cout << "dim = " << Bundle_U<Real>::size() << "  iter = " << i << "   CGiter = " << CGiter << "  CONVERGED!\n";
  for (unsigned j = 0; j < Bundle_U<Real>::size(); ++j) {
    Bundle_U<Real>::setDualVariable(j,dualVariables[j]);
  }
  return i;
}

template<typename Real>
void Bundle_U_AS<Real>::project(std::vector<Real> &x, const std::vector<Real> &v) const {
  const Real zero(0), one(1);
  std::vector<Real> vsort(Bundle_U<Real>::size(),0);
  vsort.assign(v.begin(),v.end());
  std::sort(vsort.begin(),vsort.end());
  Real sum(-1), lam(0);
  for (unsigned i = Bundle_U<Real>::size()-1; i > 0; i--) {
    sum += vsort[i];
    if ( sum >= static_cast<Real>(Bundle_U<Real>::size()-i)*vsort[i-1] ) {
      lam = sum/static_cast<Real>(Bundle_U<Real>::size()-i);
      break;
    }
  }
  if (lam == zero) {
    lam = (sum+vsort[0])/static_cast<Real>(Bundle_U<Real>::size());
  }
  for (unsigned i = 0; i < Bundle_U<Real>::size(); ++i) {
    x[i] = std::max(zero,v[i] - lam);
  }
}

template<typename Real>
Real Bundle_U_AS<Real>::computeCriticality(const std::vector<Real> &g, const std::vector<Real> &sol) const {
  const Real zero(0), one(1);
  std::vector<Real> x(Bundle_U<Real>::size(),0), Px(Bundle_U<Real>::size(),0);
  axpy(one,sol,x);
  axpy(-one,g,x);
  project(Px,x);
  scale(x,zero);
  axpy(one,sol,x);
  axpy(-one,Px,x);
  return norm(x);
}

} // namespace ROL

#endif
