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

#ifndef ROL_BUNDLE_H
#define ROL_BUNDLE_H

#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"

#include "Teuchos_RCP.hpp"

#include <vector>
#include <set>

/** \class ROL::Bundle
    \brief Provides the interface for and implments a bundle.
*/

namespace ROL {

template<class Real>
class Bundle {
/***********************************************************************************************/
/***************** BUNDLE STORAGE **************************************************************/
/***********************************************************************************************/
private: 
  std::vector<Teuchos::RCP<Vector<Real> > > subgradients_;
  std::vector<Real> linearizationErrors_;
  std::vector<Real> distanceMeasures_;

  std::vector<Real> dualVariables_;

  unsigned size_;

  unsigned maxSize_;
  unsigned remSize_;
  Real coeff_;

  bool isInitialized_;

  void remove(const std::vector<unsigned> &ind) {
    Real zero(0);
    for (unsigned j = ind.back()+1; j < size_; j++) {
      (subgradients_[j-1])->set(*(subgradients_[j]));
      linearizationErrors_[j-1] = linearizationErrors_[j];
      distanceMeasures_[j-1]    = distanceMeasures_[j];
      dualVariables_[j-1]       = dualVariables_[j];
    }
    (subgradients_[size_-1])->zero();
    linearizationErrors_[size_-1] = ROL_OVERFLOW<Real>();
    distanceMeasures_[size_-1]    = ROL_OVERFLOW<Real>();
    dualVariables_[size_-1]       = zero;
    for (unsigned i = ind.size()-1; i > 0; --i) {
      for (unsigned j = ind[i-1]+1; j < size_; j++) {
        (subgradients_[j-1])->set(*(subgradients_[j]));
        linearizationErrors_[j-1] = linearizationErrors_[j];
        distanceMeasures_[j-1]    = distanceMeasures_[j];
        dualVariables_[j-1]       = dualVariables_[j];
      }
    }
    size_ -= ind.size();
  }

  void add(const Vector<Real> &g, const Real le, const Real dm) {
    Real zero(0);
    (subgradients_[size_])->set(g);
    linearizationErrors_[size_] = le;
    distanceMeasures_[size_]    = dm;
    dualVariables_[size_]       = zero;
    size_++;
  }

protected:
  Teuchos::RCP<Vector<Real> > tG_;
  Teuchos::RCP<Vector<Real> > eG_;
  Teuchos::RCP<Vector<Real> > yG_;
  
/***********************************************************************************************/
/***************** BUNDLE MODIFICATION AND ACCESS ROUTINES *************************************/
/***********************************************************************************************/
public:
  virtual ~Bundle(void) {}
  Bundle(const unsigned maxSize = 10, const Real coeff = 0.0, const unsigned remSize = 2) 
    : size_(0), maxSize_(maxSize), remSize_(remSize), coeff_(coeff), isInitialized_(false) {
    Real zero(0);
    remSize_ = ((remSize_ < 2) ? 2 : ((remSize_ > maxSize_-1) ? maxSize_-1 : remSize_));
    subgradients_.clear();
    subgradients_.assign(maxSize,Teuchos::null);
    linearizationErrors_.clear();
    linearizationErrors_.assign(maxSize_,ROL_OVERFLOW<Real>());
    distanceMeasures_.clear();
    distanceMeasures_.assign(maxSize_,ROL_OVERFLOW<Real>());
    dualVariables_.clear();
    dualVariables_.assign(maxSize_,zero);
  }

  void initialize(const Vector<Real> &g) {
    if ( !isInitialized_ ) {
      Real zero(0), one(1);
      for (unsigned i = 0; i < maxSize_; i++) {
        subgradients_[i] = g.clone();
      }
      (subgradients_[0])->set(g);
      linearizationErrors_[0] = zero;
      distanceMeasures_[0]    = zero;
      dualVariables_[0]       = one;
      size_++;
      isInitialized_ = true;
      gx_ = g.clone();
      ge_ = g.clone();
      tG_ = g.clone();
      yG_ = g.clone();
      eG_ = g.clone();
    }
  }

  const Real linearizationError(const unsigned i) const {
    return linearizationErrors_[i];
  }

  const Real distanceMeasure(const unsigned i) const {
    return distanceMeasures_[i];
  }

  const Vector<Real> & subgradient(const unsigned i) const {
    return *(subgradients_[i]);
  }

  const Real computeAlpha(const Real dm, const Real le) const {
    Real alpha = le, two(2);
    if ( coeff_ > ROL_EPSILON<Real>() ) {
      alpha = std::max(coeff_*std::pow(dm,two),le);
    }
    return alpha;
  }

  const Real alpha(const unsigned i) const {
    return computeAlpha(distanceMeasures_[i],linearizationErrors_[i]);
  }

  unsigned size(void) const {
    return size_;
  }

  void aggregate(Vector<Real> &aggSubGrad, Real &aggLinErr, Real &aggDistMeas) const {
    Real zero(0), one(1);
    aggSubGrad.zero(); aggLinErr = zero; aggDistMeas = zero; eG_->zero();
    Real eLE(0), eDM(0), yLE(0), yDM(0), tLE(0), tDM(0);
    for (unsigned i = 0; i < size_; i++) {
      // Compute aggregate subgradient using Kahan's compensated sum
      tG_->set(aggSubGrad);
      yG_->set(*subgradients_[i]); yG_->scale(dualVariables_[i]); yG_->plus(*eG_);
      aggSubGrad.set(*tG_); aggSubGrad.plus(*yG_);
      eG_->set(*tG_); eG_->axpy(-one,aggSubGrad); eG_->plus(*yG_);
      // Compute aggregate linearization error using Kahan's compensated sum
      tLE = aggLinErr;
      yLE = dualVariables_[i]*linearizationErrors_[i] + eLE;
      aggLinErr = tLE + yLE;
      eLE = (tLE - aggLinErr) + yLE;
      // Compute aggregate distance measure using Kahan's compensated sum
      tDM = aggDistMeas;
      yDM = dualVariables_[i]*distanceMeasures_[i] + eDM;
      aggDistMeas = tDM + yDM;
      eDM = (tDM - aggDistMeas) + yDM;
    }
  }

  void reset(const Vector<Real> &g, const Real le, const Real dm) {
    if (size_ == maxSize_) {
      // Find indices to remove
      unsigned loc = size_, cnt = 0;
      std::vector<unsigned> ind(remSize_,0);
      for (unsigned i = size_; i > 0; --i) {
        if ( std::abs(linearizationErrors_[i-1]) < ROL_EPSILON<Real>() ) {
          loc = i-1;
          break;
        }
      }
      for (unsigned i = 0; i < size_; i++) {
        if ( i < loc || i > loc ) {
          ind[cnt] = i;
          cnt++;
        }
        if (cnt == remSize_) {
          break;
        }
      }
      // Remove indices
      remove(ind);
      // Add aggregate subgradient
      add(g,le,dm);
    }
  }

  void update(const bool flag, const Real linErr, const Real distMeas,
              const Vector<Real> &g, const Vector<Real> &s) {
    Real zero(0);
    if ( flag ) {
      // Serious step taken: Update linearlization errors and distance measures
      for (unsigned i = 0; i < size_; i++) {
        linearizationErrors_[i] += linErr - subgradients_[i]->dot(s.dual());
        distanceMeasures_[i]    += distMeas;
      }
      linearizationErrors_[size_] = zero;
      distanceMeasures_[size_]    = zero;
    }
    else {
      // Null step taken: 
      linearizationErrors_[size_] = linErr;
      distanceMeasures_[size_]    = distMeas;
    }
    // Update (sub)gradient bundle
    (subgradients_[size_])->set(g);
    // Update dual variables
    dualVariables_[size_] = zero;
    // Update bundle size
    size_++;
  }

  // TT: adding access routines for derived classes
protected:
  
  Real getDualVariables (const unsigned i){
    return dualVariables_[i];
  }
  
  void setDualVariables(const unsigned i, const Real val) {
    dualVariables_[i] = val;
  }

  void resetDualVariables(void){
    Real zero(0);
    dualVariables_.assign(size_,zero);
  }

/***********************************************************************************************/
/***************** DUAL CUTTING PLANE PROBLEM ROUTINES *****************************************/
/***********************************************************************************************/
protected:
  Teuchos::RCP<Vector<Real> > gx_;
  Teuchos::RCP<Vector<Real> > ge_;

private:
  std::set<unsigned> workingSet_;
  std::set<unsigned> nworkingSet_;

  void initializeDualSolver(void) {
//    for (unsigned i = 0; i < maxSize_; i++) {
//      dualVariables_[i] = 0.0;
//    }
//    dualVariables_[0] = 1.0;
//    for (unsigned i = 0; i < maxSize_; i++) {
//      dualVariables_[i] = ((i<size_) ? 1.0/(Real)size_ : 0.0);
//    }
//    nworkingSet_.clear();
//    workingSet_.clear();
//    for (unsigned i = 0; i < size_; i++) {
//      nworkingSet_.insert(i);
//    }
    Real sum(0), err(0), tmp(0), y(0), zero(0);
    for (unsigned i = 0; i < size_; i++) {
      // Compute sum of dualVariables_ using Kahan's compensated sum
      tmp = sum;
      y   = dualVariables_[i] + err;
      sum = tmp + y;
      err = (tmp - sum) + y;
    }
    for (unsigned i = 0; i < size_; i++) {
      dualVariables_[i] /= sum;
    }
    nworkingSet_.clear();
    workingSet_.clear();
    for (unsigned i = 0; i < size_; i++) {
      if ( dualVariables_[i] == zero ) {
        workingSet_.insert(i);
      }
      else {
        nworkingSet_.insert(i);
      }
    }
  }

  Real evaluateObjective(std::vector<Real> &g, const std::vector<Real> &x, const Real t) const {
    Real one(1), half(0.5);
    gx_->zero(); eG_->zero();
    for (unsigned i = 0; i < size_; i++) {
      // Compute Gx using Kahan's compensated sum
      tG_->set(*gx_);
      yG_->set(*eG_); yG_->axpy(x[i],*(subgradients_[i]));
      gx_->set(*tG_); gx_->plus(*yG_);
      eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->plus(*yG_);
    }
    Real Hx(0), val(0), err(0), tmp(0), y(0);
    for (unsigned i = 0; i < size_; i++) {
      // Compute < g_i, Gx >
      Hx   = gx_->dot(*(subgradients_[i]));
      // Add to the objective function value using Kahan's compensated sum
      tmp  = val;
      y    = x[i]*(half*Hx + alpha(i)/t) + err;
      val  = tmp + y;
      err  = (tmp - val) + y;
      // Add gradient component
      g[i] = Hx + alpha(i)/t;
    }
    return val;
  }

  void applyFullMatrix(std::vector<Real> &Hx, const std::vector<Real> &x) const {
    Real one(1);
    gx_->zero(); eG_->zero();
    for (unsigned i = 0; i < size_; i++) {
      // Compute Gx using Kahan's compensated sum
      tG_->set(*gx_);
      yG_->set(*eG_); yG_->axpy(x[i],*(subgradients_[i]));
      gx_->set(*tG_); gx_->plus(*yG_);
      eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->plus(*yG_);
    }
    for (unsigned i = 0; i < size_; i++) {
      // Compute < g_i, Gx >
      Hx[i] = subgradients_[i]->dot(*gx_);
    }
  }

  void applyMatrix(std::vector<Real> &Hx, const std::vector<Real> &x) const {
    Real one(1);
    gx_->zero(); eG_->zero();
    unsigned n = nworkingSet_.size();
    typename std::set<unsigned>::iterator it = nworkingSet_.begin(); 
    for (unsigned i = 0; i < n; i++) {
      // Compute Gx using Kahan's compensated sum
      tG_->set(*gx_);
      yG_->set(*eG_); yG_->axpy(x[i],*(subgradients_[*it]));
      gx_->set(*tG_); gx_->plus(*yG_);
      eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->plus(*yG_);
      it++;
    }
    it = nworkingSet_.begin();
    for (unsigned i = 0; i < n; i++) {
      // Compute < g_i, Gx >
      Hx[i] = subgradients_[*it]->dot(*gx_); it++;
    }
  }

  void computeLagMult(std::vector<Real> &lam, const Real mu, const std::vector<Real> &g) const {
    Real zero(0);
    unsigned n = workingSet_.size();
    if ( n > 0 ) {
      lam.resize(n,zero);
      typename std::set<unsigned>::iterator it = workingSet_.begin();
      for (unsigned i = 0; i < n; i++) {
        lam[i] = g[*it] - mu; it++;
      }
    }
    else {
      lam.clear();
    }
  }
 
  bool isNonnegative(unsigned &ind, const std::vector<Real> &x) const {
    bool flag = true;
    unsigned n = workingSet_.size(); ind = size_;
    if ( n > 0 ) {
      Real min = ROL_OVERFLOW<Real>();
      typename std::set<unsigned>::iterator it = workingSet_.begin();
      for (unsigned i = 0; i < n; i++) {
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

  Real computeAlpha(unsigned &ind, const std::vector<Real> &x, const std::vector<Real> &p) const {
    Real alpha(1), tmp(0), zero(0); ind = size_;
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

  unsigned solveEQPsubproblem(std::vector<Real> &s, Real &mu,
                        const std::vector<Real> &g, const Real tol) const {
    // Build reduced QP information
    Real zero(0);
    unsigned n = nworkingSet_.size(), cnt = 0;
    mu = zero;
    s.assign(size_,zero);
    if ( n > 0 ) {
      std::vector<Real> gk(n,zero);
      typename std::set<unsigned>::iterator it = nworkingSet_.begin();
      for (unsigned i = 0; i < n; i++) {
        gk[i] = g[*it]; it++;
      }
      std::vector<Real> sk(n,zero);
      cnt = projectedCG(sk,mu,gk,tol);
      it  = nworkingSet_.begin();
      for (unsigned i = 0; i < n; i++) {
        s[*it] = sk[i]; it++;
      }
    }
    return cnt;
  }

  void applyPreconditioner(std::vector<Real> &Px, const std::vector<Real> &x) const {
    Real zero(0);
    int type = 0;
    std::vector<Real> tmp(Px.size(),zero);
    switch (type) {
      case 0: applyPreconditioner_Identity(tmp,x); break;
      case 1: applyPreconditioner_Jacobi(tmp,x);   break;
      case 2: applyPreconditioner_SymGS(tmp,x);    break;
    }
    applyPreconditioner_Identity(Px,tmp);
  }

  void applyG(std::vector<Real> &Gx, const std::vector<Real> &x) const {
    int type = 0;
    switch (type) {
      case 0: applyG_Identity(Gx,x); break;
      case 1: applyG_Jacobi(Gx,x);   break;
      case 2: applyG_SymGS(Gx,x);    break;
    }
  }

  void applyPreconditioner_Identity(std::vector<Real> &Px, const std::vector<Real> &x) const {
    unsigned dim = nworkingSet_.size();
    Real sum(0), err(0), tmp(0), y(0);
    for (unsigned i = 0; i < dim; i++) {
      // Compute sum of x using Kahan's compensated sum
      tmp = sum;
      y   = x[i] + err;
      sum = tmp + y;
      err = (tmp - sum) + y;
    }
    sum /= (Real)dim;
    for (unsigned i = 0; i < dim; i++) {
      Px[i] = x[i] - sum;
    }
  }

  void applyG_Identity(std::vector<Real> &Gx, const std::vector<Real> &x) const {
    Gx.assign(x.begin(),x.end());
  }

  void applyPreconditioner_Jacobi(std::vector<Real> &Px, const std::vector<Real> &x) const {
    unsigned dim = nworkingSet_.size();
    Real eHe(0), sum(0), one(1), zero(0);
    Real errX(0), tmpX(0), yX(0), errE(0), tmpE(0), yE(0);
    std::vector<Real> gg(dim,zero);
    typename std::set<unsigned>::iterator it = nworkingSet_.begin(); 
    for (unsigned i = 0; i < dim; i++) {
      gg[i] = one/std::abs(subgradients_[*it]->dot(*(subgradients_[*it]))); it++;
      // Compute sum of inv(D)x using Kahan's aggregated sum
      tmpX = sum;
      yX   = x[i]*gg[i] + errX;
      sum  = tmpX + errX;
      errX = (tmpX - sum) + yX;
      // Compute sum of inv(D)e using Kahan's aggregated sum
      tmpE = eHe;
      yE   = gg[i] + errE;
      eHe  = tmpE + yE;
      errE = (tmpE - eHe) + yE;
    }
    sum /= eHe;
    for (unsigned i = 0; i < dim; i++) {
      Px[i] = (x[i]-sum)*gg[i];
    }
  }

  void applyG_Jacobi(std::vector<Real> &Gx, const std::vector<Real> &x) const {
    unsigned dim = nworkingSet_.size();
    typename std::set<unsigned>::iterator it = nworkingSet_.begin();
    for (unsigned i = 0; i < dim; i++) {
      Gx[i] = std::abs(subgradients_[*it]->dot(*(subgradients_[*it])))*x[i]; it++;
    }
  }

  void applyPreconditioner_SymGS(std::vector<Real> &Px, const std::vector<Real> &x) const {
    int dim = nworkingSet_.size();
    //unsigned cnt = 0;
    gx_->zero(); ge_->zero();
    Real eHx(0), eHe(0), one(1);
    // Forward substitution
    std::vector<Real> x1(dim,0), e1(dim,0),gg(dim,0);
    typename std::set<unsigned>::iterator it, jt;
    it = nworkingSet_.begin(); 
    for (int i = 0; i < dim; i++) {
      gx_->zero(); ge_->zero(); jt = nworkingSet_.begin();
      for (int j = 0; j < i; j++) {
        gx_->axpy(x1[j],*(subgradients_[*jt]));
        ge_->axpy(e1[j],*(subgradients_[*jt]));
        jt++;
      }
      gg[i] = subgradients_[*it]->dot(*(subgradients_[*it]));
      x1[i] = (x[i] - gx_->dot(*(subgradients_[*it])))/gg[i];
      e1[i] = (one  - ge_->dot(*(subgradients_[*it])))/gg[i];
      it++;
    }
    // Apply diagonal
    for (int i = 0; i < dim; i++) {
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
        gx_->axpy(Hx[j],*(subgradients_[*jt]));
        ge_->axpy(He[j],*(subgradients_[*jt]));
      }
      Hx[i] = (x1[i] - gx_->dot(*(subgradients_[*it])))/gg[i];
      He[i] = (e1[i] - ge_->dot(*(subgradients_[*it])))/gg[i];
      // Compute sums
      eHx += Hx[i];
      eHe += He[i];
    }
    // Accumulate preconditioned vector
    for (int i = 0; i < dim; i++) {
      Px[i] = Hx[i] - (eHx/eHe)*He[i];
    }
  }

  void applyG_SymGS(std::vector<Real> &Gx, const std::vector<Real> &x) const {
    unsigned dim = nworkingSet_.size();
    typename std::set<unsigned>::iterator it = nworkingSet_.begin();
    for (unsigned i = 0; i < dim; i++) {
      Gx[i] = std::abs(subgradients_[*it]->dot(*(subgradients_[*it])))*x[i]; it++;
    }
  }

  void computeResidualUpdate(std::vector<Real> &r, std::vector<Real> &g) const {
    unsigned n = g.size();
    std::vector<Real> Gg(n,0);
    Real y(0), ytmp(0), yprt(0), yerr(0);
    applyPreconditioner(g,r);
    applyG(Gg,g);
    // Compute multiplier using Kahan's compensated sum
    for (unsigned i = 0; i < n; i++) {
      ytmp = y;
      yprt = (r[i] - Gg[i]) + yerr;
      y    = ytmp + yprt;
      yerr = (ytmp - y) + yprt;
    }
    y /= (Real)n;
    for (unsigned i = 0; i < n; i++) {
      r[i] -= y;
    }
    applyPreconditioner(g,r);
  }

  unsigned projectedCG(std::vector<Real> &x, Real &mu, const std::vector<Real> &b, const Real tol) const {
    Real one(1), zero(0);
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
    for (unsigned i = 0; i < n; i++) {
      tmp = mu;
      y   = r0[i] + err;
      mu  = tmp + y;
      err = (tmp - mu) + y;
    }
    mu /= (Real)n;
    // Return iteration count
    return cnt;
  }

  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) const {
    // Compute dot product of two vectors using Kahan's compensated sum
    Real val(0), err(0), tmp(0), y0(0);
    unsigned n = std::min(x.size(),y.size());
    for (unsigned i = 0; i < n; i++) {
      tmp = val;
      y0  = x[i]*y[i] + err;
      val = tmp + y0;
      err = (tmp - val) + y0;
    }
    return val;
  }

  Real norm(const std::vector<Real> &x) const {
    return std::sqrt(dot(x,x));
  }

  void axpy(const Real a, const std::vector<Real> &x, std::vector<Real> &y) const {
    unsigned n = std::min(y.size(),x.size());
    for (unsigned i = 0; i < n; i++) {
      y[i] += a*x[i];
    }
  }

  void scale(std::vector<Real> &x, const Real a) const {
    for (unsigned i = 0; i < x.size(); i++) {
      x[i] *= a;
    }
  }

  void scale(std::vector<Real> &x, const Real a, const std::vector<Real> &y) const {
    unsigned n = std::min(x.size(),y.size());
    for (unsigned i = 0; i < n; i++) {
      x[i] = a*y[i];
    }
  }

  unsigned solveDual_dim1(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    dualVariables_[0] = (Real)1;
    //std::cout << "dim = " << size_ << "  iter = " << 0 << "  CONVERGED!\n";
    return 0;
  }

  unsigned solveDual_dim2(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    Real diffg  = gx_->dot(*gx_), zero(0), one(1), half(0.5);
    gx_->set(*subgradients_[0]); gx_->axpy(-one,*subgradients_[1]);
    if ( std::abs(diffg) > ROL_EPSILON<Real>() ) {
      Real diffa  = (alpha(0)-alpha(1))/t;
      Real gdiffg = subgradients_[1]->dot(*gx_);
      dualVariables_[0] = std::min(one,std::max(zero,-(gdiffg+diffa)/diffg));
      dualVariables_[1] = one-dualVariables_[0];
    }
    else {
      if ( std::abs(alpha(0)-alpha(1)) > ROL_EPSILON<Real>() ) {
        if ( alpha(0) < alpha(1) ) {
          dualVariables_[0] = one; dualVariables_[1] = zero;
        }
        else if ( alpha(0) > alpha(1) ) {
          dualVariables_[0] = zero; dualVariables_[1] = one;
        }
      }
      else {
        dualVariables_[0] = half; dualVariables_[1] = half;
      }
    }
    //std::cout << "dim = " << size_ << "  iter = " << 0 << "  CONVERGED!\n";
    return 0;
  }

  unsigned solveDual_arbitrary(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    initializeDualSolver();
    bool nonneg = false;
    unsigned ind = 0, i = 0, CGiter = 0;
    Real snorm(0), alpha(0), mu(0), one(1), zero(0);
    std::vector<Real> s(size_,0), Hs(size_,0), g(size_,0), lam(size_+1,0);
    //Real val = evaluateObjective(g,dualVariables_,t);
    evaluateObjective(g,dualVariables_,t);
    for (i = 0; i < maxit; i++) {
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
          if ( ind < size_ ) {
            workingSet_.erase(ind);
            nworkingSet_.insert(ind);
          }
        }
      }
      else {
        alpha = computeAlpha(ind,dualVariables_,s);
        if ( alpha > zero ) {
          axpy(alpha,s,dualVariables_);
          applyFullMatrix(Hs,s);
          axpy(alpha,Hs,g);
        }
        if (ind < size_) {
          workingSet_.insert(ind);
          nworkingSet_.erase(ind);
        }
      }
      //std::cout << "iter = " << i << "  snorm = " << snorm << "  alpha = " << alpha << "\n";
    }
    //Real crit = computeCriticality(g);
    //std::cout << "Criticality Measure: " << crit << "\n";
    //std::cout << "dim = " << size_ << "  iter = " << i << "   CGiter = " << CGiter << "  CONVERGED!\n";
    return i;
  }

public:
  virtual unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    unsigned iter = 0;
    if (size_ == 1) {
      iter = solveDual_dim1(t,maxit,tol);
    }
    else if (size_ == 2) {
      iter = solveDual_dim2(t,maxit,tol);
    }
    else {
      iter = solveDual_arbitrary(t,maxit,tol);
    }
    return iter;
  }

private:
  /************************************************************************/
  /********************** PROJECTION ONTO FEASIBLE SET ********************/
  /************************************************************************/
  void project(std::vector<Real> &x, const std::vector<Real> &v) const {
    std::vector<Real> vsort(size_,0);
    vsort.assign(v.begin(),v.end());
    std::sort(vsort.begin(),vsort.end());
    Real sum(-1), lam(0), zero(0), one(1);
    for (int i = size_-1; i > 0; i--) {
      sum += vsort[i];
      if ( sum >= ((Real)(size_-i))*vsort[i-1] ) {
        lam = sum/(Real)(size_-i);
        break;
      }
    }
    if (lam == zero) {
      lam = (sum+vsort[0])/(Real)size_;
    }
    for (int i = 0; i < size_; i++) {
      x[i] = std::max(zero,v[i] - lam);
    }
  }

  Real computeCriticality(const std::vector<Real> &g) {
    Real zero(0), one(1);
    std::vector<Real> x(size_,0), Px(size_,0);
    axpy(one,dualVariables_,x);
    axpy(-one,g,x);
    project(Px,x);
    scale(x,zero);
    axpy(one,dualVariables_,x);
    axpy(-one,Px,x);
    return norm(x);
  }
}; // class Bundle 

} // namespace ROL

#endif

//  void aggregate(Vector<Real> &aggSubGrad, Real &aggLinErr, Real &aggDistMeas) const {
//    aggSubGrad.zero(); aggLinErr = 0.0; aggDistMeas = 0.0;
//    for (unsigned i = 0; i < size_; i++) {
//      //if ( dualVariables_[i] > ROL_EPSILON<Real>() ) {
//        aggSubGrad.axpy(dualVariables_[i],*(subgradients_[i]));
//        aggLinErr   += dualVariables_[i]*linearizationErrors_[i];
//        aggDistMeas += dualVariables_[i]*distanceMeasures_[i];
//      //}
//    }
//  }
