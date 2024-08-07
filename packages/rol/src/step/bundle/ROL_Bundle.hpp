// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLE_H
#define ROL_BUNDLE_H

#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"

#include "ROL_Ptr.hpp"

#include <vector>
#include <set>

/** \class ROL::Bundle
    \brief Provides the interface for and implements a bundle.
*/

namespace ROL {

template<class Real>
class Bundle {
/***********************************************************************************************/
/***************** BUNDLE STORAGE **************************************************************/
/***********************************************************************************************/
private: 
  std::vector<ROL::Ptr<Vector<Real> > > subgradients_;
  std::vector<Real> linearizationErrors_;
  std::vector<Real> distanceMeasures_;

  std::vector<Real> dualVariables_;

  ROL::Ptr<Vector<Real> > tG_;
  ROL::Ptr<Vector<Real> > eG_;
  ROL::Ptr<Vector<Real> > yG_;
  ROL::Ptr<Vector<Real> > gx_;
  ROL::Ptr<Vector<Real> > ge_;

  unsigned size_;

  unsigned maxSize_;
  unsigned remSize_;
  Real coeff_;
  Real omega_;

  bool isInitialized_;

  void remove(const std::vector<unsigned> &ind) {
    Real zero(0);
    for (unsigned j = ind.back()+1; j < size_; ++j) {
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
      for (unsigned j = ind[i-1]+1; j < size_; ++j) {
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
  
/***********************************************************************************************/
/***************** BUNDLE MODIFICATION AND ACCESS ROUTINES *************************************/
/***********************************************************************************************/
public:
  virtual ~Bundle(void) {}

  Bundle(const unsigned maxSize = 10,
         const Real coeff = 0.0,
         const Real omega = 2.0,
         const unsigned remSize = 2) 
    : size_(0), maxSize_(maxSize), remSize_(remSize),
      coeff_(coeff), omega_(omega), isInitialized_(false) {
    Real zero(0);
    remSize_ = ((remSize_ < 2) ? 2 : ((remSize_ > maxSize_-1) ? maxSize_-1 : remSize_));
    coeff_ = std::max(static_cast<Real>(0),coeff_);
    omega_ = std::max(static_cast<Real>(1),omega_);
    subgradients_.clear();
    subgradients_.assign(maxSize,ROL::nullPtr);
    linearizationErrors_.clear();
    linearizationErrors_.assign(maxSize_,ROL_OVERFLOW<Real>());
    distanceMeasures_.clear();
    distanceMeasures_.assign(maxSize_,ROL_OVERFLOW<Real>());
    dualVariables_.clear();
    dualVariables_.assign(maxSize_,zero);
  }

  virtual void initialize(const Vector<Real> &g) {
    if ( !isInitialized_ ) {
      Real zero(0), one(1);
      for (unsigned i = 0; i < maxSize_; ++i) {
        subgradients_[i] = g.clone();
      }
      (subgradients_[0])->set(g);
      linearizationErrors_[0] = zero;
      distanceMeasures_[0]    = zero;
      dualVariables_[0]       = one;
      size_++;
      isInitialized_ = true;
      tG_ = g.clone();
      yG_ = g.clone();
      eG_ = g.clone();
      gx_ = g.clone();
      ge_ = g.clone();
    }
  }

  virtual unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) = 0;

  const Real linearizationError(const unsigned i) const {
    return linearizationErrors_[i];
  }

  const Real distanceMeasure(const unsigned i) const {
    return distanceMeasures_[i];
  }

  const Vector<Real> & subgradient(const unsigned i) const {
    return *(subgradients_[i]);
  }
  
  const Real getDualVariable(const unsigned i) const {
    return dualVariables_[i];
  }
  
  void setDualVariable(const unsigned i, const Real val) {
    dualVariables_[i] = val;
  }

  void resetDualVariables(void) {
    const Real zero(0);
    dualVariables_.assign(size_,zero);
  }

  const Real computeAlpha(const Real dm, const Real le) const {
    Real alpha = le;
    if ( coeff_ > ROL_EPSILON<Real>() ) {
      alpha = std::max(coeff_*std::pow(dm,omega_),le);
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
    for (unsigned i = 0; i < size_; ++i) {
      // Compute aggregate subgradient using Kahan's compensated sum
      //aggSubGrad.axpy(dualVariables_[i],*subgradients_[i]);
      yG_->set(*subgradients_[i]); yG_->scale(dualVariables_[i]); yG_->axpy(-one,*eG_);
      tG_->set(aggSubGrad); tG_->plus(*yG_);
      eG_->set(*tG_); eG_->axpy(-one,aggSubGrad); eG_->axpy(-one,*yG_);
      aggSubGrad.set(*tG_);
      // Compute aggregate linearization error using Kahan's compensated sum
      //aggLinErr += dualVariables_[i]*linearizationErrors_[i];
      yLE = dualVariables_[i]*linearizationErrors_[i] - eLE;
      tLE = aggLinErr + yLE;
      eLE = (tLE - aggLinErr) - yLE;
      aggLinErr = tLE;
      // Compute aggregate distance measure using Kahan's compensated sum
      //aggDistMeas += dualVariables_[i]*distanceMeasures_[i];
      yDM = dualVariables_[i]*distanceMeasures_[i] - eDM;
      tDM = aggDistMeas + yDM;
      eDM = (tDM - aggDistMeas) - yDM;
      aggDistMeas = tDM;
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
      for (unsigned i = 0; i < size_; ++i) {
        if ( i != loc ) {
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
      for (unsigned i = 0; i < size_; ++i) {
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

protected:
  const Real GiGj(const unsigned i, const unsigned j) const {
    return subgradient(i).dot(subgradient(j));
  }

  const Real dotGi(const unsigned i, const Vector<Real> &x) const {
    return x.dot(subgradient(i));
  }

  void addGi(const unsigned i, const Real a, Vector<Real> &x) const {
    x.axpy(a,subgradient(i));
  }

  Real evaluateObjective(std::vector<Real> &g, const std::vector<Real> &x, const Real t) const {
    Real one(1), half(0.5);
    gx_->zero(); eG_->zero();
    for (unsigned i = 0; i < Bundle<Real>::size(); ++i) {
      // Compute Gx using Kahan's compensated sum
      //gx_->axpy(x[i],*Bundle<Real>::subgradients_[i]);
      yG_->set(subgradient(i)); yG_->scale(x[i]); yG_->axpy(-one,*eG_);
      tG_->set(*gx_); tG_->plus(*yG_);
      eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->axpy(-one,*yG_);
      gx_->set(*tG_);
    }
    Real Hx(0), val(0), err(0), tmp(0), y(0);
    for (unsigned i = 0; i < size(); ++i) {
      // Compute < g_i, Gx >
      Hx   = dotGi(i,*gx_);
      // Add to the objective function value using Kahan's compensated sum
      //val += x[i]*(half*Hx + Bundle<Real>::alpha(i)/t);
      y    = x[i]*(half*Hx + alpha(i)/t) - err;
      tmp  = val + y;
      err  = (tmp - val) - y;
      val  = tmp;
      // Add gradient component
      g[i] = Hx + alpha(i)/t;
    }
    return val;
  }

  unsigned solveDual_dim1(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    setDualVariable(0,static_cast<Real>(1));
    //std::cout << "dim = " << Bundle<Real>::size() << "  iter = " << 0 << "  CONVERGED!\n";
    return 0;
  }

  unsigned solveDual_dim2(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) {
    Real diffg  = gx_->dot(*gx_), zero(0), one(1), half(0.5);
    gx_->set(subgradient(0)); addGi(1,-one,*gx_);
    if ( std::abs(diffg) > ROL_EPSILON<Real>() ) {
      Real diffa  = (alpha(0)-alpha(1))/t;
      Real gdiffg = dotGi(1,*gx_);
      setDualVariable(0,std::min(one,std::max(zero,-(gdiffg+diffa)/diffg)));
      setDualVariable(1,one-getDualVariable(0));
    }
    else {
      if ( std::abs(alpha(0)-alpha(1)) > ROL_EPSILON<Real>() ) {
        if ( alpha(0) < alpha(1) ) {
          setDualVariable(0,one); setDualVariable(1,zero);
        }
        else if ( alpha(0) > alpha(1) ) {
          setDualVariable(0,zero); setDualVariable(1,one);
        }
      }
      else {
        setDualVariable(0,half); setDualVariable(1,half);
      }
    }
    //std::cout << "dim = " << Bundle<Real>::size() << "  iter = " << 0 << "  CONVERGED!\n";
    return 0;
  }

}; // class Bundle 

} // namespace ROL

#endif

//  void aggregate(Vector<Real> &aggSubGrad, Real &aggLinErr, Real &aggDistMeas) const {
//    aggSubGrad.zero(); aggLinErr = 0.0; aggDistMeas = 0.0;
//    for (unsigned i = 0; i < size_; ++i) {
//      //if ( dualVariables_[i] > ROL_EPSILON<Real>() ) {
//        aggSubGrad.axpy(dualVariables_[i],*(subgradients_[i]));
//        aggLinErr   += dualVariables_[i]*linearizationErrors_[i];
//        aggDistMeas += dualVariables_[i]*distanceMeasures_[i];
//      //}
//    }
//  }
