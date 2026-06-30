// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTEDLAGRANGIANOBJECTIVE2_H
#define ROL_AUGMENTEDLAGRANGIANOBJECTIVE2_H

#include "ROL_AugmentedLagrangianPenalty.hpp"
#include "ROL_ScalarController.hpp"

namespace ROL {

template <class Real>
class AugmentedLagrangianObjective2 : public Objective<Real> {

private:

  // Required for Augmented Lagrangian definition
  const Ptr<Objective<Real>>  obj_;
  const std::vector<Ptr<AugmentedLagrangianPenalty<Real>>> pvec_;

  // Auxiliary storage
  Ptr<Vector<Real>> dualOptVector_;

  // Objective evaluations
  Ptr<ScalarController<Real,int>> fval_;
  Ptr<VectorController<Real,int>> gradient_;

  // Objective function
  Real fscale_;

  // Evaluation counters
  int nfval_;
  int ngval_;

  // User defined options
  bool scaleLagrangian_;
  // int HessianApprox_;

public:

  AugmentedLagrangianObjective2(const Ptr<Objective<Real>> &obj,
                                const std::vector<Ptr<AugmentedLagrangianPenalty<Real>>> &pen,
                                const Vector<Real> &dualOptVec,
                                ParameterList &parlist)
    : obj_(obj), pvec_(pen), fscale_(1), nfval_(0), ngval_(0) {

    fval_     = makePtr<ScalarController<Real,int>>();
    gradient_ = makePtr<VectorController<Real,int>>();

    dualOptVector_ = dualOptVec.clone();

    ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    scaleLagrangian_ = sublist.get("Use Scaled Augmented Lagrangian", false);
  }

  AugmentedLagrangianObjective2(const Ptr<Objective<Real>> &obj,
                                const std::vector<Ptr<AugmentedLagrangianPenalty<Real>>> &pen,
                                const Vector<Real> &dualOptVec,
                                const bool scaleLagrangian)
    : obj_(obj), pvec_(pen), fscale_(1), nfval_(0), ngval_(0),
      scaleLagrangian_(scaleLagrangian) {

    fval_     = makePtr<ScalarController<Real,int>>();
    gradient_ = makePtr<VectorController<Real,int>>();

    dualOptVector_ = dualOptVec.clone();
  }

  virtual void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    obj_->update(x,type,iter);
    for (unsigned i = 0; i < unsigned(getNumberConstraints()); ++i) pvec_[i]->update(x,type,iter);
    fval_->objectiveUpdate(type);
    gradient_->objectiveUpdate(type);
  }

  void setScaling(const Real fscale = 1.0) {
    fscale_ = fscale;
  }

  void setScaling(const Real cscale, const int k) {
    pvec_[k]->setScaling(cscale);
  }

  Real getScaling(const int k) {
    return pvec_[k]->getScaling();
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
    // Compute objective function value
    Real val = getObjectiveValue(x,tol);
    val *= fscale_;
    // Compute penalty term
    for (unsigned i = 0; i < unsigned(getNumberConstraints()); ++i)
      val += pvec_[i]->value(x,tol);
    return val;
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Compute objective function gradient
    g.set(*getObjectiveGradient(x,tol));
    g.scale(fscale_);
    for (unsigned i = 0; i < unsigned(getNumberConstraints()); ++i) {
      pvec_[i]->gradient(*dualOptVector_,x,tol);
      // std::cout << "Gradient penalty term " << i << ": " << std::scientific << std::setprecision(8) << dualOptVector_->norm() << std::endl;
      g.plus(*dualOptVector_);
    }
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Apply objective Hessian to a vector
    obj_->hessVec(hv,v,x,tol);
    hv.scale(fscale_);
    for (unsigned i = 0; i < unsigned(getNumberConstraints()); ++i) {
      pvec_[i]->hessVec(*dualOptVector_,v,x,tol);
      // std::cout << "HV penalty term " << i << ": " << std::scientific << std::setprecision(8) << dualOptVector_->norm() << std::endl;
      hv.plus(*dualOptVector_);
    }
  }

  Real dualNorm( const Vector<Real> &x, Real &tol, const int k ) {
    return pvec_[k]->dualNorm(x,tol);
  }

  Real dualResidual( const Vector<Real> &x, Real &tol, const int k ) {
    return pvec_[k]->dualResidual(x,tol);
  }

  Real feasibility( const Vector<Real> &x, Real &tol, const int k ) {
    return pvec_[k]->feasibility(x,tol);
  }

  void setPenaltyParameter( const Real penaltyParameter, const int k ) {
    pvec_[k]->setPenaltyParameter(penaltyParameter);
  }

  Real getPenaltyParameter( const int k ) {
    return pvec_[k]->getPenaltyParameter();
  }
  void setMultiplier( const Vector<Real> &multiplier, const int k ) {
    pvec_[k]->setMultiplier(multiplier);
  }

  void updateMultiplier( const Vector<Real> &x, Real &tol, const int k ) {
    pvec_[k]->updateMultiplier(x,tol);
  }

  // Return objective function value
  Real getObjectiveValue (const Vector<Real> &x, Real &tol ) {
    Real val(0);
    int key(0);
    bool isComputed = fval_->get(val,key);
    if ( !isComputed ) {
      val = obj_->value(x,tol); nfval_++;
      fval_->set(val,key);
    }
    return val;
  }

  // Compute objective function gradient
  const Ptr<const Vector<Real>> getObjectiveGradient( const Vector<Real> &x, Real &tol ) {
    int key(0);
    if (!gradient_->isComputed(key)) {
      if (gradient_->isNull(key)) gradient_->allocate(*dualOptVector_,key);
      obj_->gradient(*gradient_->set(key),x,tol); ngval_++;
    }
    return gradient_->get(key);
  }

  // Return constraint value
  const Ptr<const Vector<Real>> getConstraintVec( const Vector<Real> &x, Real &tol, const int k ) {
    return pvec_[k]->getConstraintVec(x,tol);
  }

  int getNumberConstraints(void) const {
    return pvec_.size();
  }

  // Return total number of constraint evaluations
  int getNumberConstraintEvaluations(void) const {
    int ncval = 0;
    for (unsigned i = 0; i < unsigned(getNumberConstraints()); ++i)
      ncval += pvec_[i]->getNumberConstraintEvaluations();
    return ncval;
  }

  // Return total number of objective evaluations
  int getNumberFunctionEvaluations(void) const {
    return nfval_;
  }

  // Return total number of gradient evaluations
  int getNumberGradientEvaluations(void) const {
    return ngval_;
  }

  // Reset counters
  void reset(void) {
    nfval_ = 0; ngval_ = 0;
    for (unsigned i = 0; i < unsigned(getNumberConstraints()); ++i) pvec_[i]->reset();
  }

}; // class AugmentedLagrangianObjective2

}  // namespace ROL

#endif
