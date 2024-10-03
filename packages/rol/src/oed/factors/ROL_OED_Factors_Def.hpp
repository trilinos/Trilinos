// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_FACTORS_DEF_HPP
#define ROL_OED_FACTORS_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
void Factors<Real>::evaluateModel(Vector<Real> &g, const std::vector<Real> &param) const {
  startTimer("evaluateModel");
  bool isComputed = storage_ ? g_storage_->get(g,param) : false;
  if (!isComputed) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    model_->setParameter(param);
    model_->applyAdjointJacobian(g,*c_,*theta_,tol);
    if (storage_) g_storage_->set(g,param);
  }
  stopTimer("evaluateModel");
}

template<typename Real>
void Factors<Real>::setFactors() {
  startTimer("setFactors");
  const int nsamples = sampler_->numMySamples();
  g_storage_ = makePtr<SampledVector<Real>>();
  X_.clear(); X_.resize(nsamples,nullPtr);
  for (int i = 0; i < nsamples; ++i) {
    X_[i] = theta_->dual().clone();
    evaluateModel(*X_[i],sampler_->getMyPoint(i));
  }
  stopTimer("setFactors");
}

template<typename Real>
Factors<Real>::Factors(const Ptr<Constraint<Real>>      &model,
        const Ptr<Vector<Real>>          &theta,
        const Ptr<Vector<Real>>          &obs,
        const Ptr<SampleGenerator<Real>> &sampler,
        bool                              storage,
        const Ptr<Vector<Real>>          &c)
  : ProfiledClass<Real,std::string>("OED::Factors"),
    model_(model), theta_(theta), obs_(obs), obs0_(obs_->clone()),
    c_(c == nullPtr ? obs_->dual().clone() : c), sampler_(sampler),
    storage_(storage), obs1d_(obs_->dimension()==1) {
  if (c == nullPtr) c_->setScalar(static_cast<Real>(1));
  setFactors();
}

template<typename Real>
Factors<Real>::Factors(const Ptr<Objective<Real>>       &model,
        const Ptr<Vector<Real>>          &theta,
        const Ptr<SampleGenerator<Real>> &sampler,
        bool                              storage)
  : Factors<Real>(makePtr<ConstraintFromObjective<Real>>(model),
      theta,makePtr<SingletonVector<Real>>(),sampler,storage) {}

template<typename Real>
void Factors<Real>::setPredictionVector(const Vector<Real> &c) {
  if (!obs1d_) {
    c_->set(c);
    setFactors();
  }
}

template<typename Real>
void Factors<Real>::getPredictionVector(Vector<Real> &c) const {
  c.set(*c_);
}

// Create a vector in the parameter space
template<typename Real>
Ptr<Vector<Real>> Factors<Real>::createParameterVector(bool dual) const {
  return dual ? theta_->dual().clone() : theta_->clone();
}

// Create a vector in the observation space
template<typename Real>
Ptr<Vector<Real>> Factors<Real>::createObservationVector(bool dual) const {
  if (obs1d_) return makePtr<SingletonVector<Real>>();
  return dual ? c_->clone() : obs_->clone();
}

// Compute c^T F[k] x
template<typename Real>
Real Factors<Real>::apply(const Vector<Real> &x, int k) const {
  startTimer("apply");
  if (k < 0 || k >= sampler_->numMySamples())
    throw Exception::NotImplemented(">>> ROL::OED::Factors::apply : Index is out of bounds!");
  Real val = x.dot(X_[k]->dual());
  stopTimer("apply");
  return val;
}

// Compute Fx = F[k] x
template<typename Real>
void Factors<Real>::apply(Vector<Real> &Fx, const Vector<Real> &x, int k) const {
  startTimer("apply");
  if (k < 0 || k >= sampler_->numMySamples())
    throw Exception::NotImplemented(">>> ROL::OED::Factors::apply : Index is out of bounds!");
  if (obs1d_) {
    Fx.setScalar(apply(x,k));
  }
  else {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    Fx.zero();
    model_->setParameter(sampler_->getMyPoint(k));
    model_->applyJacobian(Fx,x,*theta_,tol);
  }
  stopTimer("apply");
}

// Compute Mx = F[k]^T R F[k] x
template<typename Real>
void Factors<Real>::applyProduct(Vector<Real> &Mx, const Vector<Real> &x, int k) const {
  startTimer("applyProduct");
  if (k < 0 || k >= sampler_->numMySamples())
    throw Exception::NotImplemented(">>> ROL::OED::Factors::applyProduct : Index is out of bounds!");
  if (obs1d_) {
    Mx.set(*get(k));
    Mx.scale(apply(x,k));
  }
  else {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    Mx.zero(); obs_->zero();
    model_->setParameter(sampler_->getMyPoint(k));
    model_->applyJacobian(*obs_,x,*theta_,tol);
    model_->applyAdjointJacobian(Mx,obs_->dual(),*theta_,tol);
  }
  stopTimer("applyProduct");
}

// Compute y^T F[k]^T R F[k] x
template<typename Real>
Real Factors<Real>::applyProduct2(const Vector<Real> &x, const Vector<Real> &y, int k) const {
  startTimer("applyProduct2");
  if (k < 0 || k >= sampler_->numMySamples())
    throw Exception::NotImplemented(">>> ROL::OED::Factors::applyProduct2 : Index is out of bounds!");
  Real val(0);
  if (obs1d_) {
    Real Fx = get(k)->dot(x.dual());
    Real Fy = get(k)->dot(y.dual());
    val = Fx * Fy;
  }
  else {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    obs_->zero(); obs0_->zero();
    model_->setParameter(sampler_->getMyPoint(k));
    model_->applyJacobian(*obs_,x,*theta_,tol);
    model_->applyJacobian(*obs0_,y,*theta_,tol);
    val = obs_->dot(*obs0_);
  }
  stopTimer("applyProduct2");
  return val;
}

// Get F[k]^T c
template<typename Real>
const Ptr<const Vector<Real>> Factors<Real>::get(int k) const {
  startTimer("get");
  if (k < 0 || k >= sampler_->numMySamples())
    throw Exception::NotImplemented(">>> ROL::OED::Factors::get : Index is out of bounds!");
  stopTimer("get");
  return X_[k];
}

// Compute F(param)^T c
template<typename Real>
void Factors<Real>::evaluate(Vector<Real> &F, const std::vector<Real> &param) const {
  startTimer("evaluate");
  evaluateModel(F,param);
  stopTimer("evaluate");
}

template<typename Real>
int Factors<Real>::numFactors() const {
  return theta_->dimension();
}

template<typename Real>
int Factors<Real>::numMySamples() const {
  return sampler_->numMySamples();
}

template<typename Real>
std::vector<Real> Factors<Real>::getSample(int k) const {
  return sampler_->getMyPoint(k);
}

} // End OED Namespace
} // End ROL Namespace

#endif
