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

namespace ROL::OED {

template<typename Real>
Factors<Real>::Factors(const Ptr<Constraint<Real>>& model,
                       const Ptr<Vector<Real>>& theta,
                       const Ptr<Vector<Real>>& obs,
		       const Ptr<SampleGenerator<Real>>& sampler)
  : ProfiledClass<Real,std::string>("OED::Factors"),
    model_(model), theta_(theta), obs_(obs), sampler_(sampler),
    fdim_(theta_->dimension()), odim_(obs_->dimension()) {}

template<typename Real>
Factors<Real>::Factors(const Ptr<Objective<Real>>& model,
                       const Ptr<Vector<Real>>& theta,
		       const Ptr<SampleGenerator<Real>>& sampler)
  : Factors<Real>(makePtr<ConstraintFromObjective<Real>>(model),
      theta,makePtr<SingletonVector<Real>>(),sampler) {}

// Create a vector in the parameter space
template<typename Real>
Ptr<Vector<Real>> Factors<Real>::createParameterVector(bool dual) const {
  return dual ? theta_->dual().clone() : theta_->clone();
}

// Create a vector in the observation space
template<typename Real>
Ptr<Vector<Real>> Factors<Real>::createObservationVector(bool dual) const {
  if (odim_==1) return makePtr<SingletonVector<Real>>();
  return dual ? obs_->dual().clone() : obs_->clone();
}

// Compute F[k] x
template<typename Real>
void Factors<Real>::apply(Vector<Real>& Jx, const Vector<Real> &x, int k) const {
  apply(Jx,x,getSample(k));
  //startTimer("apply");
  //Real tol = std::sqrt(ROL_EPSILON<Real>());
  //model_->setParameter(getSample(k));
  //model_->applyJacobian(Jx,x,*theta_,tol);
  //stopTimer("apply");
}

template<typename Real>
void Factors<Real>::apply(Vector<Real>& Jx, const Vector<Real> &x, const std::vector<Real>& pt) const {
  startTimer("apply");
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  model_->setParameter(pt);
  model_->applyJacobian(Jx,x,*theta_,tol);
  stopTimer("apply");
}

// Compute F[k]* x
template<typename Real>
void Factors<Real>::applyAdjoint(Vector<Real>& Jx, const Vector<Real> &x, int k) const {
  applyAdjoint(Jx,x,getSample(k));
  //startTimer("applyAdjoint");
  //Real tol = std::sqrt(ROL_EPSILON<Real>());
  //model_->setParameter(getSample(k));
  //model_->applyAdjointJacobian(Jx,x,*theta_,tol);
  //stopTimer("applyAdjoint");
}

template<typename Real>
void Factors<Real>::applyAdjoint(Vector<Real>& Jx, const Vector<Real> &x, const std::vector<Real>& pt) const {
  startTimer("applyAdjoint");
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  model_->setParameter(pt);
  model_->applyAdjointJacobian(Jx,x,*theta_,tol);
  stopTimer("applyAdjoint");
}

template<typename Real>
int Factors<Real>::numFactors() const {
  return fdim_;
}

template<typename Real>
int Factors<Real>::numObservations() const {
  return odim_;
}

template<typename Real>
int Factors<Real>::numMySamples() const {
  return sampler_->numMySamples();
}

template<typename Real>
std::vector<Real> Factors<Real>::getSample(int k) const {
  return sampler_->getMyPoint(k);
}

template<typename Real>
const Ptr<const Vector<Real>> Factors<Real>::getTheta() const {
  return theta_;
}

} // End ROL::OED Namespace

#endif
