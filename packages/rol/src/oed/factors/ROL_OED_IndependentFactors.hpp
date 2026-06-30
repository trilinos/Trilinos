// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_INDEPENDENT_FACTORS_HPP
#define ROL_OED_INDEPENDENT_FACTORS_HPP

#include "ROL_OED_Factors.hpp"
#include <vector>

namespace ROL::OED {

template<typename Real>
class IndependentFactors : public Factors<Real> {
protected:
  std::vector<Ptr<Factors<Real>>> factors_;
  unsigned numOp_, nobs_, nfact_;

  //using ProfiledClass<Real,std::string>::startTimer;
  //using ProfiledClass<Real,std::string>::stopTimer;

public:
  IndependentFactors(const Ptr<Objective<Real>>& model,
                     const Ptr<const Vector<Real>>& theta,
                     const Ptr<SampleGenerator<Real>>& sampler);
  IndependentFactors(const std::vector<Ptr<Factors<Real>>>& factors);

  // Create a vector in the parameter space
  Ptr<Vector<Real>> createParameterVector(bool dual=false) const override;

  // Create a vector in the observation space
  Ptr<Vector<Real>> createObservationVector(bool dual=false) const override;

  // Compute Fx = F[k] x
  void apply(Vector<Real> &Fx, const Vector<Real> &x, int k) const override;
  void apply(Vector<Real> &Fx, const Vector<Real> &x, const std::vector<Real>& pt) const override;

  // Compute Fx = F[k]* x
  void applyAdjoint(Vector<Real> &Fx, const Vector<Real> &x, int k) const override;
  void applyAdjoint(Vector<Real> &Fx, const Vector<Real> &x, const std::vector<Real>& pt) const override;

  int numFactors() const override { return nfact_; }

  int numObservations() const override { return nobs_; }

  int numMySamples() const override { return factors_[0]->numMySamples(); }

  std::vector<Real> getSample(int k) const override { return factors_[0]->getSample(k); }

  const Ptr<const Vector<Real>> getTheta() const override { return factors_[0]->getTheta(); }

  const Ptr<Factors<Real>> getFactors(unsigned i) const { return factors_[i]; }
}; // class Factors

} // End ROL::OED Namespace

#include "ROL_OED_IndependentFactors_Def.hpp"

#endif
