// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_FACTORS_HPP
#define ROL_OED_FACTORS_HPP

#include "ROL_Ptr.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_SampledVector.hpp"
#include "ROL_SingletonVector.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_Types.hpp"
#include "ROL_OED_ProfiledClass.hpp"
#include <vector>

namespace ROL::OED {

template<typename Real>
class Factors : public ProfiledClass<Real,std::string> {
protected:
  const Ptr<Constraint<Real>>      model_;
  const Ptr<Vector<Real>>          theta_;
  const Ptr<Vector<Real>>          obs_;
  const Ptr<SampleGenerator<Real>> sampler_;
  const int                        fdim_;
  const int                        odim_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  Factors() : model_(nullPtr), theta_(nullPtr), obs_(nullPtr),
    sampler_(nullPtr), fdim_(0), odim_(0) {}

  Factors(const Ptr<Constraint<Real>>      &model,
          const Ptr<Vector<Real>>          &theta,
          const Ptr<Vector<Real>>          &obs,
          const Ptr<SampleGenerator<Real>> &sampler);

  Factors(const Ptr<Objective<Real>>       &model,
          const Ptr<Vector<Real>>          &theta,
          const Ptr<SampleGenerator<Real>> &sampler);

  // Create a vector in the parameter space
  virtual Ptr<Vector<Real>> createParameterVector(bool dual=false) const;

  // Create a vector in the observation space
  virtual Ptr<Vector<Real>> createObservationVector(bool dual=false) const;

  // Compute Fx = F[k] x
  virtual void apply(Vector<Real> &Fx, const Vector<Real> &x, int k) const;
  virtual void apply(Vector<Real> &Fx, const Vector<Real> &x, const std::vector<Real>& pt) const;

  // Compute Fx = F[k]* x
  virtual void applyAdjoint(Vector<Real> &Fx, const Vector<Real> &x, int k) const;
  virtual void applyAdjoint(Vector<Real> &Fx, const Vector<Real> &x, const std::vector<Real>& pt) const;

  virtual int numFactors() const;

  virtual int numObservations() const;

  virtual int numMySamples() const;

  virtual std::vector<Real> getSample(int k) const;

  virtual const Ptr<const Vector<Real>> getTheta() const;

}; // class Factors

} // End ROL::OED Namespace

#include "ROL_OED_Factors_Def.hpp"

#endif
