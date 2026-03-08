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

namespace ROL {
namespace OED {

template<typename Real>
class Factors : public ProfiledClass<Real,std::string> {
protected:
  const Ptr<Constraint<Real>>        model_;
  const Ptr<Vector<Real>>            theta_;
  const Ptr<Vector<Real>>            obs_;
  const Ptr<Vector<Real>>            obs0_;
  const Ptr<Vector<Real>>            c_;
  const Ptr<SampleGenerator<Real>>   sampler_;
  const bool                         storage_;
  const bool                         obs1d_;
  const bool                         ortho_;
  std::vector<Ptr<Vector<Real>>>     X_;
  Ptr<SampledVector<Real>>           g_storage_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

  void evaluateModel(Vector<Real> &g, const std::vector<Real> &param) const;
  void setFactors();
  void mgs2(const std::vector<Ptr<Vector<Real>>> &Y) const;

public:
  Factors(const Ptr<Constraint<Real>>      &model,
          const Ptr<Vector<Real>>          &theta,
          const Ptr<Vector<Real>>          &obs,
          const Ptr<SampleGenerator<Real>> &sampler,
          bool                              storage = true,
          const Ptr<Vector<Real>>          &c = nullPtr,
          bool                              ortho = false);
  Factors(const Ptr<Objective<Real>>       &model,
          const Ptr<Vector<Real>>          &theta,
          const Ptr<SampleGenerator<Real>> &sampler,
          bool                              storage = true,
          bool                              ortho = false);
  Factors() : ProfiledClass<Real,std::string>("OED::Factors"),
    model_(nullPtr), theta_(nullPtr), obs_(nullPtr),
    obs0_(nullPtr), c_(nullPtr), sampler_(nullPtr),
    storage_(false), obs1d_(false), ortho_(false) {}

  virtual void setPredictionVector(const Vector<Real> &c);
  virtual void getPredictionVector(Vector<Real> &c) const;

  // Create a vector in the parameter space
  virtual Ptr<Vector<Real>> createParameterVector(bool dual=false) const;

  // Create a vector in the observation space
  virtual Ptr<Vector<Real>> createObservationVector(bool dual=false) const;

  // Compute c^T F[k] x
  virtual Real apply(const Vector<Real> &x, int k) const;

  // Compute Fx = F[k] x
  virtual void apply(Vector<Real> &Fx, const Vector<Real> &x, int k) const;

  // Compute Mx = F[k]^T R F[k] x
  virtual void applyProduct(Vector<Real> &Mx, const Vector<Real> &x, int k) const;

  // Compute y^T F[k]^T R F[k] x
  virtual Real applyProduct2(const Vector<Real> &x, const Vector<Real> &y, int k) const;

  // Get F[k]^T c
  virtual const Ptr<const Vector<Real>> get(int k) const;

  // Compute F(param)^T c
  virtual void evaluate(Vector<Real> &F, const std::vector<Real> &param) const;

  // void sumAll(Real *in, Real *out, int size) const;

  virtual int numFactors() const;

  virtual int numMySamples() const;

  virtual std::vector<Real> getSample(int k) const;

}; // class Factors

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_Factors_Def.hpp"

#endif
