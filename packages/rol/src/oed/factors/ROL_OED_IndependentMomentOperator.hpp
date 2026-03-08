// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_INDEPENDENT_COVARIANCE_OPERATOR_HPP
#define ROL_OED_INDEPENDENT_COVARIANCE_OPERATOR_HPP

#include "ROL_OED_MomentOperator.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class IndependentMomentOperator : public MomentOperator<Real> {
private:
  const std::vector<Ptr<MomentOperator<Real>>> Mvec_;
  const unsigned numOp_;

protected:
  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  IndependentMomentOperator(const std::vector<Ptr<MomentOperator<Real>>> &Mvec);

  virtual Ptr<MomentOperator<Real>> clone() const;
  void setMatrixNumber(int matNum);
  virtual void update(const Vector<Real> &p, UpdateType type, int iter = -1); 
  virtual void setFactors(const Ptr<Factors<Real>> &factors);
  virtual void generateFactors(const Ptr<Constraint<Real>>      &model,
                               const Ptr<Vector<Real>>          &theta,
                               const Ptr<Vector<Real>>          &obs,
                               const Ptr<SampleGenerator<Real>> &sampler,
                               bool                              storage = true,
                               const Ptr<Vector<Real>>          &c = nullPtr,
                               bool                              ortho = false);
  virtual void generateFactors(const Ptr<Objective<Real>>       &model,
                               const Ptr<Vector<Real>>          &theta,
                               const Ptr<SampleGenerator<Real>> &sampler,
                               bool                              storage = true,
                               bool                              ortho = false);
  virtual void setPerturbation(const Ptr<LinearOperator<Real>> &pOp);

  // Compute M(p)x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
  virtual void apply(Vector<Real> &Mx,
                     const Vector<Real> &x,
                     const Vector<Real> &p);

  // Compute M(q)x where M(q) = q_1 X_1 S_1 X_1' + ... + q_N X_N S_N X_N'
  // This function is distinct from apply to allow the user to store M(p)
  virtual void applyDeriv(Vector<Real> &Mx,
                          const Vector<Real> &x,
                          const Vector<Real> &q);

  // Compute inv(M(p))x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
  virtual void applyInverse(Vector<Real> &Mx,
                            const Vector<Real> &x,
                            const Vector<Real> &p);

  virtual void applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v);

  void setNoise(const Ptr<Noise<Real>> &noise, bool isHom = false);

  Real getNoise(int k) const;

  void getRegressionInfo(RegressionType &regType, bool &homNoise,
                         Ptr<Noise<Real>> &noise) const;

  virtual Real logDeterminant(const Vector<Real> &z);

}; // class IndependentMomentOperator

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_IndependentMomentOperator_Def.hpp"

#endif
