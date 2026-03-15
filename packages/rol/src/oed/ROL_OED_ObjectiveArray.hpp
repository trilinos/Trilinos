// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_OBJECTIVEARRAY_HPP
#define ROL_OED_OBJECTIVEARRAY_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Objective.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_ProbabilityVector.hpp"
#include "ROL_DiagonalOperator.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class ObjectiveArray {
private:
  const Ptr<Objective<Real>>        modelO_;
  const Ptr<Constraint<Real>>       modelC_;
  Ptr<Vector<Real>>                  theta_;
  const Ptr<Vector<Real>>              obs_;
  const Ptr<MomentOperator<Real>>      cov_;
  const Ptr<SampleGenerator<Real>> sampler_;
  const bool                    useStorage_;
  unsigned                           nfact_;

  Ptr<const LinearOperator<Real>>     PFop_;
  Ptr<Vector<Real>>                   PFvc_;
  bool                               usePF_;

  std::vector<Ptr<Objective<Real>>> objVec_;
  std::vector<Real>                weights_;
  bool                               isHom_;
  std::string                         type_;

  Ptr<Objective<Real>> buildHomObjective(ParameterList& plist,
                                         const Ptr<Factors<Real>>& factors,
                                         const Ptr<MomentOperator<Real>>& cov0,
                                         const Ptr<Vector<Real>>& theta,
                                         const Ptr<SampleGenerator<Real>>& sampler = nullPtr,
                                         const Ptr<Objective<Real>>& predFun = nullPtr);
  Ptr<Objective<Real>> buildHomObjective(const Ptr<Vector<Real>>& c,
                                         const Ptr<Factors<Real>>& factors,
                                         const Ptr<MomentOperator<Real>>& cov0,
                                         const Ptr<Vector<Real>>& theta);
  Ptr<Objective<Real>> buildHetObjective(ParameterList& plist,
                                         const Ptr<Factors<Real>>& factors,
                                         const Ptr<MomentOperator<Real>>& cov0,
                                         const Ptr<MomentOperator<Real>>& cov1,
                                         const Ptr<Vector<Real>>& theta,
                                         const Ptr<SampleGenerator<Real>>& sampler = nullPtr,
                                         const Ptr<Objective<Real>>& predFun = nullPtr);
  Ptr<Objective<Real>> buildHetObjective(const Ptr<Vector<Real>>& c,
                                         const Ptr<Factors<Real>>& factors,
                                         const Ptr<MomentOperator<Real>>& cov0,
                                         const Ptr<MomentOperator<Real>>& cov1,
                                         const Ptr<Vector<Real>>& theta);

public:
  ObjectiveArray(const Ptr<Constraint<Real>>& model,
                 const Ptr<Vector<Real>>& obs,
                 const Ptr<SampleGenerator<Real>>& sampler,
                 const Ptr<MomentOperator<Real>>& cov,
                 ParameterList& plist);

  ObjectiveArray(const Ptr<Objective<Real>>& model,
                 const Ptr<SampleGenerator<Real>>& sampler,
                 const Ptr<MomentOperator<Real>>& cov,
                 ParameterList& plist);

  void setProbabilityVector(const Ptr<Vector<Real>>& p);

  void setTheta(const Ptr<Vector<Real>>& theta);

  void addObjective(const Ptr<Vector<Real>>& theta,
                    ParameterList& plist,
                    Real weight = Real(1),
                    const Ptr<SampleGenerator<Real>>& predSamp = nullPtr,
                    const Ptr<Objective<Real>>& predFunc = nullPtr);

  void addObjective(const Ptr<Vector<Real>>& theta,
                    const Ptr<Vector<Real>>& c,
                    Real weight = Real(1));

  const Ptr<Objective<Real>> getObjective(unsigned k) const;
  Real getWeight(unsigned k) const;
  std::vector<Real> getWeights() const;
  const Ptr<Factors<Real>> getFactors(const Ptr<Vector<Real>> &theta) const;
  const Ptr<MomentOperator<Real>> getBaseMomentOperator() const;

  unsigned numObjectives() const { return nfact_; }

  const Ptr<Vector<Real>> buildDesignVector() const {
    return makePtr<ProbabilityVector<Real>>(
           makePtr<std::vector<Real>>(sampler_->numMySamples(),0),
           sampler_->getBatchManager());
  }

}; // class ObjectiveArray

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_ObjectiveArray_Def.hpp"

#endif
