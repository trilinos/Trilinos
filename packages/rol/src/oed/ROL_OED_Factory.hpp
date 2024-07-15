// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_FACTORY_HPP
#define ROL_OED_FACTORY_HPP

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Types.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_ScaledObjective.hpp"
#include "ROL_ScalarLinearConstraint.hpp"
#include "ROL_LinearCombinationObjective.hpp"

#include "ROL_OED_BilinearConstraint.hpp"
#include "ROL_OED_ProbabilityConstraint.hpp"
#include "ROL_OED_LinearObjective.hpp"
#include "ROL_OED_A_HomObjective.hpp"
#include "ROL_OED_C_HomObjective.hpp"
#include "ROL_OED_D_HomObjective.hpp"
#include "ROL_OED_I_HomObjective.hpp"
#include "ROL_OED_Itrace_HomObjective.hpp"
#include "ROL_OED_A_HetObjective.hpp"
#include "ROL_OED_C_HetObjective.hpp"
#include "ROL_OED_D_HetObjective.hpp"
#include "ROL_OED_I_HetObjective.hpp"
#include "ROL_OED_Itrace_HetObjective.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_Radamacher.hpp"
#include "ROL_OED_L1Penalty.hpp"
#include "ROL_OED_DoubleWellPenalty.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class Factory {
private:
  // Input through constructor
  const Ptr<SampleGenerator<Real>> sampler_;
  const Ptr<Vector<Real>>          theta_;
  const Ptr<Vector<Real>>          obs_;

  // Parsed from parameter list
  bool useBudget_, useStorage_, useScale_, useL1_, useDWP_;
  Real objScale_, conScale_, L1penParam_, DWPparam_;

  Ptr<Vector<Real>> p_, w_;
  Ptr<Vector<Real>> ones_, zeros_;
  Ptr<Vector<Real>> cost_;
  Real budget_;

  Ptr<Factors<Real>>            factors_;
  Ptr<Factors<Real>>            factorsPV_;
  Ptr<TraceSampler<Real>>       traceSampler_;
  Ptr<MomentOperator<Real>>     cov0_;
  Ptr<MomentOperator<Real>>     cov1_;
  Ptr<BilinearConstraint<Real>> covar0_;
  Ptr<BilinearConstraint<Real>> covar1_;
  Ptr<LinearObjective<Real>>    lobj_;
  Ptr<QuadraticObjective<Real>> qobj_;

  Ptr<SampledVector<Real>> storage_;

  Ptr<Objective<Real>>       obj_;
  Ptr<Vector<Real>>          vec_;
  Ptr<BoundConstraint<Real>> bnd_;
  Ptr<Constraint<Real>>      econ_;
  Ptr<Vector<Real>>          emul_;
  Ptr<Constraint<Real>>      icon_;
  Ptr<Vector<Real>>          imul_;
  Ptr<BoundConstraint<Real>> ibnd_;

  std::string type_;
  bool isHom_;

public:
  Factory(const Ptr<Objective<Real>>           &model,
          const Ptr<SampleGenerator<Real>>     &sampler,
          const Ptr<Vector<Real>>              &theta,
          const Ptr<OED::MomentOperator<Real>> &cov,
                ParameterList                  &list);
  Factory(const Ptr<Constraint<Real>>          &model,
          const Ptr<SampleGenerator<Real>>     &sampler,
          const Ptr<Vector<Real>>              &theta,
          const Ptr<Vector<Real>>              &obs,
          const Ptr<OED::MomentOperator<Real>> &cov,
                ParameterList                  &list);

  void setBudgetConstraint(const Ptr<Vector<Real>> &cost, Real budget);

  Ptr<Problem<Real>> get(const Ptr<Vector<Real>> &c);
  Ptr<Problem<Real>> get(ParameterList &list,
      const Ptr<SampleGenerator<Real>> &sampler = nullPtr,
      const Ptr<Objective<Real>> &predFun = nullPtr);

  void check(std::ostream &stream = std::cout) const;

  void setDesign(Real val);
  void setDesign(const Vector<Real> &p);
  int loadDesign(const std::string &file, int dim, int n);
  const Ptr<const Vector<Real>> getDesign() const;
  const Ptr<Factors<Real>> getFactors() const;
  void printDesign(const std::string &name, const std::string &ext = ".txt") const;
  void printPredictionVariance(const Ptr<SampleGenerator<Real>> &sampler,
                               const std::string &name, const std::string &ext = ".txt") const;

  void profile(std::ostream &stream,
         const Ptr<BatchManager<Real>> &bman = nullPtr) const;
  void reset();

private:
  void buildHomObjective(ParameterList &list,
                         const Ptr<SampleGenerator<Real>> &sampler,
                         const Ptr<Objective<Real>> &predFun);
  void buildHomObjective(const Ptr<Vector<Real>> &c);
  void buildHetObjective(ParameterList &list,
                         const Ptr<SampleGenerator<Real>> &sampler,
                         const Ptr<Objective<Real>> &predFun);
  void buildHetObjective(const Ptr<Vector<Real>> &c);
  void buildVector();
  void buildBoundConstraint();
  void buildEqualityConstraint();
  void buildInequalityConstraint();

  void computeObjectiveScaling(const Ptr<Objective<Real>> &obj);

  void checkConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
                       std::ostream &stream) const;
  void checkObjective(const Ptr<Objective_SimOpt<Real>> &obj,
                      std::ostream &stream) const;
};

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_Factory_Def.hpp"

#endif
