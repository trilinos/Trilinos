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

#include "ROL_OED_ObjectiveArray.hpp"
#include "ROL_OED_ProbabilityConstraint.hpp"
#include "ROL_OED_L1Penalty.hpp"
#include "ROL_OED_DoubleWellPenalty.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class Factory {
private:
  // Input through constructor
  const Ptr<ObjectiveArray<Real>>  objArray_;
  const Ptr<SampleGenerator<Real>> sampler_;
  const Ptr<Vector<Real>>          theta_;
  const Ptr<Vector<Real>>          obs_;

  // Parsed from parameter list
  const bool useStorage_, useScale_, useL1_, useDWP_;
  const Real L1penParam_, DWPparam_;
  const unsigned DWPtype_;
  Real objScale_, conScale_;
  bool useUpperBound_, useBudget_, equality_;

  Ptr<Vector<Real>> p_, w_;
  Ptr<Vector<Real>> cost_;
  Real budget_;

  Ptr<Vector<Real>>          vec_;
  Ptr<Objective<Real>>       obj_;
  Ptr<BoundConstraint<Real>> bnd_;
  Ptr<Constraint<Real>>      econ_;
  Ptr<Vector<Real>>          emul_;
  Ptr<Constraint<Real>>      icon_;
  Ptr<Vector<Real>>          imul_;
  Ptr<BoundConstraint<Real>> ibnd_;

  std::string type_;
  bool isHom_, usePF_, useDropout_;

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

  void setBudgetConstraint(const Ptr<Vector<Real>> &cost, Real budget, bool equality=false);
  void setProbabilityVector(const Ptr<Vector<Real>> &p);

  // Get robust C-optimality problem
  Ptr<Problem<Real>> get(const std::vector<Ptr<Vector<Real>>> &thetas,
                         const std::vector<Real> &weights,
                         const Ptr<Vector<Real>> &c,
                         ParameterList &list);
  // Get robust C-optimality problem
  Ptr<Problem<Real>> get(const std::vector<Ptr<Vector<Real>>> &thetas,
                         const Ptr<Vector<Real>> &c,
                         ParameterList &list);
  // Get C-optimality problem
  Ptr<Problem<Real>> get(const Ptr<Vector<Real>> &c);
  // Get general robust problem
  Ptr<Problem<Real>> get(const std::vector<Ptr<Vector<Real>>> &thetas,
                         const std::vector<Real> &weights,
                         ParameterList &list,
                         const Ptr<SampleGenerator<Real>> &sampler = nullPtr,
                         const Ptr<Objective<Real>> &predFun = nullPtr);
  // Get general robust problem
  Ptr<Problem<Real>> get(const std::vector<Ptr<Vector<Real>>> &thetas,
                         ParameterList &list,
                         const Ptr<SampleGenerator<Real>> &sampler = nullPtr,
                         const Ptr<Objective<Real>> &predFun = nullPtr);
  // Get general problem
  Ptr<Problem<Real>> get(ParameterList &list,
                         const Ptr<SampleGenerator<Real>> &sampler = nullPtr,
                         const Ptr<Objective<Real>> &predFun = nullPtr);

  //void check(std::ostream &stream = std::cout) const;

  void setDesign(Real val);
  void setDesign(const Vector<Real> &p);
  int loadDesign(const std::string &file, int dim, int n);
  const Ptr<const Vector<Real>> getDesign() const;
  const Ptr<const Vector<Real>> getOptimizationVector() const;
  const Ptr<Factors<Real>> getFactors(const Ptr<Vector<Real>>& theta) const;
  Ptr<Vector<Real>> createDesignVector() const;
  void printDesign(const std::string &name, const std::string &ext = ".txt") const;
  void printPredictionVariance(const Ptr<SampleGenerator<Real>> &sampler,
                               const std::string &name, const std::string &ext = ".txt") const;

  void profile(std::ostream &stream,
         const Ptr<BatchManager<Real>> &bman = nullPtr) const;
  void reset();

private:
  void buildVector(bool useMinMax = false);
  void buildBoundConstraint(bool useMinMax = false);
  void buildEqualityConstraint(bool useMinMax = false);
  void buildInequalityConstraint(bool useMinMax = false);
  void computeObjectiveScaling(const Ptr<Objective<Real>> &obj);
  //void checkConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
  //                     std::ostream &stream) const;
  //void checkObjective(const Ptr<Objective_SimOpt<Real>> &obj,
  //                    std::ostream &stream) const;
};

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_Factory_Def.hpp"

#endif
