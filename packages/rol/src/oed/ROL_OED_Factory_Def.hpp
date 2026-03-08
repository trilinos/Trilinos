// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_FACTORY_DEF_HPP
#define ROL_OED_FACTORY_DEF_HPP

#include "ROL_OED_MeanValueObjective.hpp"
#include "ROL_OED_RobustObjective.hpp"
#include "ROL_OED_RobustConstraint.hpp"
#include "ROL_OED_AugmentedObjective.hpp"
#include "ROL_OED_AugmentedConstraint.hpp"

namespace ROL {
namespace OED {

template<typename Real>
Factory<Real>::Factory(const Ptr<Objective<Real>>           &model,
                       const Ptr<SampleGenerator<Real>>     &sampler,
                       const Ptr<Vector<Real>>              &theta,
                       const Ptr<OED::MomentOperator<Real>> &cov,
                             ParameterList                  &list)
  : objArray_(makePtr<ObjectiveArray<Real>>(model,sampler,cov,list)),
    sampler_(sampler), theta_(theta),
    useStorage_(list.sublist("OED").get("Use Storage",true)),
    useScale_(list.sublist("OED").get("Use Scaling",true)),
    useL1_(list.sublist("OED").get("Use L1 Penalty",false)),
    useDWP_(list.sublist("OED").get("Use Double-Well Penalty",false)),
    L1penParam_(list.sublist("OED").get("L1 Penalty Parameter",1.0)),
    DWPparam_(list.sublist("OED").get("Double-Well Penalty Parameter",1.0)),
    DWPtype_(list.sublist("OED").get("Double-Well Penalty Type",1)),
    objScale_(list.sublist("OED").get("Objective Scaling",-1.0)),
    conScale_(list.sublist("OED").get("Constraint Scaling",-1.0)),
    useUpperBound_(list.sublist("OED").get("Use Design Upper Bound", false)),
    useBudget_(false), equality_(false), usePF_(false), useDropout_(false) {
  if (useDWP_ && DWPtype_==1u) useUpperBound_ = true;
  RegressionType regType;
  Ptr<Noise<Real>> noise;
  cov->getRegressionInfo(regType,isHom_,noise);
}

template<typename Real>
Factory<Real>::Factory(const Ptr<Constraint<Real>>          &model,
                       const Ptr<SampleGenerator<Real>>     &sampler,
                       const Ptr<Vector<Real>>              &theta,
                       const Ptr<Vector<Real>>              &obs,
                       const Ptr<OED::MomentOperator<Real>> &cov,
                             ParameterList                  &list)
  : objArray_(makePtr<ObjectiveArray<Real>>(model,obs,sampler,cov,list)),
    sampler_(sampler), theta_(theta), obs_(obs),
    useStorage_(list.sublist("OED").get("Use Storage",true)),
    useScale_(list.sublist("OED").get("Use Scaling",true)),
    useL1_(list.sublist("OED").get("Use L1 Penalty",false)),
    useDWP_(list.sublist("OED").get("Use Double-Well Penalty",false)),
    L1penParam_(list.sublist("OED").get("L1 Penalty Parameter",1.0)),
    DWPparam_(list.sublist("OED").get("Double-Well Penalty Parameter",1.0)),
    DWPtype_(list.sublist("OED").get("Double-Well Penalty Type",1)),
    objScale_(list.sublist("OED").get("Objective Scaling",-1.0)),
    conScale_(list.sublist("OED").get("Constraint Scaling",-1.0)),
    useUpperBound_(list.sublist("OED").get("Use Design Upper Bound", false)),
    useBudget_(false), equality_(false), usePF_(false), useDropout_(false) {
  if (useDWP_ && DWPtype_==1u) useUpperBound_ = true;
  RegressionType regType;
  Ptr<Noise<Real>> noise;
  cov->getRegressionInfo(regType,isHom_,noise);
}

template<typename Real>
void Factory<Real>::setBudgetConstraint(const Ptr<Vector<Real>> &cost, Real budget, bool equality) {
  cost_      = cost;
  budget_    = budget;
  useBudget_ = true;
  equality_  = equality;
}

template<typename Real>
void Factory<Real>::setProbabilityVector(const Ptr<Vector<Real>> &p) {
  objArray_->setProbabilityVector(p);
  usePF_ = true;
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(const std::vector<Ptr<Vector<Real>>> &thetas,
                                      const std::vector<Real> &weights,
                                      const Ptr<Vector<Real>> &c,
                                      ParameterList &list) {
  bool useMinMax = list.sublist("OED").get("Use Minimax Formulation", false);
  useDropout_ = list.sublist("OED").get("Use Drop Out Sampling", false);

  if (useDropout_) objArray_->setTheta(theta_);
  for (unsigned i = 0u; i < thetas.size(); ++i)
    objArray_->addObjective(thetas[i],c,weights[i]);

  buildVector(useMinMax);
  buildBoundConstraint(useMinMax);
  buildEqualityConstraint(useMinMax);
  buildInequalityConstraint(useMinMax);

  // Build objective function
  std::vector<Ptr<Objective<Real>>> ovec;
  std::vector<Real> wvec;
  if (useL1_ || useDWP_) {
    if (useL1_) {
      ovec.push_back(makePtr<L1Penalty<Real>>());
      wvec.push_back(L1penParam_);
    }
    if (useDWP_) {
      ovec.push_back(makePtr<DoubleWellPenalty<Real>>(DWPtype_));
      wvec.push_back(DWPparam_);
    }
  }
  if (useMinMax) {
    if (useL1_ || useDWP_)
      obj_ = makePtr<AugmentedObjective<Real>>(ovec,wvec);
    else
      obj_ = makePtr<RobustObjective<Real>>();
  }
  else {
    //throw Exception::NotImplemented(">>> OED::Factory : Mean value robust OED not implemented!");
    auto obj0 = makePtr<MeanValueObjective<Real>>(objArray_);
    obj_ = obj0;
    if (useL1_ || useDWP_) {
      ovec.push_back(obj0);
      wvec.push_back(static_cast<Real>(1));
      obj_ = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
    }
  }

  // Build optimization problem
  Ptr<Problem<Real>> problem = makePtr<Problem<Real>>(obj_, vec_);

  // Add usual OED constraints
  problem->addBoundConstraint(bnd_);
  if (icon_!=nullPtr && imul_!=nullPtr && ibnd_!=nullPtr)
    problem->addLinearConstraint("Budget",icon_,imul_,ibnd_);
  if (econ_!=nullPtr && emul_!=nullPtr) {
    auto name = (useBudget_ && equality_) ? "Budget" : "Probability";
    problem->addLinearConstraint(name,econ_,emul_);
  }

  // Add MinMax constraint
  if (useMinMax) {
    std::vector<Real> ibnd_vec(thetas.size(),static_cast<Real>(0));
    auto                       icon = makePtr<RobustConstraint<Real>>(objArray_);
    Ptr<BoundConstraint<Real>> ibnd = makePtr<StdBoundConstraint<Real>>(ibnd_vec,false);
    Ptr<Vector<Real>>          imul = icon->buildRangeVector();
    problem->addConstraint("Robust Constraint",icon,imul,ibnd);
  }

  return problem;
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(const std::vector<Ptr<Vector<Real>>> &thetas,
                                      const Ptr<Vector<Real>> &c,
                                      ParameterList &list) {

  unsigned N = thetas.size();
  std::vector<Real> weights(N, static_cast<Real>(1) / static_cast<Real>(N));
  return Factory<Real>::get(thetas,weights,c,list);
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(const Ptr<Vector<Real>> &c) {
  objArray_->addObjective(theta_,c,static_cast<Real>(1));
  buildVector();
  buildBoundConstraint();
  buildEqualityConstraint();
  buildInequalityConstraint();

  Ptr<Objective<Real>> obj0, obj;
  obj0 = objArray_->getObjective(0u);
  computeObjectiveScaling(obj0);
  if (useScale_) obj = makePtr<ScaledObjective<Real>>(obj0,objScale_);
  else           obj = obj0;

  Ptr<Problem<Real>> problem;
  if (useL1_ || useDWP_) {
    std::vector<Ptr<Objective<Real>>> ovec;
    std::vector<Real> wvec;
    ovec.push_back(obj);
    wvec.push_back(static_cast<Real>(1));
    if (useL1_) {
      ovec.push_back(makePtr<L1Penalty<Real>>());
      wvec.push_back(L1penParam_);
    }
    if (useDWP_) {
      ovec.push_back(makePtr<DoubleWellPenalty<Real>>(DWPtype_));
      wvec.push_back(DWPparam_);
    }
    Ptr<Objective<Real>> penObj = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
    problem = makePtr<Problem<Real>>(penObj,vec_);
  }
  else {
    problem = makePtr<Problem<Real>>(obj,vec_);
  }
  problem->addBoundConstraint(bnd_);
  if (icon_!=nullPtr && imul_!=nullPtr && ibnd_!=nullPtr)
    problem->addLinearConstraint("Budget",icon_,imul_,ibnd_);
  if (econ_!=nullPtr && emul_!=nullPtr) {
    auto name = (useBudget_ && equality_) ? "Budget" : "Probability";
    problem->addLinearConstraint(name,econ_,emul_);
  }
  return problem;
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(const std::vector<Ptr<Vector<Real>>> &thetas,
                                      const std::vector<Real> &weights,
                                      ParameterList &list,
                                      const Ptr<SampleGenerator<Real>> &sampler,
                                      const Ptr<Objective<Real>> &predFun) {
  bool useMinMax = list.sublist("OED").get("Use Minimax Formulation", false);
  useDropout_ = list.sublist("OED").get("Use Drop Out Sampling", false);

  if (useDropout_) objArray_->setTheta(theta_);
  for (unsigned i = 0u; i < thetas.size(); ++i)
    objArray_->addObjective(thetas[i],list,weights[i],sampler,predFun);

  buildVector(useMinMax);
  buildBoundConstraint(useMinMax);
  buildEqualityConstraint(useMinMax);
  buildInequalityConstraint(useMinMax);

  // Build objective function
  std::vector<Ptr<Objective<Real>>> ovec;
  std::vector<Real> wvec;
  if (useL1_ || useDWP_) {
    if (useL1_) {
      ovec.push_back(makePtr<L1Penalty<Real>>());
      wvec.push_back(L1penParam_);
    }
    if (useDWP_) {
      ovec.push_back(makePtr<DoubleWellPenalty<Real>>(DWPtype_));
      wvec.push_back(DWPparam_);
    }
  }
  if (useMinMax) {
    if (useL1_ || useDWP_)
      obj_ = makePtr<AugmentedObjective<Real>>(ovec,wvec);
    else
      obj_ = makePtr<RobustObjective<Real>>();
  }
  else {
    //throw Exception::NotImplemented(">>> OED::Factory : Mean value robust OED not implemented!");
    auto obj0 = makePtr<MeanValueObjective<Real>>(objArray_);
    obj_ = obj0;
    if (useL1_ || useDWP_) {
      ovec.push_back(obj0);
      wvec.push_back(static_cast<Real>(1));
      obj_ = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
    }
  }

  // Build optimization problem
  Ptr<StochasticProblem<Real>> problem = makePtr<StochasticProblem<Real>>(obj_, vec_);

  // Add usual OED constraints
  problem->addBoundConstraint(bnd_);
  if (icon_!=nullPtr && imul_!=nullPtr && ibnd_!=nullPtr)
    problem->addLinearConstraint("Budget",icon_,imul_,ibnd_);
  if (econ_!=nullPtr && emul_!=nullPtr) {
    auto name = (useBudget_ && equality_) ? "Budget" : "Probability";
    problem->addLinearConstraint(name,econ_,emul_);
  }

  // Add MinMax constraint
  if (useMinMax) {
    std::vector<Real> ibnd_vec(thetas.size(),static_cast<Real>(0));
    auto                       icon = makePtr<RobustConstraint<Real>>(objArray_);
    Ptr<BoundConstraint<Real>> ibnd = makePtr<StdBoundConstraint<Real>>(ibnd_vec,false);
    Ptr<Vector<Real>>          imul = icon->buildRangeVector();
    problem->addConstraint("Robust Constraint",icon,imul,ibnd);
  }

  //std::string type = list.sublist("OED").get("Optimality Type","I");
  type_ = list.sublist("OED").get("Optimality Type","I");
  bool useTrace = list.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type_=="I" && !useTrace) || type_=="R") {
    if (sampler != nullPtr) {
      ParameterList slist;
      std::string sct = (type_=="I" ? "Risk Neutral" : "Risk Averse");
      slist.sublist("SOL").sublist("Objective").set("Type",                             sct);
      slist.sublist("SOL").sublist("Objective").set("Store Sampled Value and Gradient", useStorage_);
      bool usePD = false;
      if (type_=="R")
        usePD = list.sublist("OED").sublist("R-Optimality").get("Use Primal-Dual Algorithm",false);
      if (usePD && useMinMax)
        throw Exception::NotImplemented(">>> OED::Factory : Primal-dual risk is currently not supported for minimax problem!");
      if (!usePD || useMinMax) {
        if (type_=="R") {
          std::string rm = "Generalized Moreau-Yosida CVaR";
          Real alpha = list.sublist("OED").sublist("R-Optimality").get("Confidence Level",             0.9);
          Real lam   = list.sublist("OED").sublist("R-Optimality").get("Convex Combination Parameter", 0.0);
          Real eps   = list.sublist("OED").sublist("R-Optimality").get("Smoothing Parameter",          1e-4);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name",rm);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Confidence Level",             alpha);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Convex Combination Parameter", lam);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Smoothing Parameter",          eps);
        }
        if (useMinMax) problem->makeConstraintStochastic("Robust Constraint",slist,sampler);
        else //          problem->makeObjectiveStochastic(slist,sampler);
          throw Exception::NotImplemented(">>> OED::Factory : Trace-based type I and type R are currently not supported!");
      }
      else {
        if (type_=="R") {
          const Real one(1);
          std::string rm = "CVaR";
          Real alpha = list.sublist("OED").sublist("R-Optimality").get("Confidence Level",                0.9);
          Real lam   = list.sublist("OED").sublist("R-Optimality").get("Convex Combination Parameter",    0.0);
          list.sublist("SOL").sublist("Objective").set("Type", sct);
          list.sublist("SOL").sublist("Objective").sublist("Risk Measure").get("Name",rm);
          list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").get("Confidence Level",             alpha);
          list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").get("Convex Combination Parameter", one-lam);
        }
      }
    }
    else {
      throw Exception::NotImplemented(">>> OED::Factory : Sampler is a null pointer!");
    }
  }
  else if (type_=="A" || type_=="C" || type_=="D" || (type_=="I" && useTrace)) {
    // Do nothing
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }

  return problem;
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(const std::vector<Ptr<Vector<Real>>> &thetas,
                                      ParameterList &list,
                                      const Ptr<SampleGenerator<Real>> &sampler,
                                      const Ptr<Objective<Real>> &predFun) {
  unsigned N = thetas.size();
  std::vector<Real> weights(N, static_cast<Real>(1) / static_cast<Real>(N));
  return Factory<Real>::get(thetas,weights,list,sampler,predFun);
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(ParameterList &list,
                                const Ptr<SampleGenerator<Real>> &sampler,
                                const Ptr<Objective<Real>> &predFun) {
  objArray_->addObjective(theta_,list,static_cast<Real>(1),sampler,predFun);
  buildVector();
  buildBoundConstraint();
  buildEqualityConstraint();
  buildInequalityConstraint();

  Ptr<Objective<Real>> obj0, obj;
  obj0 = objArray_->getObjective(0u);
  computeObjectiveScaling(obj0);
  if (useScale_) obj = makePtr<ScaledObjective<Real>>(obj0,objScale_);
  else           obj = obj0;

  Ptr<StochasticProblem<Real>> problem;
  if (useL1_ || useDWP_) {
    std::vector<Ptr<Objective<Real>>> ovec;
    std::vector<Real> wvec;
    ovec.push_back(objArray_->getObjective(0u));
    wvec.push_back(static_cast<Real>(1));
    if (useL1_) {
      ovec.push_back(makePtr<L1Penalty<Real>>());
      wvec.push_back(L1penParam_);
    }
    if (useDWP_) {
      ovec.push_back(makePtr<DoubleWellPenalty<Real>>(DWPtype_));
      wvec.push_back(DWPparam_);
    }
    Ptr<Objective<Real>> penObj = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
    problem = makePtr<StochasticProblem<Real>>(penObj,vec_);
  }
  else {
    problem = makePtr<StochasticProblem<Real>>(objArray_->getObjective(0u),vec_);
  }
  problem->addBoundConstraint(bnd_);
  if (icon_!=nullPtr && imul_!=nullPtr && ibnd_!=nullPtr)
    problem->addLinearConstraint("Budget",icon_,imul_,ibnd_);
  if (econ_!=nullPtr && emul_!=nullPtr) {
    auto name = (useBudget_ && equality_) ? "Budget" : "Probability";
    problem->addLinearConstraint(name,econ_,emul_);
  }

  //std::string type = list.sublist("OED").get("Optimality Type","I");
  type_ = list.sublist("OED").get("Optimality Type","I");
  bool useTrace = list.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type_=="I" && !useTrace) || type_=="R") {
    if (sampler != nullPtr) {
      ParameterList slist;
      std::string sct = (type_=="I" ? "Risk Neutral" : "Risk Averse");
      slist.sublist("SOL").sublist("Objective").set("Type",                             sct);
      slist.sublist("SOL").sublist("Objective").set("Store Sampled Value and Gradient", useStorage_);
      bool usePD = false;
      if (type_=="R")
        usePD = list.sublist("OED").sublist("R-Optimality").get("Use Primal-Dual Algorithm",false);
      if (!usePD) {
        if (type_=="R") {
          std::string rm = "Generalized Moreau-Yosida CVaR";
          Real alpha = list.sublist("OED").sublist("R-Optimality").get("Confidence Level",             0.9);
          Real lam   = list.sublist("OED").sublist("R-Optimality").get("Convex Combination Parameter", 0.0);
          Real eps   = list.sublist("OED").sublist("R-Optimality").get("Smoothing Parameter",          1e-4);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name",rm);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Confidence Level",             alpha);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Convex Combination Parameter", lam);
          slist.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Smoothing Parameter",          eps);
        }
        problem->makeObjectiveStochastic(slist,sampler);
      }
      else {
        if (type_=="R") {
          const Real one(1);
          std::string rm = "CVaR";
          Real alpha = list.sublist("OED").sublist("R-Optimality").get("Confidence Level",                0.9);
          Real lam   = list.sublist("OED").sublist("R-Optimality").get("Convex Combination Parameter",    0.0);
          list.sublist("SOL").sublist("Objective").set("Type", sct);
          list.sublist("SOL").sublist("Objective").sublist("Risk Measure").get("Name",rm);
          list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").get("Confidence Level",             alpha);
          list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").get("Convex Combination Parameter", one-lam);
        }
      }
    }
    else {
      throw Exception::NotImplemented(">>> OED::Factory : Sampler is a null pointer!");
    }
  }
  else if (type_=="A" || type_=="C" || type_=="D" || (type_=="I" && useTrace)) {
    // Do nothing
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }
  return problem;
}

//template<typename Real>
//void Factory<Real>::check(std::ostream &stream) const {
//  if (covar0_ != nullPtr)            checkConstraint(covar0_,stream);
//  if (covar1_ != nullPtr && !isHom_) checkConstraint(covar1_,stream);
//  if (lobj_ != nullPtr && isHom_)    checkObjective(lobj_,stream);
//  if (qobj_ != nullPtr && !isHom_)   checkObjective(qobj_,stream);
//}

template<typename Real>
void Factory<Real>::setDesign(Real val) {
  p_->setScalar(val);
}

template<typename Real>
void Factory<Real>::setDesign(const Vector<Real> &p) {
  p_->set(p);
}

template<typename Real>
int Factory<Real>::loadDesign(const std::string &file, int dim, int n) {
  const Real tol(std::sqrt(ROL_EPSILON<Real>()));
  std::fstream input;
  input.open(file.c_str(), std::ios::in);
  if (!input.is_open()) {
    //throw Exception::NotImplemented(">>> OED::Factory:loadDesign : Cannot open input file!\n");
    return 1;
  }
  int gdim = sampler_->numGlobalSamples();
  int sdim = sampler_->getMyPoint(0).size();
  if (gdim != n || sdim != dim) {
    //throw Exception::NotImplemented(">>> OED::Factory::loadDesign : Dimension mismatch!\n");
    return 2;
  }
  int nsamp = sampler_->numMySamples();
  std::vector<Real> pin(nsamp), psamp(nsamp);
  Real diff(0), prob(0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < dim; ++j) {
      input >> pin[j];
    }
    input >> prob;
    for (int k = 0; k < nsamp; ++k) {
      psamp = sampler_->getMyPoint(k);
      diff  = static_cast<Real>(0);
      for (int j = 0; j < dim; ++j) {
        diff += std::pow(psamp[j]-pin[j],2);
      }
      diff = std::sqrt(diff);
      if (diff < tol) {
        dynamicPtrCast<ProbabilityVector<Real>>(p_)->setProbability(k,prob);
        break;
      }
    }
  }
  input.close();
  return 0;
}

template<typename Real>
const Ptr<const Vector<Real>> Factory<Real>::getDesign() const {
  return p_;
}

template<typename Real>
const Ptr<const Vector<Real>> Factory<Real>::getOptimizationVector() const {
  return vec_;
}

template<typename Real>
const Ptr<Factors<Real>> Factory<Real>::getFactors(const Ptr<Vector<Real>>& theta) const {
  return objArray_->getFactors(theta);
}

template<typename Real>
void Factory<Real>::printDesign(const std::string &name, const std::string &ext) const {
  int myRank = sampler_->batchID();
  int nsamp  = sampler_->numMySamples();
  int dim    = sampler_->getMyPoint(0).size();
  std::stringstream filename;
  filename << name << "_" << myRank << ext;
  std::ofstream file;
  file.open(filename.str());
  file << std::scientific << std::setprecision(15);
  for (int i = 0; i < nsamp; ++i) {
    for (int j = 0; j < dim; ++j) {
      file << std::right << std::setw(25)
           << sampler_->getMyPoint(i)[j];
    }
    file << std::right << std::setw(25)
         << dynamicPtrCast<ProbabilityVector<Real>>(p_)->getProbability(i)
         << std::endl;
  }
  file.close();
}

template<typename Real>
void Factory<Real>::printPredictionVariance(const Ptr<SampleGenerator<Real>> &sampler,
                             const std::string &name, const std::string &ext) const {
  Ptr<Objective<Real>> obj;
  auto factors = objArray_->getFactors(theta_);
  auto cov00 = objArray_->getBaseMomentOperator()->clone();
  cov00->setFactors(factors);
  cov00->setMatrixNumber(0);
  if (isHom_) {
    Ptr<BilinearConstraint<Real>> cov = makePtr<BilinearConstraint<Real>>(factors,cov00,"I");
    Ptr<LinearObjective<Real>> lobj = makePtr<LinearObjective<Real>>(factors,"I");
    obj = makePtr<Hom::I_Objective<Real>>(cov,lobj,theta_,useStorage_);
  }
  else {
    auto cov01 = objArray_->getBaseMomentOperator()->clone();
    cov01->setFactors(factors);
    cov01->setMatrixNumber(1);
    Ptr<BilinearConstraint<Real>> cov0 = makePtr<BilinearConstraint<Real>>(factors,cov00,"I");
    Ptr<BilinearConstraint<Real>> cov1 = makePtr<BilinearConstraint<Real>>(factors,cov01,"I");
    Ptr<QuadraticObjective<Real>> qobj = makePtr<QuadraticObjective<Real>>(cov0);
    obj = makePtr<Het::I_Objective<Real>>(cov1,qobj,theta_,useStorage_);
  }
  Real tol(1e-8);
  obj->update(*p_);
  int myRank = sampler->batchID();
  std::stringstream filename;
  filename << name << "_" << myRank << ext;
  std::ofstream file;
  file.open(filename.str());
  file << std::scientific << std::setprecision(15);
  for (int i = 0; i < sampler->numMySamples(); ++i) {
    std::vector<Real> param = sampler->getMyPoint(i);
    obj->setParameter(param);
    Real val = obj->value(*p_,tol);
    for (int j = 0; j < static_cast<int>(param.size()); ++j) {
      file << std::right << std::setw(25) << param[j];
    }
    file << std::right << std::setw(25) << val << std::endl;
  }
  file.close();
}

template<typename Real>
void Factory<Real>::profile(std::ostream &stream,
       const Ptr<BatchManager<Real>> &bman) const {
  //stream << std::endl << std::string(80,'-') << std::endl;
  //if (cov0_ != nullPtr) {
  //  cov0_->summarize(stream,bman);
  //}
  //if (covar0_ != nullPtr) {
  //  stream << std::string(80,'-') << std::endl;
  //  covar0_->summarize(stream,bman);
  //}
  //if (lobj_ != nullPtr) {
  //  stream << std::string(80,'-') << std::endl;
  //  lobj_->summarize(stream,bman);
  //}
  //if (factors_ != nullPtr) {
  //  stream << std::string(80,'-') << std::endl;
  //  factors_->summarize(stream,bman);
  //}
  //if (traceSampler_ != nullPtr && type_ == "A") {
  //  stream << std::string(80,'-') << std::endl;
  //  traceSampler_->summarize(stream,bman);
  //}
  Ptr<Objective<Real>> obj;
  for (unsigned i = 0u; i < objArray_->numObjectives(); ++i) {
    stream << std::string(80,'=') << std::endl;
    stream << "  Objective " << i << std::endl;
    stream << std::string(80,'-') << std::endl;
    
    if (usePF_ || useDropout_)
      obj = staticPtrCast<AffineTransformObjective<Real>>(objArray_->getObjective(i))->getObjective();
    else
      obj = objArray_->getObjective(i);
    staticPtrCast<BaseObjective<Real>>(obj)->summarize(stream,bman);
  }
  stream << std::string(80,'=') << std::endl << std::endl;
}

template<typename Real>
void Factory<Real>::reset() {
  for (unsigned i = 0u; i < objArray_->numObjectives(); ++i)
    staticPtrCast<BaseObjective<Real>>(objArray_->getObjective(i))->reset();
  //if (cov0_ != nullPtr)         cov0_->reset();
  //if (covar0_ != nullPtr)       covar0_->reset();
  //if (lobj_ != nullPtr)         lobj_->reset();
  //if (factors_ != nullPtr)      factors_->reset();
  //if (traceSampler_ != nullPtr) traceSampler_->reset();
}

template<typename Real>
Ptr<Vector<Real>> Factory<Real>::createDesignVector() const {
  return objArray_->buildDesignVector();
}

template<typename Real>
void Factory<Real>::buildVector(bool useMinMax) {
  p_ = Factory<Real>::createDesignVector();
  if (useBudget_) {
    p_->set(cost_->dual());
    p_->applyUnary(Elementwise::Reciprocal<Real>());
    p_->scale(budget_);
  }
  else p_->setScalar(static_cast<Real>(1)/static_cast<Real>(p_->dimension()));
  if (useMinMax) vec_ = makePtr<DesignVector<Real>>(p_);
  else           vec_ = p_;
}

template<typename Real>
void Factory<Real>::buildBoundConstraint(bool useMinMax) {
  //bnd_  = makePtr<Bounds<Real>>(zeros_,ones_);
  auto zeros = objArray_->buildDesignVector(); zeros->zero();
  auto ones  = objArray_->buildDesignVector();
  ones->setScalar(static_cast<Real>(1));
  Ptr<Vector<Real>> lbnd, ubnd;
  if (useMinMax) {
    lbnd = makePtr<DesignVector<Real>>(zeros,ROL_NINF<Real>());
    ubnd = makePtr<DesignVector<Real>>(ones, ROL_INF<Real>()); 
  }
  else {
    lbnd = zeros;
    ubnd = ones;
  }
  if (useUpperBound_)
    bnd_ = makePtr<ROL::Bounds<Real>>(lbnd,ubnd);
  else
    bnd_ = makePtr<ROL::Bounds<Real>>(*lbnd,true);
}

template<typename Real>
void Factory<Real>::buildEqualityConstraint(bool useMinMax) {
  if (!useBudget_) {
    //Real N(sampler_->numGlobalSamples());
    auto ones = objArray_->buildDesignVector();
    ones->setScalar(static_cast<Real>(1));
    if (useMinMax) {
      econ_ = makePtr<AugmentedConstraint<Real>>(*ones,useScale_,conScale_);
      emul_ = staticPtrCast<AugmentedConstraint<Real>>(econ_)->buildRangeVector()->dual().clone();
    }
    else {
      econ_ = makePtr<ProbabilityConstraint<Real>>(*ones,useScale_,conScale_);
      emul_ = makePtr<ProbabilityConstraintRangeVector<Real>>(sampler_->getBatchManager(),1);
      //emul_ = makePtr<DualScaledStdVector<Real>>(
      //          makePtr<std::vector<Real>>(1,0),
      //          makePtr<std::vector<Real>>(1,static_cast<Real>(1)));
      //          //makePtr<std::vector<Real>>(1,static_cast<Real>(1)/(N*N)));
      //          //makePtr<std::vector<Real>>(1,static_cast<Real>(1)/std::pow(N,2.0/3.0)));
    }
  }
  else {
    if (equality_) {
      if (useScale_) {
        const Real one(1);
        cost_->scale(one/budget_);
        budget_ = one;
      }
      if (useMinMax) {
        econ_ = makePtr<AugmentedConstraint<Real>>(cost_,budget_);
        emul_ = staticPtrCast<AugmentedConstraint<Real>>(icon_)->buildRangeVector()->dual().clone();
      }
      else {
        econ_ = makePtr<ScalarLinearConstraint<Real>>(cost_,budget_);
        emul_ = makePtr<SingletonVector<Real>>();
      }
    }
    else {
      econ_ = nullPtr;
      emul_ = nullPtr;
    }
  }
}

template<typename Real>
void Factory<Real>::buildInequalityConstraint(bool useMinMax) {
  if (!useBudget_ || equality_) {
    icon_ = nullPtr;
    imul_ = nullPtr;
    ibnd_ = nullPtr;
  }
  else {
    auto lbnd = makePtr<SingletonVector<Real>>(); lbnd->zero();
    if (useScale_) {
      const Real one(1);
      cost_->scale(one/budget_);
      budget_ = one;
    }
    if (useMinMax) {
      icon_ = makePtr<AugmentedConstraint<Real>>(cost_,budget_);
      imul_ = staticPtrCast<AugmentedConstraint<Real>>(icon_)->buildRangeVector()->dual().clone();
      ibnd_ = makePtr<Bounds<Real>>(*lbnd,false);
    }
    else {
      icon_ = makePtr<ScalarLinearConstraint<Real>>(cost_,budget_);
      imul_ = makePtr<SingletonVector<Real>>();
      ibnd_ = makePtr<Bounds<Real>>(*lbnd,false);
    }
  }
}

template<typename Real>
void Factory<Real>::computeObjectiveScaling(const Ptr<Objective<Real>> &obj) {
  if (useScale_ && objScale_ <= static_cast<Real>(0)) {
    const Real one(1);
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    Ptr<Vector<Real>> pu = p_->clone();
    Ptr<Vector<Real>> gu = p_->dual().clone();
    pu->setScalar(one/static_cast<Real>(p_->dimension()));
    obj->update(*pu,UpdateType::Temp);
    obj->gradient(*gu,*pu,tol);
    Real val = gu->norm();
    objScale_ = one/std::max(one,val);
  }
}

//template<typename Real>
//void Factory<Real>::checkConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
//                     std::ostream &stream) const {
//  stream << std::endl;
//  const Real zero(0), one(1);
//  Ptr<Vector<Real>> u, z, r, p, du, dz;
//  u  = theta_->clone();        u->randomize();
//  z  = p_->clone();            z->randomize(zero,one);
//  r  = theta_->dual().clone(); r->randomize();
//  p  = theta_->clone();        p->randomize();
//  du = theta_->clone();        du->randomize();
//  dz = p_->clone();            dz->randomize(zero,one);
//
//  stream << std::endl << "Check Jacobian_1 of Linear Constraint" << std::endl;
//  con->checkApplyJacobian_1(*u,*z,*du,*r,true,stream);
//
//  stream << std::endl << "Check Jacobian_2 of Linear Constraint" << std::endl;
//  con->checkApplyJacobian_2(*u,*z,*dz,*r,true,stream);
//
//  stream << std::endl << "Check Hessian_11 of Linear Constraint" << std::endl;
//  con->checkApplyAdjointHessian_11(*u,*z,*p,*du,*r,true,stream);
//
//  stream << std::endl << "Check Hessian_12 of Linear Constraint" << std::endl;
//  con->checkApplyAdjointHessian_12(*u,*z,*p,*du,*dz,true,stream);
//
//  stream << std::endl << "Check Hessian_21 of Linear Constraint" << std::endl;
//  con->checkApplyAdjointHessian_21(*u,*z,*p,*dz,*r,true,stream);
//
//  stream << std::endl << "Check Hessian_22 of Linear Constraint" << std::endl;
//  con->checkApplyAdjointHessian_22(*u,*z,*p,*dz,*dz,true,stream);
//
//  stream << std::endl << "Check Adjoint Jacobian_1 of Linear Constraint" << std::endl;
//  con->checkAdjointConsistencyJacobian_1(*p,*du,*u,*z,true,stream);
//
//  stream << std::endl << "Check Adjoint Jacobian_2 of Linear Constraint" << std::endl;
//  con->checkAdjointConsistencyJacobian_2(*p,*dz,*u,*z,true,stream);
//
//  stream << std::endl << "Check Linear Constraint Solve" << std::endl;
//  con->checkSolve(*u,*z,*r,true,stream);
//
//  stream << std::endl << "Check Inverse Jacobian_1 of Linear Constraint" << std::endl;
//  con->checkInverseJacobian_1(*r,*du,*u,*z,true,stream);
//
//  stream << std::endl << "Check Inverse Adjoint Jacobian_1 of Linear Constraint" << std::endl;
//  con->checkInverseAdjointJacobian_1(*r,*p,*u,*z,true,stream);
//}
//
//template<typename Real>
//void Factory<Real>::checkObjective(const Ptr<Objective_SimOpt<Real>> &obj,
//                    std::ostream &stream) const {
//  stream << std::endl;
//  Ptr<Vector<Real>> u, z, du, dz;
//  u  = theta_->clone(); u->randomize();
//  z  = p_->clone();     z->randomize();
//  du = theta_->clone(); du->randomize();
//  dz = p_->clone();     dz->randomize();
//
//  std::stringstream objType;
//  if (isHom_) objType << "Linear";
//  else        objType << "Quadratic";
//  stream << std::endl << "Check Gradient_1 of " << objType.str() << " Objective" << std::endl;
//  obj->checkGradient_1(*u,*z,*du,true,stream);
//
//  stream << std::endl << "Check Gradient_2 of " << objType.str() << " Objective" << std::endl;
//  obj->checkGradient_2(*u,*z,*dz,true,stream);
//
//  stream << std::endl << "Check Hessian_11 of " << objType.str() << " Objective" << std::endl;
//  obj->checkHessVec_11(*u,*z,*du,true,stream);
//
//  stream << std::endl << "Check Hessian_12 of " << objType.str() << " Objective" << std::endl;
//  obj->checkHessVec_12(*u,*z,*dz,true,stream);
//
//  stream << std::endl << "Check Hessian_21 of " << objType.str() << " Objective" << std::endl;
//  obj->checkHessVec_21(*u,*z,*du,true,stream);
//
//  stream << std::endl << "Check Hessian_22 of " << objType.str() << " Objective" << std::endl;
//  obj->checkHessVec_22(*u,*z,*dz,true,stream);
//}

} // End OED Namespace
} // End ROL Namespace

#endif
