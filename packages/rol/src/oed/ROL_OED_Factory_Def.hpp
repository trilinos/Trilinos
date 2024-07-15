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

namespace ROL {
namespace OED {

template<typename Real>
Factory<Real>::Factory(const Ptr<Objective<Real>>           &model,
                       const Ptr<SampleGenerator<Real>>     &sampler,
                       const Ptr<Vector<Real>>              &theta,
                       const Ptr<OED::MomentOperator<Real>> &cov,
                             ParameterList                  &list)
  : sampler_(sampler), theta_(theta), useBudget_(false), useL1_(false),
    useDWP_(false), L1penParam_(1), DWPparam_(1), cov0_(cov) {
  const Real zero(0), one(1);
  useStorage_ = list.sublist("OED").get("Use Storage",true);
  useScale_   = list.sublist("OED").get("Use Scaling",true);
  objScale_   = list.sublist("OED").get("Objective Scaling",-1.0);
  conScale_   = list.sublist("OED").get("Constraint Scaling",-1.0);
  useL1_      = list.sublist("OED").get("Use L1 Penalty",false);
  useDWP_     = list.sublist("OED").get("Use Double-Well Penalty",false);
  L1penParam_ = list.sublist("OED").get("L1 Penalty Parameter",1.0);
  DWPparam_   = list.sublist("OED").get("Double-Well Penalty Parameter",1.0);

  factors_ = makePtr<Factors<Real>>(model,theta_,sampler_);
  cov0_->setFactors(factors_);
  cov0_->setMatrixNumber(0);

  RegressionType regType;
  Ptr<Noise<Real>> noise;
  cov0_->getRegressionInfo(regType,isHom_,noise);
  if (!isHom_) {
    cov1_ = cov0_->clone();
    cov1_->setFactors(factors_);
    cov1_->setMatrixNumber(1);
  }

  zeros_ = makePtr<ProbabilityVector<Real>>(
           makePtr<std::vector<Real>>(sampler->numMySamples(),zero),
           sampler_->getBatchManager());
  ones_  = makePtr<ProbabilityVector<Real>>(
           makePtr<std::vector<Real>>(sampler->numMySamples(),one),
           sampler_->getBatchManager());

  buildVector();
  buildBoundConstraint();
}

template<typename Real>
Factory<Real>::Factory(const Ptr<Constraint<Real>>          &model,
                       const Ptr<SampleGenerator<Real>>     &sampler,
                       const Ptr<Vector<Real>>              &theta,
                       const Ptr<Vector<Real>>              &obs,
                       const Ptr<OED::MomentOperator<Real>> &cov,
                             ParameterList                  &list)
  : sampler_(sampler), theta_(theta), obs_(obs), useBudget_(false), useL1_(false),
    useDWP_(false), L1penParam_(1), DWPparam_(1), cov0_(cov) {
  const Real zero(0), one(1);
  useStorage_ = list.sublist("OED").get("Use Storage",true);
  useScale_   = list.sublist("OED").get("Use Scaling",true);
  objScale_   = list.sublist("OED").get("Objective Scaling",-1.0);
  conScale_   = list.sublist("OED").get("Constraint Scaling",-1.0);
  useL1_      = list.sublist("OED").get("Use L1 Penalty",false);
  useDWP_     = list.sublist("OED").get("Use Double-Well Penalty",false);
  L1penParam_ = list.sublist("OED").get("L1 Penalty Parameter",1.0);
  DWPparam_   = list.sublist("OED").get("Double-Well Penalty Parameter",1.0);

  factors_ = makePtr<Factors<Real>>(model,theta_,obs_,sampler_);
  cov0_->setFactors(factors_);
  cov0_->setMatrixNumber(0);

  RegressionType regType;
  Ptr<Noise<Real>> noise;
  cov0_->getRegressionInfo(regType,isHom_,noise);
  if (!isHom_) {
    cov1_ = cov0_->clone();
    cov1_->setFactors(factors_);
    cov1_->setMatrixNumber(1);
  }

  zeros_ = makePtr<ProbabilityVector<Real>>(
           makePtr<std::vector<Real>>(sampler->numMySamples(),zero),
           sampler_->getBatchManager());
  ones_  = makePtr<ProbabilityVector<Real>>(
           makePtr<std::vector<Real>>(sampler->numMySamples(),one),
           sampler_->getBatchManager());

  buildVector();
  buildBoundConstraint();
}

template<typename Real>
void Factory<Real>::setBudgetConstraint(const Ptr<Vector<Real>> &cost, Real budget) {
  cost_ = cost;
  budget_ = budget;
  useBudget_ = true;
  vec_->setScalar(budget/static_cast<Real>(cost->dimension()));
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(const Ptr<Vector<Real>> &c) {
  buildEqualityConstraint();
  buildInequalityConstraint();

  Ptr<Problem<Real>> problem;
  if (isHom_) buildHomObjective(c);
  else        buildHetObjective(c);
  if (useL1_ || useDWP_) {
    std::vector<Ptr<Objective<Real>>> ovec;
    std::vector<Real> wvec;
    ovec.push_back(obj_);
    wvec.push_back(static_cast<Real>(1));
    if (useL1_) {
      ovec.push_back(makePtr<L1Penalty<Real>>());
      wvec.push_back(L1penParam_);
    }
    if (useDWP_) {
      ovec.push_back(makePtr<DoubleWellPenalty<Real>>());
      wvec.push_back(DWPparam_);
    }
    Ptr<Objective<Real>> penObj = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
    problem = makePtr<Problem<Real>>(penObj,vec_);
  }
  else {
    problem = makePtr<Problem<Real>>(obj_,vec_);
  }
  problem->addBoundConstraint(bnd_);
  if (useBudget_) problem->addLinearConstraint("Budget",icon_,imul_,ibnd_);
  else            problem->addLinearConstraint("Probability",econ_,emul_);
  return problem;
}

template<typename Real>
Ptr<Problem<Real>> Factory<Real>::get(ParameterList &list,
    const Ptr<SampleGenerator<Real>> &sampler,
    const Ptr<Objective<Real>> &predFun) {
  buildEqualityConstraint();
  buildInequalityConstraint();

  if (isHom_) buildHomObjective(list,sampler,predFun);
  else        buildHetObjective(list,sampler,predFun);
  Ptr<StochasticProblem<Real>> problem;
  if (useL1_ || useDWP_) {
    std::vector<Ptr<Objective<Real>>> ovec;
    std::vector<Real> wvec;
    ovec.push_back(obj_);
    wvec.push_back(static_cast<Real>(1));
    if (useL1_) {
      ovec.push_back(makePtr<L1Penalty<Real>>());
      wvec.push_back(L1penParam_);
    }
    if (useDWP_) {
      ovec.push_back(makePtr<DoubleWellPenalty<Real>>());
      wvec.push_back(DWPparam_);
    }
    Ptr<Objective<Real>> penObj = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
    problem = makePtr<StochasticProblem<Real>>(penObj,vec_);
  }
  else {
    problem = makePtr<StochasticProblem<Real>>(obj_,vec_);
  }
  problem->addBoundConstraint(bnd_);
  if (useBudget_) problem->addLinearConstraint("Budget",icon_,imul_,ibnd_);
  else            problem->addLinearConstraint("Probability",econ_,emul_);

  std::string type = list.sublist("OED").get("Optimality Type","I");
  bool useTrace = list.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type=="I" && !useTrace) || type_=="R") {
    if (sampler != nullPtr) {
      ParameterList slist;
      std::string sct = (type_=="I" ? "Risk Neutral" : "Risk Averse");
      slist.sublist("SOL").sublist("Objective").set("Type",                             sct);
      slist.sublist("SOL").sublist("Objective").set("Store Sampled Value and Gradient", useStorage_);
      bool usePD = false;
      if (type_=="R")
        usePD = list.sublist("OED").sublist("R-Optimality").get("Use Primal-Dual Algorithm",false);
      if (!usePD) {
        if (type=="R") {
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
        if (type=="R") {
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
  else if (type=="A" || type=="C" || type=="D" || (type=="I" && useTrace)) {
    // Do nothing
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }
  return problem;
}

template<typename Real>
void Factory<Real>::check(std::ostream &stream) const {
  if (covar0_ != nullPtr)            checkConstraint(covar0_,stream);
  if (covar1_ != nullPtr && !isHom_) checkConstraint(covar1_,stream);
  if (lobj_ != nullPtr && isHom_)    checkObjective(lobj_,stream);
  if (qobj_ != nullPtr && !isHom_)   checkObjective(qobj_,stream);
}

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
const Ptr<Factors<Real>> Factory<Real>::getFactors() const {
  return factors_;
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
  if (isHom_) {
    Ptr<BilinearConstraint<Real>> cov = makePtr<BilinearConstraint<Real>>(factors_,cov0_,"I");
    Ptr<LinearObjective<Real>> lobj = makePtr<LinearObjective<Real>>(factors_,"I");
    obj = makePtr<Hom::I_Objective<Real>>(cov,lobj,theta_,useStorage_);
  }
  else {
    Ptr<BilinearConstraint<Real>> cov0 = makePtr<BilinearConstraint<Real>>(factors_,cov0_,"I");
    Ptr<BilinearConstraint<Real>> cov1 = makePtr<BilinearConstraint<Real>>(factors_,cov1_,"I");
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
  stream << std::endl << std::string(80,'-') << std::endl;
  if (cov0_ != nullPtr) {
    cov0_->summarize(stream,bman);
  }
  if (covar0_ != nullPtr) {
    stream << std::string(80,'-') << std::endl;
    covar0_->summarize(stream,bman);
  }
  if (lobj_ != nullPtr) {
    stream << std::string(80,'-') << std::endl;
    lobj_->summarize(stream,bman);
  }
  if (factors_ != nullPtr) {
    stream << std::string(80,'-') << std::endl;
    factors_->summarize(stream,bman);
  }
  if (traceSampler_ != nullPtr && type_ == "A") {
    stream << std::string(80,'-') << std::endl;
    traceSampler_->summarize(stream,bman);
  }
  stream << std::string(80,'-') << std::endl << std::endl;
}

template<typename Real>
void Factory<Real>::reset() {
  if (cov0_ != nullPtr)         cov0_->reset();
  if (covar0_ != nullPtr)       covar0_->reset();
  if (lobj_ != nullPtr)         lobj_->reset();
  if (factors_ != nullPtr)      factors_->reset();
  if (traceSampler_ != nullPtr) traceSampler_->reset();
}

template<typename Real>
void Factory<Real>::buildHomObjective(ParameterList &list,
                       const Ptr<SampleGenerator<Real>> &sampler,
                       const Ptr<Objective<Real>> &predFun) {
  type_ = list.sublist("OED").get("Optimality Type","C");
  Ptr<Objective<Real>> obj;
  storage_ = makePtr<SampledVector<Real>>();
  bool useTrace = list.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type_ == "I" && !useTrace) || type_ == "R") {
    if (predFun == nullPtr) {
      covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_);
      lobj_   = makePtr<LinearObjective<Real>>(factors_,type_);
    }
    else {
      factorsPV_ = makePtr<Factors<Real>>(predFun,theta_,sampler);
      covar0_ = makePtr<BilinearConstraint<Real>>(factorsPV_,cov0_,type_);
      lobj_   = makePtr<LinearObjective<Real>>(factorsPV_,type_);
    }
    obj     = makePtr<Hom::I_Objective<Real>>(covar0_,lobj_,theta_,useStorage_);
    int dim = sampler_->getMyPoint(0).size();
    Real wt = 0.0;
    std::vector<Real> sum(dim,0.0), mean(dim,0.0), pt(dim,0.0);
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      pt = sampler_->getMyPoint(i);
      wt = sampler_->getMyWeight(i);
      for (int j = 0; j < dim; ++j) sum[j] += wt*pt[j];
    }
    sampler_->sumAll(&sum[0],&mean[0],dim);
    obj->setParameter(mean);
  }
  else if (type_ == "C") {
    Ptr<Vector<Real>> c = theta_->dual().clone();
    Real cval = list.sublist("OED").sublist("C-Optimality").get("C Value",1.0);
    c->setScalar(cval);
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,c);
    lobj_   = makePtr<LinearObjective<Real>>(c);
    obj     = makePtr<Hom::C_Objective<Real>>(covar0_,lobj_,theta_,useStorage_);
  }
  else if (type_ == "D") {
    traceSampler_ = makePtr<TraceSampler<Real>>(theta_);
    covar0_       = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_,traceSampler_);
    obj           = makePtr<Hom::D_Objective<Real>>(covar0_,theta_,useStorage_);
  }
  else if (type_ == "A") {
    bool useRandomTrace = list.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = list.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
    int nfactors = theta_->dimension();
    int size = (useRandomTrace ? nRandomTrace : nfactors);
    if (useRandomTrace) traceSampler_ = makePtr<Radamacher<Real>>(theta_,size);
    else                traceSampler_ = makePtr<TraceSampler<Real>>(theta_);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_,traceSampler_);
    lobj_   = makePtr<LinearObjective<Real>>(theta_,traceSampler_);
    obj     = makePtr<Hom::A_Objective<Real>>(covar0_,lobj_,theta_,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else if (type_ == "I" && useTrace) {
    bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
    int nfactors = theta_->dimension();
    int size = (useRandomTrace ? nRandomTrace : nfactors);
    if (useRandomTrace) traceSampler_ = makePtr<Radamacher<Real>>(theta_,size);
    else                traceSampler_ = makePtr<TraceSampler<Real>>(theta_);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    if (predFun == nullPtr) {
      covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_,traceSampler_);
    }
    else {
      factorsPV_ = makePtr<Factors<Real>>(predFun,theta_,sampler);
      covar0_ = makePtr<BilinearConstraint<Real>>(factorsPV_,cov0_,type_,traceSampler_);
    }
    lobj_   = makePtr<LinearObjective<Real>>(theta_,traceSampler_);
    obj     = makePtr<Hom::Itrace_Objective<Real>>(covar0_,lobj_,theta_,sampler,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }

  computeObjectiveScaling(obj);
  if (useScale_) obj_ = makePtr<ScaledObjective<Real>>(obj,objScale_);
  else           obj_ = obj;
}

template<typename Real>
void Factory<Real>::buildHomObjective(const Ptr<Vector<Real>> &c) {
  Ptr<Objective<Real>> obj;
  covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,c);
  lobj_  = makePtr<LinearObjective<Real>>(c);
  obj    = makePtr<Hom::C_Objective<Real>>(covar0_,lobj_,theta_,useStorage_);

  computeObjectiveScaling(obj);
  if (useScale_) obj_ = makePtr<ScaledObjective<Real>>(obj,objScale_);
  else           obj_ = obj;
}

template<typename Real>
void Factory<Real>::buildHetObjective(ParameterList &list,
                       const Ptr<SampleGenerator<Real>> &sampler,
                       const Ptr<Objective<Real>> &predFun) {
  type_ = list.sublist("OED").get("Optimality Type","C");
  Ptr<Objective<Real>> obj;
  storage_ = makePtr<SampledVector<Real>>();
  bool useTrace = list.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type_ == "I" && !useTrace) || type_ == "R") {
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_);
    if (predFun == nullPtr) {
      covar1_ = makePtr<BilinearConstraint<Real>>(factors_,cov1_,type_);
    }
    else {
      factorsPV_ = makePtr<Factors<Real>>(predFun,theta_,sampler);
      covar1_ = makePtr<BilinearConstraint<Real>>(factorsPV_,cov1_,type_);
    }
    qobj_   = makePtr<QuadraticObjective<Real>>(covar0_);
    obj     = makePtr<Het::I_Objective<Real>>(covar1_,qobj_,theta_,useStorage_);
    int dim = sampler_->getMyPoint(0).size();
    Real wt(0);
    std::vector<Real> sum(dim,0.0), mean(dim,0.0), pt(dim,0.0);
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      pt = sampler_->getMyPoint(i);
      wt = sampler_->getMyWeight(i);
      for (int j = 0; j < dim; ++j) {
        sum[j] += wt*pt[j];
      }
    }
    sampler_->sumAll(&sum[0],&mean[0],dim);
    obj->setParameter(mean);
  }
  else if (type_ == "C") {
    Ptr<Vector<Real>> c = theta_->dual().clone();
    Real cval = list.sublist("OED").sublist("C-Optimality").get("C Value",1.0);
    c->setScalar(cval);
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,c);
    covar1_ = makePtr<BilinearConstraint<Real>>(factors_,cov1_,c);
    qobj_   = makePtr<QuadraticObjective<Real>>(covar0_);
    obj     = makePtr<Het::C_Objective<Real>>(covar1_,qobj_,theta_,useStorage_);
  }
  else if (type_ == "D") {
    traceSampler_ = makePtr<TraceSampler<Real>>(theta_);
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_,traceSampler_);
    covar1_ = makePtr<BilinearConstraint<Real>>(factors_,cov1_,type_,traceSampler_);
    obj     = makePtr<Het::D_Objective<Real>>(covar0_,covar1_,theta_,useStorage_);
  }
  else if (type_ == "A") {
    bool useRandomTrace = list.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = list.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
    int nfactors = theta_->dimension();
    int size = (useRandomTrace ? nRandomTrace : nfactors);
    if (useRandomTrace) traceSampler_ = makePtr<Radamacher<Real>>(theta_,size);
    else                traceSampler_ = makePtr<TraceSampler<Real>>(theta_);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_,traceSampler_);
    covar1_ = makePtr<BilinearConstraint<Real>>(factors_,cov1_,type_,traceSampler_);
    qobj_   = makePtr<QuadraticObjective<Real>>(covar0_);
    obj     = makePtr<Het::A_Objective<Real>>(covar1_,qobj_,theta_,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else if (type_ == "I" && useTrace) {
    bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
    int nfactors = theta_->dimension();
    int size = (useRandomTrace ? nRandomTrace : nfactors);
    if (useRandomTrace) traceSampler_ = makePtr<Radamacher<Real>>(theta_,size);
    else                traceSampler_ = makePtr<TraceSampler<Real>>(theta_);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,type_,traceSampler_);
    if (predFun == nullPtr) {
      covar1_ = makePtr<BilinearConstraint<Real>>(factors_,cov1_,type_,traceSampler_);
    }
    else {
      factorsPV_ = makePtr<Factors<Real>>(predFun,theta_,sampler);
      covar1_ = makePtr<BilinearConstraint<Real>>(factorsPV_,cov1_,type_,traceSampler_);
    }
    qobj_   = makePtr<QuadraticObjective<Real>>(covar0_);
    obj     = makePtr<Het::Itrace_Objective<Real>>(covar1_,qobj_,theta_,sampler,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }

  computeObjectiveScaling(obj);
  if (useScale_) obj_ = makePtr<ScaledObjective<Real>>(obj,objScale_);
  else           obj_ = obj;
}

template<typename Real>
void Factory<Real>::buildHetObjective(const Ptr<Vector<Real>> &c) {
  Ptr<Objective<Real>> obj;
  covar0_ = makePtr<BilinearConstraint<Real>>(factors_,cov0_,c);
  covar1_ = makePtr<BilinearConstraint<Real>>(factors_,cov1_,c);
  qobj_   = makePtr<QuadraticObjective<Real>>(covar0_);
  obj     = makePtr<Het::C_Objective<Real>>(covar1_,qobj_,theta_,useStorage_);

  computeObjectiveScaling(obj);
  if (useScale_) obj_ = makePtr<ScaledObjective<Real>>(obj,objScale_);
  else           obj_ = obj;
}

template<typename Real>
void Factory<Real>::buildVector() {
  p_ = makePtr<ProbabilityVector<Real>>(
       makePtr<std::vector<Real>>(sampler_->numMySamples(),0),
       sampler_->getBatchManager());
  vec_ = p_;
}

template<typename Real>
void Factory<Real>::buildBoundConstraint() {
  //bnd_  = makePtr<Bounds<Real>>(zeros_,ones_);
  bnd_  = makePtr<Bounds<Real>>(*zeros_);
}

template<typename Real>
void Factory<Real>::buildEqualityConstraint() {
  if (!useBudget_) {
    //Real N(sampler_->numGlobalSamples());
    econ_ = makePtr<ProbabilityConstraint<Real>>(*ones_,useScale_,conScale_);
    emul_ = makePtr<DualScaledStdVector<Real>>(
              makePtr<std::vector<Real>>(1,0),
              makePtr<std::vector<Real>>(1,static_cast<Real>(1)));
              //makePtr<std::vector<Real>>(1,static_cast<Real>(1)/(N*N)));
              //makePtr<std::vector<Real>>(1,static_cast<Real>(1)/std::pow(N,2.0/3.0)));
  }
  else {
    econ_ = nullPtr;
    emul_ = nullPtr;
  }
}

template<typename Real>
void Factory<Real>::buildInequalityConstraint() {
  if (!useBudget_) {
    icon_ = nullPtr;
    imul_ = nullPtr;
    ibnd_ = nullPtr;
  }
  else {
    if (useScale_) {
      const Real one(1);
      cost_->scale(one/budget_);
      budget_ = one;
    }
    Ptr<SingletonVector<Real>> l = makePtr<SingletonVector<Real>>();
    l->zero();
    icon_ = makePtr<ScalarLinearConstraint<Real>>(cost_,budget_);
    imul_ = makePtr<SingletonVector<Real>>();
    ibnd_ = makePtr<Bounds<Real>>(*l,false);
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

template<typename Real>
void Factory<Real>::checkConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
                     std::ostream &stream) const {
  stream << std::endl;
  const Real zero(0), one(1);
  Ptr<Vector<Real>> u, z, r, p, du, dz;
  u  = theta_->clone();        u->randomize();
  z  = p_->clone();            z->randomize(zero,one);
  r  = theta_->dual().clone(); r->randomize();
  p  = theta_->clone();        p->randomize();
  du = theta_->clone();        du->randomize();
  dz = p_->clone();            dz->randomize(zero,one);

  stream << std::endl << "Check Jacobian_1 of Linear Constraint" << std::endl;
  con->checkApplyJacobian_1(*u,*z,*du,*r,true,stream);

  stream << std::endl << "Check Jacobian_2 of Linear Constraint" << std::endl;
  con->checkApplyJacobian_2(*u,*z,*dz,*r,true,stream);

  stream << std::endl << "Check Hessian_11 of Linear Constraint" << std::endl;
  con->checkApplyAdjointHessian_11(*u,*z,*p,*du,*r,true,stream);

  stream << std::endl << "Check Hessian_12 of Linear Constraint" << std::endl;
  con->checkApplyAdjointHessian_12(*u,*z,*p,*du,*dz,true,stream);

  stream << std::endl << "Check Hessian_21 of Linear Constraint" << std::endl;
  con->checkApplyAdjointHessian_21(*u,*z,*p,*dz,*r,true,stream);

  stream << std::endl << "Check Hessian_22 of Linear Constraint" << std::endl;
  con->checkApplyAdjointHessian_22(*u,*z,*p,*dz,*dz,true,stream);

  stream << std::endl << "Check Adjoint Jacobian_1 of Linear Constraint" << std::endl;
  con->checkAdjointConsistencyJacobian_1(*p,*du,*u,*z,true,stream);

  stream << std::endl << "Check Adjoint Jacobian_2 of Linear Constraint" << std::endl;
  con->checkAdjointConsistencyJacobian_2(*p,*dz,*u,*z,true,stream);

  stream << std::endl << "Check Linear Constraint Solve" << std::endl;
  con->checkSolve(*u,*z,*r,true,stream);

  stream << std::endl << "Check Inverse Jacobian_1 of Linear Constraint" << std::endl;
  con->checkInverseJacobian_1(*r,*du,*u,*z,true,stream);

  stream << std::endl << "Check Inverse Adjoint Jacobian_1 of Linear Constraint" << std::endl;
  con->checkInverseAdjointJacobian_1(*r,*p,*u,*z,true,stream);
}

template<typename Real>
void Factory<Real>::checkObjective(const Ptr<Objective_SimOpt<Real>> &obj,
                    std::ostream &stream) const {
  stream << std::endl;
  Ptr<Vector<Real>> u, z, du, dz;
  u  = theta_->clone(); u->randomize();
  z  = p_->clone();     z->randomize();
  du = theta_->clone(); du->randomize();
  dz = p_->clone();     dz->randomize();

  std::stringstream objType;
  if (isHom_) objType << "Linear";
  else        objType << "Quadratic";
  stream << std::endl << "Check Gradient_1 of " << objType.str() << " Objective" << std::endl;
  obj->checkGradient_1(*u,*z,*du,true,stream);

  stream << std::endl << "Check Gradient_2 of " << objType.str() << " Objective" << std::endl;
  obj->checkGradient_2(*u,*z,*dz,true,stream);

  stream << std::endl << "Check Hessian_11 of " << objType.str() << " Objective" << std::endl;
  obj->checkHessVec_11(*u,*z,*du,true,stream);

  stream << std::endl << "Check Hessian_12 of " << objType.str() << " Objective" << std::endl;
  obj->checkHessVec_12(*u,*z,*dz,true,stream);

  stream << std::endl << "Check Hessian_21 of " << objType.str() << " Objective" << std::endl;
  obj->checkHessVec_21(*u,*z,*du,true,stream);

  stream << std::endl << "Check Hessian_22 of " << objType.str() << " Objective" << std::endl;
  obj->checkHessVec_22(*u,*z,*dz,true,stream);
}

} // End OED Namespace
} // End ROL Namespace

#endif
