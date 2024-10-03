// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STOCHASTICPROBLEM_DEF_HPP
#define ROL_STOCHASTICPROBLEM_DEF_HPP

namespace ROL {

template<typename Real>
StochasticProblem<Real>::StochasticProblem(const Ptr<Objective<Real>> &obj,
                                           const Ptr<Vector<Real>>    &x,
                                           const Ptr<Vector<Real>>    &g)
  : Problem<Real>(obj,x,g), needRiskLessObj_(true) {}

template<typename Real>
void StochasticProblem<Real>::makeObjectiveStochastic(ParameterList                    &list,
                                                      const Ptr<SampleGenerator<Real>> &fsampler,
                                                      const Ptr<SampleGenerator<Real>> &gsampler,
                                                      const Ptr<SampleGenerator<Real>> &hsampler) {
  // Throw an exception if problem has been finalized
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeObjectiveStochastic: Cannot set stochastic objective after problem has been finalized!");
  // Throw an exception if the value sampler is null
  ROL_TEST_FOR_EXCEPTION(fsampler == nullPtr,std::invalid_argument,
    ">>> ROL::StochasticProblem::makeObjectiveStochastic: Objective function value sampler is null!");
  // Store original objective function for reuse later
  ORIGINAL_obj_ = INPUT_obj_;
  // Check samplers
  Ptr<SampleGenerator<Real>> _gsampler, _hsampler;
  _gsampler = (gsampler == nullPtr ?  fsampler : gsampler);
  _hsampler = (hsampler == nullPtr ? _gsampler : hsampler);
  // Determine Stochastic Objective Type
  std::string type = list.sublist("SOL").sublist("Objective").get("Type","Risk Neutral");
  if ( type == "Risk Neutral" ) {
    needRiskLessObj_ = true;
    objList_ = nullPtr;
    bool storage = list.sublist("SOL").sublist("Objective").sublist("Risk Neutral").get("Use Storage",true);
    INPUT_obj_ = makePtr<RiskNeutralObjective<Real>>(ORIGINAL_obj_,fsampler,_gsampler,_hsampler,storage);
  }
  else if ( type == "Risk Averse" || type == "Deviation" || type == "Error" ||
            type == "Regret"      || type == "Probability" ) {
    needRiskLessObj_ = false;
    objList_ = makePtr<ParameterList>();
    objList_->sublist("SOL") = list.sublist("SOL").sublist("Objective");
    INPUT_obj_ = makePtr<StochasticObjective<Real>>(ORIGINAL_obj_,*objList_,fsampler,_gsampler,_hsampler);
  }
  else if ( type == "Mean Value" ) {
    needRiskLessObj_ = true;
    objList_ = nullPtr;
    INPUT_obj_ = makePtr<MeanValueObjective<Real>>(ORIGINAL_obj_,fsampler);
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeObjectiveStochastic: Invalid stochastic optimization type!");
  }
}

template<typename Real>
void StochasticProblem<Real>::makeObjectiveStochastic(const Ptr<RandVarFunctional<Real>> &rvf,
                                                      ParameterList                      &list,
                                                      const Ptr<SampleGenerator<Real>>   &fsampler,
                                                      const Ptr<SampleGenerator<Real>>   &gsampler,
                                                      const Ptr<SampleGenerator<Real>>   &hsampler) {
  // Throw an exception if problem has been finalized
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeObjectiveStochastic: Cannot set stochastic objective after problem has been finalized!");
  // Throw an exception if the value sampler is null
  ROL_TEST_FOR_EXCEPTION(fsampler == nullPtr,std::invalid_argument,
    ">>> ROL::StochasticProblem::makeObjectiveStochastic: Objective function value sampler is null!");
  // Throw an exception if the value sampler is null
  ROL_TEST_FOR_EXCEPTION(rvf == nullPtr,std::invalid_argument,
    ">>> ROL::StochasticProblem::makeObjectiveStochastic: Risk measure is null!");
  // Store original objective function for reuse later
  ORIGINAL_obj_ = INPUT_obj_;
  // Check samplers
  Ptr<SampleGenerator<Real>> _gsampler, _hsampler;
  _gsampler = (gsampler == nullPtr ?  fsampler : gsampler);
  _hsampler = (hsampler == nullPtr ? _gsampler : hsampler);
  // Determine Stochastic Objective Type
  needRiskLessObj_ = false;
  objList_  = makePtr<ParameterList>();
  *objList_ = list;
  //objList_->sublist("SOL") = list.sublist("SOL").sublist("Objective");
  INPUT_obj_ = makePtr<StochasticObjective<Real>>(ORIGINAL_obj_,rvf,fsampler,_gsampler,_hsampler);
}

template<typename Real>
void StochasticProblem<Real>::makeConstraintStochastic(std::string                       name,
                                                       ParameterList                    &list,
                                                       const Ptr<SampleGenerator<Real>> &sampler,
                                                       const Ptr<BatchManager<Real>>    &bman) {
  // Throw an exception if problem has been finalized
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeConstraintStochastic: Cannot set stochastic constraint after problem has been finalized!");
  // Throw an exception if the value sampler is null
  ROL_TEST_FOR_EXCEPTION(sampler == nullPtr,std::invalid_argument,
    ">>> ROL::StochasticProblem::makeConstraintStochastic: Constraint sampler is null!");
  // Store original constraint for reuse later
  auto it = INPUT_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it == INPUT_con_.end(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeConstraintStochastic: Constraint does not exist!");
  ROL_TEST_FOR_EXCEPTION(ORIGINAL_con_.find(name) != ORIGINAL_con_.end(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeConstraintStochastic: Constraint already set!");
  ORIGINAL_con_.insert({name,it->second});
  // Determine Stochastic Constraint Type
  std::string type = list.sublist("SOL").sublist(name).get("Type","Risk Neutral");
  Ptr<Constraint<Real>>      con = it->second.constraint;
  Ptr<Vector<Real>>          mul = it->second.multiplier;
  Ptr<Vector<Real>>          res = it->second.residual;
  Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
  if ( type == "Risk Neutral" ) {
    ROL_TEST_FOR_EXCEPTION(bman == nullPtr,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeConstraintStochastic: Risk neutral constraints need a valid BatchManager!");
    conList_.insert({name,std::pair<Ptr<ParameterList>,bool>(nullPtr,true)});
    con = makePtr<RiskNeutralConstraint<Real>>(it->second.constraint,sampler,bman);
  }
  else if ( type == "Almost Sure" ) {
    conList_.insert({name,std::pair<Ptr<ParameterList>,bool>(nullPtr,true)});
    int nsamp = sampler->numMySamples();
    con = makePtr<AlmostSureConstraint<Real>>(sampler,it->second.constraint);
    std::vector<Ptr<Vector<Real>>> mvec(nsamp,nullPtr), rvec(nsamp,nullPtr);
    for (int j = 0; j < nsamp; ++j) {
      mvec[j] = mul->clone(); mvec[j]->set(*mul);
      rvec[j] = res->clone(); rvec[j]->set(*res);
    }
    mul = makePtr<DualSimulatedVector<Real>>(mvec,sampler->getBatchManager(),sampler);
    res = makePtr<PrimalSimulatedVector<Real>>(rvec,sampler->getBatchManager(),sampler);
    if (bnd != nullPtr)
      bnd = makePtr<SimulatedBoundConstraint<Real>>(sampler, bnd);
  }
  else if ( type == "Risk Averse" || type == "Deviation" || type == "Error" ||
            type == "Regret"      || type == "Probability" ) {
    ROL_TEST_FOR_EXCEPTION(bnd == nullPtr,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeConstraintStochastic: Stochastic constraints must be inequalities!");
    Ptr<ParameterList> clist = makePtr<ParameterList>();
    clist->sublist("SOL") = list.sublist("SOL").sublist(name);
    conList_.insert({name,std::pair<Ptr<ParameterList>,bool>(clist,false)});
    con = makePtr<StochasticConstraint<Real>>(it->second.constraint,sampler,*clist);
  }
  else if ( type == "Mean Value" ) {
    conList_.insert({name,std::pair<Ptr<ParameterList>,bool>(nullPtr,true)});
    con = makePtr<MeanValueConstraint<Real>>(it->second.constraint,sampler);
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeConstraintStochastic: Invalid stochastic optimization type!");
  }
  Problem<Real>::removeConstraint(name);
  if(bnd != nullPtr) Problem<Real>::addConstraint(name,con,mul,bnd,res);
  else               Problem<Real>::addConstraint(name,con,mul,res);
}

template<typename Real>
void StochasticProblem<Real>::makeLinearConstraintStochastic(std::string                       name,
                                                             ParameterList                    &list,
                                                             const Ptr<SampleGenerator<Real>> &sampler,
                                                             const Ptr<BatchManager<Real>>    &bman) {
  // Throw an exception if problem has been finalized
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Cannot set stochastic constraint after problem has been finalized!");
  // Throw an exception if the value sampler is null                                
  ROL_TEST_FOR_EXCEPTION(sampler == nullPtr,std::invalid_argument,
    ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Constraint sampler is null!");
  // Store original constraint for reuse later
  auto it = INPUT_linear_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it == INPUT_linear_con_.end(),std::invalid_argument,
    ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Constraint does not exist!");
  ROL_TEST_FOR_EXCEPTION(ORIGINAL_linear_con_.find(name) != ORIGINAL_linear_con_.end(),std::invalid_argument,
      ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Constraint already set!");
  ORIGINAL_linear_con_.insert({name,it->second});
  // Determine Stochastic Constraint Type
  std::string type = list.sublist("SOL").sublist(name).get("Type","Risk Neutral");
  Ptr<Constraint<Real>>      con = it->second.constraint;
  Ptr<Vector<Real>>          mul = it->second.multiplier;
  Ptr<Vector<Real>>          res = it->second.residual;
  Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
  if ( type == "Risk Neutral" ) {
    ROL_TEST_FOR_EXCEPTION(bman == nullPtr,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Risk neutral constraints need a valid BatchManager!");
    con = makePtr<RiskNeutralConstraint<Real>>(it->second.constraint,sampler,bman);
  }
  else if ( type == "Almost Sure" ) {
    int nsamp = sampler->numMySamples();
    con = makePtr<AlmostSureConstraint<Real>>(sampler,it->second.constraint);
    std::vector<Ptr<Vector<Real>>> mvec(nsamp,nullPtr), rvec(nsamp,nullPtr);
    for (int j = 0; j < nsamp; ++j) {
      mvec[j] = mul->clone(); mvec[j]->set(*mul);
      rvec[j] = res->clone(); rvec[j]->set(*res);
    }
    mul = makePtr<DualSimulatedVector<Real>>(mvec,sampler->getBatchManager(),sampler);
    res = makePtr<PrimalSimulatedVector<Real>>(rvec,sampler->getBatchManager(),sampler);
    if (bnd != nullPtr)
      bnd = makePtr<SimulatedBoundConstraint<Real>>(sampler, bnd);
  }
  else if ( type == "Mean Value" ) {
    con = makePtr<MeanValueConstraint<Real>>(it->second.constraint,sampler);
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Invalid stochastic optimization type!");
  }
  Problem<Real>::removeLinearConstraint(name);
  if(bnd != nullPtr) Problem<Real>::addLinearConstraint(name,con,mul,bnd,res);
  else               Problem<Real>::addLinearConstraint(name,con,mul,res);
}

template<typename Real>
void StochasticProblem<Real>::resetStochasticObjective(void) {
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::resetStochasticObjective: Cannot reset stochastic objective after problem has been finalized!");
  if (ORIGINAL_obj_ != nullPtr) {
    INPUT_obj_       = ORIGINAL_obj_;
    needRiskLessObj_ = true;
    objList_         = nullPtr;
  }
  ORIGINAL_obj_ = nullPtr;
}

template<typename Real>
void StochasticProblem<Real>::resetStochasticConstraint(std::string name) {
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::resetStochasticConstraint: Cannot reset stochastic constraint after problem has been finalized!");
  auto it = ORIGINAL_con_.find(name);
  if (it != ORIGINAL_con_.end()) {
    Ptr<Constraint<Real>>      con = it->second.constraint;
    Ptr<Vector<Real>>          mul = it->second.multiplier;
    Ptr<Vector<Real>>          res = it->second.residual;
    Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
    Problem<Real>::removeConstraint(name);
    if (bnd != nullPtr) Problem<Real>::addConstraint(name,con,mul,bnd,res);
    else                Problem<Real>::addConstraint(name,con,mul,res);
    conList_.erase(conList_.find(name));
    ORIGINAL_con_.erase(it);
  }
}

template<typename Real>
void StochasticProblem<Real>::resetStochasticLinearConstraint(std::string name) {
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::resetStochasticLinearConstraint: Cannot reset stochastic constraint after problem has been finalized!");
  auto it = ORIGINAL_linear_con_.find(name);
  if (it != ORIGINAL_linear_con_.end()) {
    Ptr<Constraint<Real>>      con = it->second.constraint;
    Ptr<Vector<Real>>          mul = it->second.multiplier;
    Ptr<Vector<Real>>          res = it->second.residual;
    Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
    Problem<Real>::removeLinearConstraint(name);
    if (bnd != nullPtr) Problem<Real>::addLinearConstraint(name,con,mul,bnd,res);
    else                Problem<Real>::addLinearConstraint(name,con,mul,res);
    ORIGINAL_linear_con_.erase(it);
  }
}

template<typename Real>
void StochasticProblem<Real>::resetStochastic(void) {
  ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::reset: Cannot reset stochastic problem after problem has been finalized!");
  // Reset objective
  resetStochasticObjective();
  // Reset general constraints
  std::vector<std::string> names;
  for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
    names.push_back(it->first);
  }
  for (auto it = names.begin(); it != names.end(); ++it) {
    resetStochasticConstraint(*it);
  }
  // Reset linear constraints
  names.clear();
  for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
    names.push_back(it->first);
  }
  for (auto it = names.begin(); it != names.end(); ++it) {
    resetStochasticLinearConstraint(*it);
  }
  // Reset primal optimization variables
  if (ORIGINAL_xprim_ != nullPtr) {
    INPUT_xprim_ = ORIGINAL_xprim_;
    ORIGINAL_xprim_ = nullPtr;
  }
  // Reset dual optimization variables
  if (ORIGINAL_xdual_ != nullPtr) {
    INPUT_xdual_ = ORIGINAL_xdual_;
    ORIGINAL_xdual_ = nullPtr;
  }
  // Reset bound constraint
  if (ORIGINAL_bnd_ != nullPtr) {
    INPUT_bnd_ = ORIGINAL_bnd_;
    ORIGINAL_bnd_ = nullPtr;
  }
}

template<typename Real>
std::vector<Real> StochasticProblem<Real>::getObjectiveStatistic(void) const {
  ROL_TEST_FOR_EXCEPTION(!isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::getObjectiveStatistic: Cannot get statistic if problem has not been finalized!");
  try {
    Ptr<std::vector<Real>> stat
      = dynamicPtrCast<RiskVector<Real>>(INPUT_xprim_)->getStatistic();
    if (stat != nullPtr) return *stat;
    else                 return std::vector<Real>();
  }
  catch (std::exception &e) {
    return std::vector<Real>();
  }
}

template<typename Real>
std::vector<Real> StochasticProblem<Real>::getConstraintStatistic(std::string name) const {
  ROL_TEST_FOR_EXCEPTION(!isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::getConstraintStatistic: Cannot get statistic if problem has not been finalized!");
  auto it = statMap_.find(name);
  ROL_TEST_FOR_EXCEPTION(it==statMap_.end(),std::invalid_argument,
    ">>> ROL::StochasticProblem::getConstraintStatistic: Constraint does not exist!");
  try {
    Ptr<std::vector<Real>> stat
      = dynamicPtrCast<RiskVector<Real>>(INPUT_xprim_)->getStatistic(1,it->second);
    if (stat != nullPtr) return *stat;
    else                 return std::vector<Real>();
  }
  catch (std::exception &e) {
    return std::vector<Real>();
  }
}

template<typename Real>
Real StochasticProblem<Real>::getSolutionStatistic(int comp, std::string name) const {
  ROL_TEST_FOR_EXCEPTION(!isFinalized(),std::invalid_argument,
    ">>> ROL::StochasticProblem::getConstraintStatistic: Cannot get statistic if problem has not been finalized!");
  ROL_TEST_FOR_EXCEPTION(comp>1||comp<0,std::invalid_argument,
    ">>> ROL::StochasticProblem::getSolutionStatistic: Component must be either 0 or 1!");
  Real val(0);
  if (comp == 0) {
    try {
      val = dynamicPtrCast<StochasticObjective<Real>>(INPUT_obj_)->computeStatistic(*INPUT_xprim_);
    }
    catch (std::exception &e) {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> ROL::StochasticProblem::getSolutionStatistic: Objective does not have a computeStatistic function!");
    }
  }
  else {
    auto it = statMap_.find(name);
    ROL_TEST_FOR_EXCEPTION(it==statMap_.end(),std::invalid_argument,
      ">>> ROL::StochasticProblem::getSolutionStatistic: Constraint does not exist!");
    try {
      auto it2 = INPUT_con_.find(name);
      val = dynamicPtrCast<StochasticConstraint<Real>>(it2->second.constraint)->computeStatistic(*INPUT_xprim_);
    }
    catch (std::exception &e) {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> ROL::StochasticProblem::getSolutionStatistic: Constraint does not have a computeStatistic function!");
    }
  }
  return val;
}

template<typename Real>
void StochasticProblem<Real>::finalize(bool lumpConstraints, bool printToStream, std::ostream &outStream) {
  if (!Problem<Real>::isFinalized()) {
    std::vector<Ptr<ParameterList>> conList;
    bool flag(true);
    risk_ = !needRiskLessObj_;
    size_t cnt(0);
    needRiskLessCon_.clear();
    statMap_.clear();
    for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
      auto it2 = conList_.find(it->first);
      if (it2==conList_.end()) {
        conList.push_back(nullPtr);
        needRiskLessCon_.push_back(true);
      }
      else {
        conList.push_back(std::get<0>(it2->second));
        needRiskLessCon_.push_back(std::get<1>(it2->second));
        flag = std::get<1>(it2->second);
        if (!flag) {
          dynamicPtrCast<StochasticConstraint<Real>>(it->second.constraint)->setIndex(cnt);
          risk_ = true;
        }
      }
      statMap_.insert({it->first,cnt});
      cnt++;
    }
    // Set objective function
    if (risk_) {
      if (needRiskLessObj_) {
        Ptr<Objective<Real>> obj = INPUT_obj_;
        INPUT_obj_ = makePtr<RiskLessObjective<Real>>(obj);
      }
      // Set risk vector
      ORIGINAL_xprim_ = INPUT_xprim_;
      INPUT_xprim_ = makePtr<RiskVector<Real>>(objList_,conList,ORIGINAL_xprim_);
      ORIGINAL_xdual_ = INPUT_xdual_;
      INPUT_xdual_ = makePtr<RiskVector<Real>>(objList_,conList,ORIGINAL_xdual_);
      if (objList_ != nullPtr) {
        Real statObj = objList_->sublist("SOL").get("Initial Statistic",1.0);
        dynamicPtrCast<RiskVector<Real>>(INPUT_xprim_)->setStatistic(statObj,0);
      }
      for (size_t i = 0; i < conList.size(); ++i) {
        if (conList[i] != nullPtr) {
          Real statCon = conList[i]->sublist("SOL").get("Initial Statistic",1.0);
          dynamicPtrCast<RiskVector<Real>>(INPUT_xprim_)->setStatistic(statCon,1,i);
        }
      }
      // Set risk bound constraint
      if (INPUT_bnd_ != nullPtr) {
        ORIGINAL_bnd_ = INPUT_bnd_;
        INPUT_bnd_ = makePtr<RiskBoundConstraint<Real>>(objList_,conList,ORIGINAL_bnd_);
      }
      // Set appropriate general constraints to be risk less
      cnt = 0;
      std::unordered_map<std::string,ConstraintData<Real>> riskless_con;
      for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
        if (needRiskLessCon_[cnt]) {
          Ptr<Constraint<Real>>      con = makePtr<RiskLessConstraint<Real>>(it->second.constraint);
          Ptr<Vector<Real>>          mul = it->second.multiplier;
          Ptr<Vector<Real>>          res = it->second.residual;
          Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
          riskless_con.insert({it->first,ConstraintData<Real>(con,mul,res,bnd)});
          if (ORIGINAL_con_.count(it->first) == size_t(0))
            ORIGINAL_con_.insert({it->first,ConstraintData<Real>(it->second.constraint,mul,res,bnd)});
        }
        cnt++;
      }
      for (auto it = riskless_con.begin(); it != riskless_con.end(); ++it) {
        Ptr<Constraint<Real>>      con = it->second.constraint;
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        Problem<Real>::removeConstraint(it->first);
        if (bnd != nullPtr) Problem<Real>::addConstraint(it->first,con,mul,bnd,res);
        else                Problem<Real>::addConstraint(it->first,con,mul,res);
      }
      // Set all linear constraints to be risk less
      riskless_con.clear();
      for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
        Ptr<Constraint<Real>>      con = makePtr<RiskLessConstraint<Real>>(it->second.constraint);
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        riskless_con.insert({it->first,ConstraintData<Real>(con,mul,res,bnd)});
        if (ORIGINAL_linear_con_.count(it->first) == size_t(0))
          ORIGINAL_linear_con_.insert({it->first,ConstraintData<Real>(it->second.constraint,mul,res,bnd)});
      }
      for (auto it = riskless_con.begin(); it != riskless_con.end(); ++it) {
        Ptr<Constraint<Real>>      con = it->second.constraint;
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        Problem<Real>::removeLinearConstraint(it->first);
        if (bnd != nullPtr) Problem<Real>::addLinearConstraint(it->first,con,mul,bnd,res);
        else                Problem<Real>::addLinearConstraint(it->first,con,mul,res);
      }
    }
    // Call default finalize
    Problem<Real>::finalize(lumpConstraints,printToStream,outStream);
  }
}

template<typename Real>
void StochasticProblem<Real>::edit(void) {
  Problem<Real>::edit();

  if (risk_) {
    if (needRiskLessObj_ && ORIGINAL_obj_ != nullPtr) {
      INPUT_obj_ = ORIGINAL_obj_;
      ORIGINAL_obj_ = nullPtr;
    }
    if (ORIGINAL_xprim_ != nullPtr) INPUT_xprim_ = ORIGINAL_xprim_;
    if (ORIGINAL_xdual_ != nullPtr) INPUT_xdual_ = ORIGINAL_xdual_;
    if (ORIGINAL_bnd_   != nullPtr) INPUT_bnd_   = ORIGINAL_bnd_;
    ORIGINAL_xprim_ = nullPtr;
    ORIGINAL_xdual_ = nullPtr;
    ORIGINAL_bnd_   = nullPtr;
    size_t cnt = 0;

    std::unordered_map<std::string,ConstraintData<Real>> riskless_con;
    for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
      if (needRiskLessCon_[cnt]) {
        Ptr<Constraint<Real>>      con = it->second.constraint;
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        riskless_con.insert({it->first,ConstraintData<Real>(con,mul,res,bnd)});
      }
      cnt++;
    }
    for (auto it = riskless_con.begin(); it != riskless_con.end(); ++it) {
      auto it2 = ORIGINAL_con_.find(it->first);
      if (it2 != ORIGINAL_con_.end()) {
        Ptr<Constraint<Real>>      con = it2->second.constraint;
        Ptr<Vector<Real>>          mul = it2->second.multiplier;
        Ptr<Vector<Real>>          res = it2->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it2->second.bounds;
        Problem<Real>::removeConstraint(it2->first);
        if (bnd != nullPtr) Problem<Real>::addConstraint(it2->first,con,mul,bnd,res);
        else                Problem<Real>::addConstraint(it2->first,con,mul,res);
        ORIGINAL_con_.erase(it2);
      }
    }
    // Set all linear constraints to be risk less
    riskless_con.clear();
    for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
      Ptr<Constraint<Real>>      con = it->second.constraint;
      Ptr<Vector<Real>>          mul = it->second.multiplier;
      Ptr<Vector<Real>>          res = it->second.residual;
      Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
      riskless_con.insert({it->first,ConstraintData<Real>(con,mul,res,bnd)});
    }
    for (auto it = riskless_con.begin(); it != riskless_con.end(); ++it) {
      auto it2 = ORIGINAL_linear_con_.find(it->first);
      if (it2 != ORIGINAL_linear_con_.end()) {
        Ptr<Constraint<Real>>      con = it2->second.constraint;
        Ptr<Vector<Real>>          mul = it2->second.multiplier;
        Ptr<Vector<Real>>          res = it2->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it2->second.bounds;
        Problem<Real>::removeLinearConstraint(it2->first);
        if (bnd != nullPtr) Problem<Real>::addLinearConstraint(it2->first,con,mul,bnd,res);
        else                Problem<Real>::addLinearConstraint(it2->first,con,mul,res);
        ORIGINAL_linear_con_.erase(it2);
      }
    }
  }
  risk_ = false;
  needRiskLessCon_.clear();
  statMap_.clear();
}

}  // namespace ROL

#endif // ROL_STOCHASTICPROBLEM_DEF_HPP
