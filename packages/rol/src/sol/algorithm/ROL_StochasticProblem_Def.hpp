// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_STOCHASTICPROBLEM_DEF_HPP
#define ROL_STOCHASTICPROBLEM_DEF_HPP

namespace ROL {

template<typename Real>
StochasticProblem<Real>::StochasticProblem(const Ptr<Objective<Real>> &obj,
                                           const Ptr<Vector<Real>>    &x,
                                           const Ptr<Vector<Real>>    &g)
  : NewOptimizationProblem<Real>(obj,x,g), needRiskLessObj_(true) {}

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
  ORIGINAL_con_.insert({it->first,it->second});
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
    bool storage = list.sublist("SOL").sublist("Objective").sublist("Risk Neutral").get("Use Storage",true);
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
    mul = makePtr<DualSimulatedVector<Real>>(mvec,sampler->getBatchManager,sampler);
    res = makePtr<PrimalSimulatedVector<Real>>(rvec,sampler->getBatchManager,sampler);
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
    con = makePtr<StochasticConstraint<Real>>(con,sampler,clist);
  }
  else if ( type == "Mean Value" ) {
    conList_.insert({name,std::pair<Ptr<ParameterList>,bool>(nullPtr,true)});
    con = makePtr<MeanValueConstraint<Real>>(con,sampler);
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeConstraintStochastic: Invalid stochastic optimization type!");
  }
  if(bnd != nullPtr) NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,bnd,res,true);
  else               NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,res,true);
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
  ORIGINAL_linear_con_.insert({it->first,it->second});
  // Determine Stochastic Constraint Type
  std::string type = list.sublist("SOL").sublist(name).get("Type","Risk Neutral");
  Ptr<Constraint<Real>>      con = it->second.constraint;
  Ptr<Vector<Real>>          mul = it->second.multiplier;
  Ptr<Vector<Real>>          res = it->second.residual;
  Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
  if ( type == "Risk Neutral" ) {
    ROL_TEST_FOR_EXCEPTION(bman == nullPtr,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Risk neutral constraints need a valid BatchManager!");
    bool storage = list.sublist("SOL").sublist("Objective").sublist("Risk Neutral").get("Use Storage",true);
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
    mul = makePtr<DualSimulatedVector<Real>>(mvec,sampler->getBatchManager,sampler);
    res = makePtr<PrimalSimulatedVector<Real>>(rvec,sampler->getBatchManager,sampler);
    if (bnd != nullPtr)
      bnd = makePtr<SimulatedBoundConstraint<Real>>(sampler, bnd);
  }
  else if ( type == "Mean Value" ) {
    con = makePtr<MeanValueConstraint<Real>>(con,sampler);
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> ROL::StochasticProblem::makeLinearConstraintStochastic: Invalid stochastic optimization type!");
  }
  if(bnd != nullPtr) NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,bnd,res,true);
  else               NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,res,true);
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
    if (bnd != nullPtr) NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,bnd,res,true);
    else                NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,res,true);
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
    if (bnd != nullPtr) NewOptimizationProblem<Real>::setLinearConstraint(it->first,con,mul,res,bnd,true);
    else                NewOptimizationProblem<Real>::setLinearConstraint(it->first,con,mul,res,true);
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
void StochasticProblem<Real>::finalize(bool lumpConstraints, bool printToStream, std::ostream &outStream) {
  if (!NewOptimizationProblem<Real>::isFinalized()) {
    std::vector<Ptr<ParameterList>> conList;
    std::vector<bool> riskLessCon;
    bool flag(true), risk(!needRiskLessObj_);
    int cnt(0);
    for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
      auto it2 = conList_.find(it->first);
      if (it2==conList_.end()) {
        conList.push_back(nullPtr);
        riskLessCon.push_back(true);
      }
      else {
        conList.push_back(std::get<0>(it2->second));
        riskLessCon.push_back(std::get<1>(it2->second));
        flag = std::get<1>(it2->second);
        if (!flag) {
          dynamicPtrCast<StochasticConstraint<Real>>(it->second.constraint)->setIndex(cnt);
          risk = true;
        }
      }
      cnt++;
    }
    // Set objective function
    if (risk) {
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
      for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
        if (riskLessCon[cnt]) {
          Ptr<Constraint<Real>>      con = makePtr<RiskLessConstraint<Real>>(it->second.constraint);
          Ptr<Vector<Real>>          mul = it->second.multiplier;
          Ptr<Vector<Real>>          res = it->second.residual;
          Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
          if (bnd != nullPtr) NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,bnd,res,true);
          else                NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,res,true);
        }
        cnt++;
      }
      // Set all linear constraints to be risk less
      for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
        Ptr<Constraint<Real>>      con = makePtr<RiskLessConstraint<Real>>(it->second.constraint);
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        if (bnd != nullPtr) NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,bnd,res,true);
        else                NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,res,true);
      }
    }
    // Call default finalize
    NewOptimizationProblem<Real>::finalize(lumpConstraints,printToStream,outStream);
  }
}

}  // namespace ROL

#endif // ROL_STOCHASTICPROBLEM_DEF_HPP
