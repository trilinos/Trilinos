// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MULTIOBJECTIVEFACTORY_DEF_HPP
#define ROL_MULTIOBJECTIVEFACTORY_DEF_HPP

#include "ROL_Solver.hpp"
#include "ROL_LinearCombinationObjective.hpp"
#include "ROL_ScaledObjective.hpp"
#include "ROL_ConstraintNBI.hpp"

namespace ROL {

template<typename Real>
MultiObjectiveFactory<Real>::MultiObjectiveFactory(const Ptr<Vector<Real>> &x,
                                                   const Ptr<Vector<Real>> &g)
  : hasBounds_(false),
    hasEquality_(false), hasInequality_(false),
    hasLinearEquality_(false), hasLinearInequality_(false),
    hasProximableObjective_(false), hasProjectionList_(false),
    cnt_obj_(0), cnt_econ_(0), cnt_icon_(0), cnt_linear_econ_(0), cnt_linear_icon_(0),
    isObjInit_(false) {
  INPUT_obj_.clear();
  INPUT_nobj_  = nullPtr; 
  INPUT_xprim_ = x;
  INPUT_bnd_   = nullPtr;
  INPUT_con_.clear();
  INPUT_linear_con_.clear();
  if (g==nullPtr) INPUT_xdual_ = x->dual().clone();
  else            INPUT_xdual_ = g;
  solution_vec_.clear();
  scale_vec_.clear();
  shift_vec_.clear();
}

template<typename Real>
void MultiObjectiveFactory<Real>::addObjective(std::string name, const Ptr<Objective<Real>> &obj, bool reset) {
  if (reset) INPUT_obj_.clear();

  auto it = INPUT_obj_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_obj_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Objective names must be distinct!");

  INPUT_obj_.insert({name,obj});
  cnt_obj_++;
}

template<typename Real>
void MultiObjectiveFactory<Real>::removeObjective(std::string name) {
  auto it = INPUT_obj_.find(name);
  if (it!=INPUT_obj_.end()) {
    INPUT_obj_.erase(it);
    cnt_obj_--;
  }
}

template<typename Real>
void MultiObjectiveFactory<Real>::addBoundConstraint(const Ptr<BoundConstraint<Real>> &bnd) {
  INPUT_bnd_ = bnd;
  hasBounds_ = true;
}

template<typename Real>
void MultiObjectiveFactory<Real>::removeBoundConstraint() {
  INPUT_bnd_ = nullPtr;
  hasBounds_ = false;
}

template<typename Real>
void MultiObjectiveFactory<Real>::addProximableObjective(const Ptr<Objective<Real>> &nobj) {
  INPUT_nobj_ = nobj;
  hasProximableObjective_ = true;
}

template<typename Real>
void MultiObjectiveFactory<Real>::removeProximableObjective() {
  INPUT_nobj_ = nullPtr;
  hasProximableObjective_ = false;
}
template<typename Real>
void MultiObjectiveFactory<Real>::addConstraint( std::string                  name,
                                                 const Ptr<Constraint<Real>> &econ,
                                                 const Ptr<Vector<Real>>     &emul,
                                                 const Ptr<Vector<Real>>     &eres,
                                                 bool                         reset) {
  if (reset) INPUT_con_.clear();

  auto it = INPUT_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_con_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Constraint names must be distinct!");

  INPUT_con_.insert({name,ConstraintData<Real>(econ,emul,eres)});
  hasEquality_ = true;
  cnt_econ_++;
}

template<typename Real>
void MultiObjectiveFactory<Real>::addConstraint( std::string                       name,
                                                 const Ptr<Constraint<Real>>      &icon,
                                                 const Ptr<Vector<Real>>          &imul,
                                                 const Ptr<BoundConstraint<Real>> &ibnd,
                                                 const Ptr<Vector<Real>>          &ires,
                                                 bool                              reset) {
  if (reset) INPUT_con_.clear();

  auto it = INPUT_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_con_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Constraint names must be distinct!");

  INPUT_con_.insert({name,ConstraintData<Real>(icon,imul,ires,ibnd)});
  hasInequality_ = true;
  cnt_icon_++;
}

template<typename Real>
void MultiObjectiveFactory<Real>::removeConstraint(std::string name) {
  auto it = INPUT_con_.find(name);
  if (it!=INPUT_con_.end()) {
    if (it->second.bounds==nullPtr) cnt_econ_--;
    else                            cnt_icon_--;
    INPUT_con_.erase(it);
  }
  if (cnt_econ_==0) hasEquality_   = false;
  if (cnt_icon_==0) hasInequality_ = false;
}

template<typename Real>
void MultiObjectiveFactory<Real>::addLinearConstraint( std::string                  name,
                                                       const Ptr<Constraint<Real>> &linear_econ,
                                                       const Ptr<Vector<Real>>     &linear_emul,
                                                       const Ptr<Vector<Real>>     &linear_eres,
                                                       bool                         reset) {
  if (reset) INPUT_linear_con_.clear();

  auto it = INPUT_linear_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_linear_con_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Linear constraint names must be distinct!");

  INPUT_linear_con_.insert({name,ConstraintData<Real>(linear_econ,linear_emul,linear_eres)});
  hasLinearEquality_ = true;
  cnt_linear_econ_++;
}

template<typename Real>
void MultiObjectiveFactory<Real>::addLinearConstraint( std::string                       name,
                                                       const Ptr<Constraint<Real>>      &linear_icon,
                                                       const Ptr<Vector<Real>>          &linear_imul,
                                                       const Ptr<BoundConstraint<Real>> &linear_ibnd,
                                                       const Ptr<Vector<Real>>          &linear_ires,
                                                       bool                              reset) {
  if (reset) INPUT_linear_con_.clear();

  auto it = INPUT_linear_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_linear_con_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Linear constraint names must be distinct!");

  INPUT_linear_con_.insert({name,ConstraintData<Real>(linear_icon,linear_imul,linear_ires,linear_ibnd)});
  hasLinearInequality_ = true;
  cnt_linear_icon_++;
}

template<typename Real>
void MultiObjectiveFactory<Real>::removeLinearConstraint(std::string name) {
  auto it = INPUT_linear_con_.find(name);
  if (it!=INPUT_linear_con_.end()) {
    if (it->second.bounds==nullPtr) cnt_linear_econ_--;
    else                            cnt_linear_icon_--;
    INPUT_linear_con_.erase(it);
  }
  if (cnt_linear_econ_==0) hasLinearEquality_ = false;
  if (cnt_linear_icon_==0) hasLinearInequality_ = false;
}

template<typename Real>
void MultiObjectiveFactory<Real>::setProjectionAlgorithm(ParameterList &list) {
  ppa_list_ = list;
  hasProjectionList_ = true;
}

template<typename Real>
void MultiObjectiveFactory<Real>::addConstraintsToProblem(Ptr<Problem<Real>> &problem) {
  if (hasBounds_) problem->addBoundConstraint(INPUT_bnd_);
  if (hasProximableObjective_) problem->addProximableObjective(INPUT_nobj_);
  if (hasEquality_ || hasInequality_) {
    for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
      auto con = it->second.constraint;
      auto mul = it->second.multiplier;
      auto res = it->second.residual;
      auto bnd = it->second.bounds;
      if (con != nullPtr) {
        if (con->isActivated()) {
          if (bnd != nullPtr)
            problem->addConstraint(it->first,con,mul,bnd,res);
          else
            problem->addConstraint(it->first,con,mul,res);
        }
      }
    }
  }
  if (hasLinearEquality_ || hasLinearInequality_) {
    for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
      auto con = it->second.constraint;
      auto mul = it->second.multiplier;
      auto res = it->second.residual;
      auto bnd = it->second.bounds;
      if (con != nullPtr) {
        if (con->isActivated()) {
          if (bnd != nullPtr)
            problem->addLinearConstraint(it->first,con,mul,bnd,res);
          else
            problem->addLinearConstraint(it->first,con,mul,res);
        }
      }
    }
  }
  if (hasProjectionList_) problem->setProjectionAlgorithm(ppa_list_);
}

template<typename Real>
void MultiObjectiveFactory<Real>::computeUtopia(ParameterList &parlist, std::ostream &outStream,
                                                const Ptr<StatusTest<Real>>& status,
                                                bool combineStatus) {
  if (values_.size() != cnt_obj_) {
    solution_vec_.clear();
    scale_vec_.resize(cnt_obj_,static_cast<Real>(1));
    shift_vec_.resize(cnt_obj_,static_cast<Real>(0));
    values_.resize(cnt_obj_);
    unsigned cnt = 0u;
    for (auto it = INPUT_obj_.begin(); it != INPUT_obj_.end(); ++it) {
      // Build single-objective optimization problem
      auto problem = getScalarProblem(it->first,parlist,outStream);
      problem->finalize(false,true,outStream);
      // Solve single-objective optimization problem
      auto solver = makePtr<Solver<Real>>(problem,parlist);
      solver->solve(outStream,status,combineStatus);
      // Store Pareto data
      auto x = INPUT_xprim_->clone();
      x->set(*INPUT_xprim_);
      std::vector<Real> lam(cnt_obj_,static_cast<Real>(0));
      lam[cnt] = static_cast<Real>(1);
      values_[cnt] = evaluateObjectiveVector(*x);
      ParetoData<Real> pd(x,lam,values_[cnt],solver->getAlgorithmState()->statusFlag);
      solution_vec_.push_back(pd);
      // Update shift values
      shift_vec_[cnt] = -values_[cnt][cnt];
      cnt++;
    }
    nvalues_.assign(values_.begin(),values_.end());
  }
}

template<typename Real>
void MultiObjectiveFactory<Real>::initializeObjectives(ParameterList &parlist, std::ostream &outStream,
                                                       const Ptr<StatusTest<Real>>& status,
                                                       bool combineStatus) {
  if (!isObjInit_) {
    const bool normalize = parlist.sublist("Multi-Objective").get("Normalize Objectives",false);
    computeUtopia(parlist,outStream,status,combineStatus);
    if (normalize) {
      const Real tol = 1e-2*std::sqrt(ROL_EPSILON<Real>());
      // Compute normalization data
      for (unsigned i = 0u; i < cnt_obj_; ++i) {
        Real maxobj = values_[0][i];
        for (unsigned j = 1u; j < cnt_obj_; ++j)
          maxobj = maxobj < values_[j][i] ? values_[j][i] : maxobj;
        outStream << maxobj << "  " << values_[i][i] << "  " << maxobj-values_[i][i] << std::endl;
        if (std::abs(maxobj-values_[i][i]) > tol*maxobj) {
          scale_vec_[i] /= maxobj-values_[i][i];
          shift_vec_[i] /= maxobj-values_[i][i];
        }
        for (unsigned j = 0u; j < cnt_obj_; ++j)
          nvalues_[j][i] *= scale_vec_[i];
      }
    }
    unsigned cnt = 0u;
    obj_.clear();
    for (auto it = INPUT_obj_.begin(); it != INPUT_obj_.end(); ++it) {
      obj_.push_back(makePtr<ScaledObjective<Real>>(it->second,scale_vec_[cnt],shift_vec_[cnt]));
      cnt++;
    }
    isObjInit_ = true;
  }
}

template<typename Real>
Ptr<Problem<Real>> MultiObjectiveFactory<Real>::getScalarProblem(unsigned ind, ParameterList &parlist,
                                                                 std::ostream &outStream,
                                                                 bool initGuess, const Ptr<Vector<Real>>& x0) {
  auto it = INPUT_obj_.begin(); //+static_cast<int>(ind);
  for (unsigned i = 0; i < ind; ++i) ++it;
  ROL_TEST_FOR_EXCEPTION(it == INPUT_obj_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Objective does not exist!");
  return getScalarProblem(it->first,parlist,outStream,initGuess,x0);
}

template<typename Real>
Ptr<Problem<Real>> MultiObjectiveFactory<Real>::getScalarProblem(std::string name, ParameterList &parlist,
                                                                 std::ostream &outStream,
                                                                 bool initGuess, const Ptr<Vector<Real>>& x0) {
  auto it = INPUT_obj_.find(name);
  ROL_TEST_FOR_EXCEPTION(it == INPUT_obj_.end(),std::invalid_argument,
    ">>> ROL::MultiObjectiveFactory: Objective does not exist!");

  if (initGuess && x0 != nullPtr) INPUT_xprim_->set(*x0);
  auto problem = makePtr<Problem<Real>>(it->second,INPUT_xprim_,INPUT_xdual_);
  addConstraintsToProblem(problem);
  return problem;
}

template<typename Real>
Ptr<Problem<Real>> MultiObjectiveFactory<Real>::makeScalarProblem(const std::vector<Real> &lam, ParameterList &parlist,
                                                                  std::ostream &outStream,
                                                                  bool initGuess, const Ptr<Vector<Real>>& x0,
                                                                  const Ptr<StatusTest<Real>>& status,
                                                                  bool combineStatus) {
  Ptr<Objective<Real>>  obj;
  Ptr<Vector<Real>>     x;
  Ptr<Constraint<Real>> nbicon;
  Ptr<Vector<Real>>     nbimul;
  if (cnt_obj_ == 1u) {
    obj = INPUT_obj_.begin()->second;
    x   = INPUT_xprim_;
  }
  else {
    // Normalize the objective functions based on utopia points
    initializeObjectives(parlist,outStream,status,combineStatus);
    // Set initial guess
    x = INPUT_xprim_;
    if (initGuess && x0 != nullPtr) {
      x->set(*x0);
    }
    else {
      x->zero();
      for (unsigned i = 0u; i < cnt_obj_; ++i) x->axpy(lam[i],*solution_vec_[i].solution);
    }
    // Set up objective functions and constraints
    auto name = parlist.sublist("Multi-Objective").get("Scalarization Type","Convex Combination");
    auto type = StringToEMOType(name);
    switch (type) {
      case MOTYPE_CC: {
        obj = makePtr<LinearCombinationObjective<Real>>(lam,obj_);
        break;
      }
      case MOTYPE_NBI: {
        computeUtopia(parlist,outStream,status,combineStatus);
        obj = obj_[0];
        nbicon = makePtr<ConstraintNBI<Real>>(obj_,lam,nvalues_);
        nbimul = makePtr<StdVector<Real>>(cnt_obj_-1u);
        break;
      }
      case MOTYPE_LAST:
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>> ROL::MultiObjectiveFactory: Invalid scalar problem type!");
        break;
    }
  }
  auto problem = makePtr<Problem<Real>>(obj,x,INPUT_xdual_);
  addConstraintsToProblem(problem);
  if (nbicon != nullPtr && nbimul != nullPtr)
    problem->addConstraint("NBI",nbicon,nbimul);
  return problem;
}

}  // namespace ROL

#endif // ROL_MULTIOBJECTIVEFACTORY_DEF_HPP
