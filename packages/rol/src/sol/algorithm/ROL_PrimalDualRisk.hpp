// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PRIMALDUALRISK_H
#define ROL_PRIMALDUALRISK_H

#include "ROL_Solver.hpp"
#include "ROL_StochasticObjective.hpp"
#include "ROL_PD_MeanSemiDeviation.hpp"
#include "ROL_PD_MeanSemiDeviationFromTarget.hpp"
#include "ROL_PD_CVaR.hpp"
#include "ROL_PD_BPOE.hpp"
#include "ROL_PD_HMCR2.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_RiskBoundConstraint.hpp"
#include "ROL_RiskLessConstraint.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class PrimalDualRisk {
private:
  const Ptr<Problem<Real>> input_;
  const Ptr<SampleGenerator<Real>> sampler_;
  Ptr<PD_RandVarFunctional<Real>> rvf_;
  ParameterList parlist_;
  // General algorithmic parameters
  int maxit_;
  bool print_;
  Real gtolmin_;
  Real ctolmin_;
  Real ltolmin_;
  Real ltolupdate_;
  Real tolupdate0_;
  Real tolupdate1_;
  Real lalpha_;
  Real lbeta_;
  // Subproblem solver tolerances
  Real gtol_;
  Real ctol_;
  Real ltol_;
  // Penalty parameter information
  Real penaltyParam_;
  Real maxPen_;
  // Forced udpate information
  Real update_;
  int  freq_;

  Ptr<StochasticObjective<Real>> pd_objective_;
  Ptr<Vector<Real>>              pd_vector_;
  Ptr<BoundConstraint<Real>>     pd_bound_;
  Ptr<Constraint<Real>>          pd_constraint_;
  Ptr<Constraint<Real>>          pd_linear_constraint_;
  Ptr<Problem<Real>>             pd_problem_;

  int iter_, nfval_, ngrad_, ncval_;
  bool converged_;
  Real lnorm_;
  std::string name_;

public:
  PrimalDualRisk(const Ptr<Problem<Real>> &input,
                 const Ptr<SampleGenerator<Real>> &sampler,
                 ParameterList &parlist)
    : input_(input), sampler_(sampler), parlist_(parlist),
      iter_(0), converged_(true), lnorm_(ROL_INF<Real>()) {
    // Get status test information
    maxit_     = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Iteration Limit",100);
    print_     = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Print Subproblem Solve History",false);
    gtolmin_   = parlist.sublist("Status Test").get("Gradient Tolerance", 1e-8);
    ctolmin_   = parlist.sublist("Status Test").get("Constraint Tolerance", 1e-8);
    ltolmin_   = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Dual Tolerance",1e-6);
    gtolmin_   = (gtolmin_ <= static_cast<Real>(0) ?     std::sqrt(ROL_EPSILON<Real>()) : gtolmin_);
    ctolmin_   = (ctolmin_ <= static_cast<Real>(0) ?     std::sqrt(ROL_EPSILON<Real>()) : ctolmin_);
    ltolmin_   = (ltolmin_ <= static_cast<Real>(0) ? 1e2*std::sqrt(ROL_EPSILON<Real>()) : ltolmin_);
    // Get solver tolerances
    gtol_       = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Initial Gradient Tolerance", 1e-4);
    ctol_       = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Initial Constraint Tolerance", 1e-4);
    ltol_       = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Initial Dual Tolerance", 1e-2);
    ltolupdate_ = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Dual Tolerance Update Scale", 1e-1);
    tolupdate0_ = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Solver Tolerance Decrease Scale", 9e-1);
    tolupdate1_ = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Solver Tolerance Update Scale", 1e-1);
    lalpha_     = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Dual Tolerance Decrease Exponent", 1e-1);
    lbeta_      = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Dual Tolerance Update Exponent", 9e-1);
    gtol_       = std::max(gtol_, gtolmin_);
    ctol_       = std::max(ctol_, ctolmin_);
    ltol_       = std::max(ltol_, ltolmin_);
    // Get penalty parameter
    penaltyParam_ = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Initial Penalty Parameter", 10.0);
    maxPen_       = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Maximum Penalty Parameter", -1.0);
    update_       = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Penalty Update Scale", 10.0);
    maxPen_       = (maxPen_ <= static_cast<Real>(0) ? ROL_INF<Real>() : maxPen_);
    penaltyParam_ = std::min(penaltyParam_,maxPen_);
    // Get update parameters
    freq_         = parlist.sublist("SOL").sublist("Primal Dual Risk").get("Update Frequency", 0);
    // Create risk vector and risk-averse objective
    ParameterList olist; olist.sublist("SOL") = parlist.sublist("SOL").sublist("Objective");
    std::string type = olist.sublist("SOL").get<std::string>("Type");
    if (type == "Risk Averse") {
      name_ = olist.sublist("SOL").sublist("Risk Measure").get<std::string>("Name");
    }
    else if (type == "Probability") {
      name_ = olist.sublist("SOL").sublist("Probability"). get<std::string>("Name");
    }
    else {
      std::stringstream message;
      message << ">>> " << type << " is not implemented!";
      throw Exception::NotImplemented(message.str());
    }
    Ptr<ParameterList> parlistptr = makePtrFromRef<ParameterList>(olist);
    if (name_ == "CVaR") {
      parlistptr->sublist("SOL").set("Type","Risk Averse");
      parlistptr->sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
      Real alpha = olist.sublist("SOL").sublist("Risk Measure").sublist("CVaR").get("Convex Combination Parameter", 1.0);
      Real beta  = olist.sublist("SOL").sublist("Risk Measure").sublist("CVaR").get("Confidence Level",             0.9);
      rvf_ = makePtr<PD_CVaR<Real>>(alpha, beta);
    }
    else if (name_ == "Mean Plus Semi-Deviation") {
      parlistptr->sublist("SOL").set("Type","Risk Averse");
      parlistptr->sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Semi-Deviation");
      Real coeff = olist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Semi-Deviation").get("Coefficient", 0.5);
      rvf_ = makePtr<PD_MeanSemiDeviation<Real>>(coeff);
    }
    else if (name_ == "Mean Plus Semi-Deviation From Target") {
      parlistptr->sublist("SOL").set("Type","Risk Averse");
      parlistptr->sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Semi-Deviation From Target");
      Real coeff  = olist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Semi-Deviation From Target").get("Coefficient", 0.5);
      Real target = olist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Semi-Deviation From Target").get("Target", 1.0);
      rvf_ = makePtr<PD_MeanSemiDeviationFromTarget<Real>>(coeff, target);
    }
    else if (name_ == "HMCR") {
      parlistptr->sublist("SOL").set("Type","Risk Averse");
      parlistptr->sublist("SOL").sublist("Risk Measure").set("Name","HMCR");
      //Real alpha = olist.sublist("SOL").sublist("Risk Measure").sublist("HMCR").get("Convex Combination Parameter", 1.0);
      Real beta  = olist.sublist("SOL").sublist("Risk Measure").sublist("HMCR").get("Confidence Level",             0.9);
      rvf_ = makePtr<PD_HMCR2<Real>>(beta);
    }
    else if (name_ == "bPOE") {
      parlistptr->sublist("SOL").set("Type","Probability");
      parlistptr->sublist("SOL").sublist("Probability").set("Name","bPOE");
      Real thresh = olist.sublist("SOL").sublist("Probability").sublist("bPOE").get("Threshold", 1.0);
      rvf_ = makePtr<PD_BPOE<Real>>(thresh);
    }
    else {
      std::stringstream message;
      message << ">>> " << name_ << " is not implemented!";
      throw Exception::NotImplemented(message.str());
    }
    pd_vector_    = makePtr<RiskVector<Real>>(parlistptr,
                                              input_->getPrimalOptimizationVector());
    rvf_->setData(*sampler_, penaltyParam_);
    pd_objective_ = makePtr<StochasticObjective<Real>>(input_->getObjective(),
                                                       rvf_, sampler_, true);
    // Create risk bound constraint
    pd_bound_     = makePtr<RiskBoundConstraint<Real>>(parlistptr,
                                                       input_->getBoundConstraint());
    // Create riskless constraint
    pd_constraint_ = nullPtr;
    if (input_->getConstraint() != nullPtr) {
      pd_constraint_ = makePtr<RiskLessConstraint<Real>>(input_->getConstraint());
    }
    pd_linear_constraint_ = nullPtr;
    if (input_->getPolyhedralProjection() != nullPtr) {
      pd_linear_constraint_ = makePtr<RiskLessConstraint<Real>>(input_->getPolyhedralProjection()->getLinearConstraint());
    }
    // Build primal-dual subproblems
    pd_problem_ = makePtr<Problem<Real>>(pd_objective_, pd_vector_);
    if (pd_bound_->isActivated()) {
      pd_problem_->addBoundConstraint(pd_bound_);
    }
    if (pd_constraint_ != nullPtr) {
      pd_problem_->addConstraint("PD Constraint",pd_constraint_,input_->getMultiplierVector());
    }
    if (pd_linear_constraint_ != nullPtr) {
      pd_problem_->addLinearConstraint("PD Linear Constraint",pd_linear_constraint_,input_->getPolyhedralProjection()->getMultiplier());
      pd_problem_->setProjectionAlgorithm(parlist);
    }
  }

  void check(std::ostream &outStream = std::cout) {
    pd_problem_->check(true,outStream);
  }

  void run(std::ostream &outStream = std::cout) {
    const Real one(1);
    Real theta(1);
    int spiter(0);
    iter_ = 0; converged_ = true; lnorm_ = ROL_INF<Real>();
    nfval_ = 0; ncval_ = 0; ngrad_ = 0;
    // Output
    printHeader(outStream);
    Ptr<Solver<Real>> solver;
    for (iter_ = 0; iter_ < maxit_; ++iter_) {
      parlist_.sublist("Status Test").set("Gradient Tolerance",   gtol_);
      parlist_.sublist("Status Test").set("Constraint Tolerance", ctol_);
      solver = makePtr<Solver<Real>>(pd_problem_, parlist_);
      if (print_) solver->solve(outStream);
      else        solver->solve();
      converged_ = (solver->getAlgorithmState()->statusFlag == EXITSTATUS_CONVERGED
                  ||solver->getAlgorithmState()->statusFlag == EXITSTATUS_USERDEFINED
                   ? true : false);
      spiter += solver->getAlgorithmState()->iter;
      nfval_ += solver->getAlgorithmState()->nfval;
      ngrad_ += solver->getAlgorithmState()->ngrad;
      ncval_ += solver->getAlgorithmState()->ncval;
      lnorm_  = rvf_->computeDual(*sampler_);
      // Output
      print(*solver->getAlgorithmState(),outStream);
      // Check termination conditions
      if (checkStatus(*solver->getAlgorithmState(),outStream)) break;
      // Update penalty parameter and solver tolerances
      rvf_->updateDual(*sampler_);
      if (converged_) {
        if (lnorm_ > penaltyParam_*ltol_ || (freq_ > 0 && iter_%freq_ == 0)) {
          penaltyParam_  = std::min(update_*penaltyParam_, maxPen_);
          rvf_->updatePenalty(penaltyParam_);
          theta = std::min(one/penaltyParam_,one);
          ltol_ = std::max(ltolupdate_*std::pow(theta,lalpha_), ltolmin_);
          gtol_ = std::max(tolupdate0_*gtol_, gtolmin_);
          ctol_ = std::max(tolupdate0_*ctol_, ctolmin_);
        }
        else {
          theta = std::min(one/penaltyParam_,one);
          ltol_ = std::max(ltol_*std::pow(theta,lbeta_), ltolmin_);
          gtol_ = std::max(tolupdate1_*gtol_, gtolmin_);
          ctol_ = std::max(tolupdate1_*ctol_, ctolmin_);
        }
      }
    }
    input_->getPrimalOptimizationVector()->set(
      *dynamicPtrCast<RiskVector<Real>>(pd_problem_->getPrimalOptimizationVector())->getVector());
    // Output reason for termination
    if (iter_ >= maxit_) {
      outStream << "Maximum number of iterations exceeded" << std::endl;
    }
    outStream << "Primal Dual Risk required " << spiter
              << " subproblem iterations" << std::endl;
  }

private:
  void printHeader(std::ostream &outStream) const {
    std::ios_base::fmtflags flags = outStream.flags();
    outStream << std::scientific << std::setprecision(6);
    outStream << std::endl << "Primal Dual Risk Minimization using "
              << name_ << std::endl << "  "
              << std::setw(8)  << std::left << "iter"
              << std::setw(15) << std::left << "value";
    if (pd_constraint_ != nullPtr) {
      outStream << std::setw(15) << std::left << "cnorm";
    }
    outStream << std::setw(15) << std::left << "gnorm"
              << std::setw(15) << std::left << "lnorm"
              << std::setw(15) << std::left << "penalty";
    if (pd_constraint_ != nullPtr) {
      outStream << std::setw(15) << std::left << "ctol";
    }
    outStream << std::setw(15) << std::left << "gtol"
              << std::setw(15) << std::left << "ltol"
              << std::setw(10) << std::left << "nfval"
              << std::setw(10) << std::left << "ngrad";
    if (pd_constraint_ != nullPtr) {
      outStream << std::setw(10) << std::left << "ncval";
    }
    outStream << std::setw(10) << std::left << "subiter"
              << std::setw(8)  << std::left << "success"
              << std::endl;
    outStream.setf(flags);
  }

  void print(const AlgorithmState<Real> &state, std::ostream &outStream) const {
    std::ios_base::fmtflags flags = outStream.flags();
    outStream << std::scientific << std::setprecision(6);
    outStream << "  "
              << std::setw(8)  << std::left << iter_+1
              << std::setw(15) << std::left << state.value;
    if (pd_constraint_ != nullPtr) {
      outStream << std::setw(15) << std::left << state.cnorm;
    }
    outStream << std::setw(15) << std::left << state.gnorm
              << std::setw(15) << std::left << lnorm_
              << std::setw(15) << std::left << penaltyParam_;
    if (pd_constraint_ != nullPtr) {
      outStream << std::setw(15) << std::left << ctol_;
    }
    outStream << std::setw(15) << std::left << gtol_
              << std::setw(15) << std::left << ltol_
              << std::setw(10) << std::left << nfval_
              << std::setw(10) << std::left << ngrad_;
    if (pd_constraint_ != nullPtr) {
      outStream << std::setw(10) << std::left << ncval_;
    }
    outStream << std::setw(10) << std::left << state.iter
              << std::setw(8)  << std::left << converged_
              << std::endl;
    outStream.setf(flags);
  }

  bool checkStatus(const AlgorithmState<Real> &state, std::ostream &outStream) const {
    bool flag = false;
//    if (converged_ && state.iter==0 && lnorm_ < tol) {
//      outStream << "Subproblem solve converged in zero iterations"
//                << " and the difference in the multipliers was less than "
//                << tol1 << std::endl;
//      flag = true;
//    }
    if (pd_constraint_ == nullPtr) {
      if (state.gnorm < gtolmin_ && lnorm_/penaltyParam_ < ltolmin_) {
        outStream << "Solver tolerance met"
                  << " and the difference in the multipliers was less than "
                  << ltolmin_ << std::endl;
        flag = true;
      }
    }
    else {
      if (state.gnorm < gtolmin_ && state.cnorm < ctolmin_ && lnorm_/penaltyParam_ < ltolmin_) {
        outStream << "Solver tolerance met"
                  << " and the difference in the multipliers was less than "
                  << ltolmin_ << std::endl;
        flag = true;
      }
    }
    return flag;
  }

}; // class PrimalDualRisk

} // namespace ROL

#endif
