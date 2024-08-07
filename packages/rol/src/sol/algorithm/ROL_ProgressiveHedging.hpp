// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROGRESSIVEHEDGING_H
#define ROL_PROGRESSIVEHEDGING_H

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Solver.hpp"
#include "ROL_PH_Objective.hpp"
#include "ROL_PH_StatusTest.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_RiskBoundConstraint.hpp"
#include "ROL_RiskLessConstraint.hpp"

/** @ingroup algo_group
    \class ROL::ProgressiveHedging
    \brief Provides the interface to solve a stochastic program using progressive hedging.

    The progressive hedging algorithm was introduced in

    Rockafellar, R. T., and Wets, R. J-B. Scenarios and policy aggregation in
    optimization under uncertainty. Mathematics of Operations Research, 16
    (1991), 119-147.

    This algorithm solves deterministic optimization problems for each sample and then
    aggregates the optimization variables to produce an ``implementable'' solution.
    Progressive hedging has recently been applied to solve risk-averse and buffered
    probability optimization problems in

    Rockafellar, R. T., Solving stochastic programming problems with risk
    measures by progessive hedging, Set-valued and Variational Analysis,
    published online 2017.

    Rockafellar, R. T., and Uryasev, S. Minimizing buffered probability of
    exceedance by progressive hedging, Mathematical Programming B,
    published online 2018.

    This implementation can also minimize deviation, error and regret arising
    from the expectation risk quadrangle as well as the smoothed probability
    of exceedance.
    ---
*/


namespace ROL {

template <class Real>
class ProgressiveHedging {
private:
  const Ptr<Problem<Real>> input_;
  const Ptr<SampleGenerator<Real>> sampler_;
  ParameterList parlist_;
  bool usePresolve_;
  bool useInexact_;
  Real penaltyParam_;
  Real maxPen_;
  Real update_;
  int  freq_;
  Real ztol_;
  int maxit_;
  bool print_;

  bool hasStat_;
  Ptr<PH_Objective<Real>>           ph_objective_;
  Ptr<Vector<Real>>                 ph_vector_;
  Ptr<BoundConstraint<Real>>        ph_bound_;
  Ptr<Constraint<Real>>             ph_constraint_;
  Ptr<Problem<Real>>                ph_problem_;
  Ptr<Solver<Real>>                 ph_solver_;
  Ptr<PH_StatusTest<Real>>          ph_status_;
  Ptr<Vector<Real>> z_psum_, z_gsum_;
  std::vector<Ptr<Vector<Real>>> wvec_;

  void presolve(void) {
    Solver<Real> solver(input_,parlist_);
    for (int j = 0; j < sampler_->numMySamples(); ++j) {
      input_->getObjective()->setParameter(sampler_->getMyPoint(j));
      if (input_->getConstraint() != nullPtr) {
        input_->getConstraint()->setParameter(sampler_->getMyPoint(j));
      }
      solver.solve();
      z_psum_->axpy(sampler_->getMyWeight(j),*input_->getPrimalOptimizationVector());
      solver.reset();
    }
    // Aggregation
    z_gsum_->zero();
    sampler_->sumAll(*z_psum_,*z_gsum_);
  }

public:
  ProgressiveHedging(const Ptr<Problem<Real>> &input,
                     const Ptr<SampleGenerator<Real>> &sampler,
                     ParameterList &parlist)
    : input_(input), sampler_(sampler), parlist_(parlist), hasStat_(false) {
    // Get algorithmic parameters
    usePresolve_  = parlist.sublist("SOL").sublist("Progressive Hedging").get("Use Presolve",true);
    useInexact_   = parlist.sublist("SOL").sublist("Progressive Hedging").get("Use Inexact Solve",true);
    penaltyParam_ = parlist.sublist("SOL").sublist("Progressive Hedging").get("Initial Penalty Parameter",10.0);
    maxPen_       = parlist.sublist("SOL").sublist("Progressive Hedging").get("Maximum Penalty Parameter",-1.0);
    update_       = parlist.sublist("SOL").sublist("Progressive Hedging").get("Penalty Update Scale",10.0);
    freq_         = parlist.sublist("SOL").sublist("Progressive Hedging").get("Penalty Update Frequency",0);
    ztol_         = parlist.sublist("SOL").sublist("Progressive Hedging").get("Nonanticipativity Constraint Tolerance",1e-4);
    maxit_        = parlist.sublist("SOL").sublist("Progressive Hedging").get("Iteration Limit",100);
    print_        = parlist.sublist("SOL").sublist("Progressive Hedging").get("Print Subproblem Solve History",false);
    maxPen_       = (maxPen_ <= static_cast<Real>(0) ? ROL_INF<Real>() : maxPen_);
    penaltyParam_ = std::min(penaltyParam_,maxPen_);
    // Create progressive hedging vector
    ParameterList olist; olist.sublist("SOL") = parlist.sublist("SOL").sublist("Objective");
    std::string type = olist.sublist("SOL").get("Type","Risk Neutral");
    std::string prob = olist.sublist("SOL").sublist("Probability").get("Name","bPOE");
    hasStat_ = ((type=="Risk Averse") ||
                (type=="Deviation")   ||
                (type=="Probability" && prob=="bPOE"));
    Ptr<ParameterList> parlistptr = makePtrFromRef<ParameterList>(olist);
    if (hasStat_) {
      ph_vector_  = makePtr<RiskVector<Real>>(parlistptr,
                                              input_->getPrimalOptimizationVector());
    }
    else {
      ph_vector_  = input_->getPrimalOptimizationVector();
    }
    // Create progressive hedging objective function
    ph_objective_ = makePtr<PH_Objective<Real>>(input_->getObjective(),
                                                ph_vector_,
                                                penaltyParam_,
                                                olist);
    // Create progressive hedging bound constraint
    if (hasStat_) {
      ph_bound_   = makePtr<RiskBoundConstraint<Real>>(parlistptr,
                                                       input_->getBoundConstraint());
    }
    else {
      ph_bound_   = input_->getBoundConstraint();
    }
    // Create progressive hedging constraint
    ph_constraint_ = nullPtr;
    if (input_->getConstraint() != nullPtr) {
      if (hasStat_) {
        ph_constraint_ = makePtr<RiskLessConstraint<Real>>(input_->getConstraint());
      }
      else {
        ph_constraint_ = input_->getConstraint();
      }
    }
    // Build progressive hedging subproblems
    ph_problem_ = makePtr<Problem<Real>>(ph_objective_, ph_vector_);
    if (ph_bound_ != nullPtr) {
      if (ph_bound_->isActivated()) {
        ph_problem_->addBoundConstraint(ph_bound_);
      }
    }
    if (ph_constraint_ != nullPtr) {
      ph_problem_->addConstraint("PH Constraint",ph_constraint_,
                                 input_->getMultiplierVector());
    }
    // Build progressive hedging subproblem solver
    ph_solver_    = makePtr<Solver<Real>>(ph_problem_, parlist);
    // Build progressive hedging status test for inexact solves
    if (useInexact_) {
      ph_status_  = makePtr<PH_StatusTest<Real>>(parlist,
                                                 *ph_vector_);
    }
    else {
      ph_status_  = nullPtr;
    }
    // Initialize vector storage
    z_psum_       = ph_problem_->getPrimalOptimizationVector()->clone();
    z_gsum_       = ph_problem_->getPrimalOptimizationVector()->clone();
    z_gsum_->set(*ph_problem_->getPrimalOptimizationVector());
    wvec_.resize(sampler_->numMySamples());
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      wvec_[i] = z_psum_->clone(); wvec_[i]->zero();
    }
    if (usePresolve_) {
      presolve();
    }
  }

  void check(std::ostream &outStream = std::cout, const int numSamples = 1) {
    int nsamp = std::min(sampler_->numMySamples(),numSamples);
    for (int i = 0; i < nsamp; ++i) {
      ph_objective_->setParameter(sampler_->getMyPoint(i));
      ph_objective_->setData(z_gsum_,wvec_[i],penaltyParam_);
      if (ph_constraint_ != nullPtr) {
        ph_constraint_->setParameter(sampler_->getMyPoint(i));
      }
      ph_problem_->check(true,outStream);
    }
  }

  void run(std::ostream &outStream = std::cout) {
    const Real zero(0);
    std::vector<Real> vec_p(2), vec_g(2);
    Real znorm(ROL_INF<Real>()), zdotz(0);
    int iter(0), tspiter(0), flag = 1;
    bool converged = true;
    // Output
    outStream << std::scientific << std::setprecision(6);
    outStream << std::endl << "Progressive Hedging"
              << std::endl << "  "
              << std::setw(8)  << std::left << "iter"
              << std::setw(15) << std::left << "||z-Ez||"
              << std::setw(15) << std::left << "penalty"
              << std::setw(10) << std::left << "subiter"
              << std::setw(8)  << std::left << "success"
              << std::endl;
    for (iter = 0; iter < maxit_; ++iter) {
      z_psum_->zero(); vec_p[0] = zero; vec_p[1] = zero;
      ph_problem_->getPrimalOptimizationVector()->set(*z_gsum_);
      // Solve concurrent optimization problems
      for (int j = 0; j < sampler_->numMySamples(); ++j) {
        ph_objective_->setData(z_gsum_,wvec_[j],penaltyParam_);
        ph_problem_->getObjective()->setParameter(sampler_->getMyPoint(j));
        if (useInexact_) {
          ph_status_->setData(iter,z_gsum_);
        }
        if (ph_problem_->getConstraint() != nullPtr) {
          ph_problem_->getConstraint()->setParameter(sampler_->getMyPoint(j));
        }
        if (print_) {
          ph_solver_->solve(outStream,ph_status_,true);
        }
        else {
          ph_solver_->solve(ph_status_,true);
        }
        wvec_[j]->axpy(penaltyParam_,*ph_problem_->getPrimalOptimizationVector());
        vec_p[0] += sampler_->getMyWeight(j)
                   *ph_problem_->getPrimalOptimizationVector()->dot(
                   *ph_problem_->getPrimalOptimizationVector());
        vec_p[1] += static_cast<Real>(ph_solver_->getAlgorithmState()->iter);
        z_psum_->axpy(sampler_->getMyWeight(j),*ph_problem_->getPrimalOptimizationVector());
        converged = (ph_solver_->getAlgorithmState()->statusFlag == EXITSTATUS_CONVERGED
                   ||ph_solver_->getAlgorithmState()->statusFlag == EXITSTATUS_USERDEFINED
                    ? converged : false);
        ph_solver_->reset();
      }
      // Aggregation
      z_gsum_->zero(); vec_g[0] = zero; vec_g[1] = zero;
      sampler_->sumAll(*z_psum_,*z_gsum_);
      sampler_->sumAll(&vec_p[0],&vec_g[0],2);
      // Multiplier Update
      for (int j = 0; j < sampler_->numMySamples(); ++j) {
        wvec_[j]->axpy(-penaltyParam_,*z_gsum_);
      }
      zdotz  = z_gsum_->dot(*z_gsum_);
      znorm  = std::sqrt(std::abs(vec_g[0] - zdotz));
      tspiter += static_cast<int>(vec_g[1]);
      // Output
      outStream << "  "
                << std::setw(8)  << std::left << iter
                << std::setw(15) << std::left << znorm
                << std::setw(15) << std::left << penaltyParam_
                << std::setw(10) << std::left << static_cast<int>(vec_g[1])
                << std::setw(8)  << std::left << converged
                << std::endl;
      // Check termination criterion
      if (znorm <= ztol_ && converged) {
        flag = 0;
        outStream << "Converged: Nonanticipativity constraint tolerance satisfied!" << std::endl;
        break;
      }
      converged = true;
      // Update penalty parameter
      if (freq_ > 0 && iter%freq_ == 0) {
        penaltyParam_ *= update_;
      }
      penaltyParam_ = std::min(penaltyParam_,maxPen_);
    }
    if (hasStat_) {
      input_->getPrimalOptimizationVector()->set(*dynamicPtrCast<RiskVector<Real>>(z_gsum_)->getVector());
    }
    else {
      input_->getPrimalOptimizationVector()->set(*z_gsum_);
    }
    // Output reason for termination
    if (iter >= maxit_ && flag != 0) {
      outStream << "Maximum number of iterations exceeded" << std::endl;
    }
    outStream << "Total number of subproblem iterations per sample: "
              << tspiter << " / " << sampler_->numGlobalSamples()
              << " ~ " << static_cast<int>(std::ceil(static_cast<Real>(tspiter)/static_cast<Real>(sampler_->numGlobalSamples())))
              << std::endl;
  }

}; // class ProgressiveHedging

} // namespace ROL

#endif
